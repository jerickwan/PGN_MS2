'''
Generator Classes
env = chem
'''
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools
from pathlib import Path

import itertools
import functools
import operator
import re

import pickle
import pandas as pd
import numpy as np

import copy
import yaml

from Molecules import AMINO_ACID_DB, GLYCAN_DB,\
    ALL_AMINO_ACID_LST, ALL_GLC_LST, ALL_MUR_LST
from Molecules import Molecule, AminoAcid, Peptide, Glycan, Peptidoglycan
from Molecules import NULL_AA
from Molecules import run_two_pdt_rxn, run_polymerisation_rxn, create_mol_cache_key
from Exceptions import NoReactionError, TooManyPdtsError, InputError
from Illustrator import Illustrator
from Common import BOND_CHAR, BRANCH_CHAR,\
    BATCH_SIZE, MEMORY, TIME_STRING, OUTPUT_ADDUCTS
from Common import Counter
from Common import make_dir, flatten

# %% Globals

Lac = Molecule("OC(C)C(O)=O")
gly_Lac = Glycan("OC(C)C(O)=O", "Lac", None, None)  # Used for lactoyl peptides
GEN_MAX_PEPTIDE_LENGTH = 8
GEN_MIN_PRODUCT_MASS = 250
GEN_WARNING_THRESHOLD = 20000  # Ask for confirmation if no. > threshold
ILL_FONTSIZE = 12
ILL_SIZE = "A4"


# %% Cached

@MEMORY.cache(ignore=["generator", "glycans", "peptides"])
def cache_PGN_generation(generator, total, cache_key, glycans, peptides):
    '''Returns list of PGN generated from given glycans and peptides'''
    to_add = []
    generator.counter.set_counter(total, "monomers")
    combs = itertools.product(glycans, peptides)
    for p, g in combs:
        to_add.append(p & g)
        generator.counter.update_counter()
    return to_add


@MEMORY.cache(ignore=["generator", "valid_P1", "valid_P2"])
def cache_polymer_generation(generator, total, name,
                               P1_bond, P2_bond,
                               cache_key,
                               valid_P1, valid_P2):
    '''Returns list of polymerised PGN generated from given PGN and bonding criteria'''
    to_add = []
    generator.counter.set_counter(total, name)
    combs = itertools.product(valid_P1, valid_P2)
    for P1, P2 in combs:
        pdt = run_polymerisation_rxn(name, P1_bond, P2_bond, P1, P2)
        if pdt not in to_add:
            to_add.append(pdt)
        else:
            print(f"Duplicates of {pdt.name} exist!")
        generator.counter.update_counter()
    return to_add

# %% Generator


class Generator():

    # %%% Main

    valid_AAs = [*ALL_AMINO_ACID_LST, None]
    valid_glycans = {"Glc": [*ALL_GLC_LST, None],
                     "Mur": [*ALL_MUR_LST, None]}

    def __init__(self,
                 name=TIME_STRING,
                 modifications=None,
                 glycan_lmin=0,
                 glycan_lmax=0,
                 peptide_lmin=0,
                 peptide_lmax=0,
                 num_polymerisations=0,
                 diffcalc_units_min=99,
                 diffcalc_maxdiff_per_unit=None,
                 warning_threshold=GEN_WARNING_THRESHOLD,
                 output_adducts=OUTPUT_ADDUCTS
                 ):
        """


        Parameters
        ----------
        name : string, optional
            The name of Generator (used in filenames). The default is TIME_STRING.
        modifications : dict, optional
            Modifications applied to Peptidoglycan.
            The default is None (no mods).
        glycan_lmin : int, optional
            Min glycan length. The default is 0.
        glycan_lmax : int, optional
            Max glycan length. The default is 0.
        peptide_lmin : int, optional
            Min peptide length. The default is 0.
        peptide_lmax : int, optional
            Max peptide length. The default is 0.
        num_polymerisations : int, optional
            Max polymerisations per molecule. The default is 0.
        diffcalc_units_min : int, optional
            Minimum no. of PGN units to apply diffcalc for polymerisations.
            The default is 99 (not used).
        diffcalc_maxdiff_per_unit : list, optional
            No. of differences allowed per PGN unit for each. First integer in
            list is used for 2-mers, second 3-mers,... etc.
            The default is None (not used).
        output_adducts : list, optional
            List of output adducts. The default is OUTPUT_ADDUCTS.
        warning_threshold : int, optional
            Generation of n-mers above threshold will prompt for confirmation
            first. The default is GEN_WARNING_THRESHOLD.

        Returns
        -------
        None.

        """
        # Naming and Directory
        self.name = None
        self.dir = None
        self.refresh_name(name)

        # PGN Modifications
        if modifications is None:
            self.modifications = {
                mod: False for mod in Peptidoglycan.modifications}
        else:
            self.modifications = modifications

        # PGN Glycan
        self.glycan_lmin = glycan_lmin
        self.glycan_lmax = glycan_lmax
        self.glycan_units = {}
        self.glycan_combs = None

        # PGN Peptide
        self.peptide_lmin = peptide_lmin
        self.peptide_lmax = peptide_lmax
        self.peptide_residues = {}
        self.peptide_combs = None

        # PGN Bridge Peptides
        self.bridge_peptides = {}  # idx : grp: [peptides]
        self.bridge_peptides_processed = None

        # PGN Polymers
        self.num_polymerisations = num_polymerisations  # number of polymerisations
        self.polymerisation_types = []  # criteria for polymerisations
        self.polymerisation_types_processed = None

        # Diffcalc
        self.diffcalc_units_min = diffcalc_units_min
        if diffcalc_maxdiff_per_unit is None:
            self.diffcalc_maxdiff_per_unit = [99]
        else:
            self.diffcalc_maxdiff_per_unit = diffcalc_maxdiff_per_unit
        self.diffcalc = lambda x: 0  # diff calculator
        self.diffcalc_params = {}

        # Misc
        self.warning_threshold = warning_threshold
        # use dummy AAs as temporary bridge peptides to reduce complexity
        self.dummy_AAs = {
            i: AMINO_ACID_DB[f"{i}-Lat"] for i in range(1, 6)}
        self.counter = Counter()
        self.master_dataframe = None
        self.refresh_parameters()
        self.PGN_dict = {}
        self.reference_dict = {}
        self.output_adducts = output_adducts

        print(f"\tOutput adducts:\t{self.output_adducts}")
        print(f"\tMin. mass:\t{GEN_MIN_PRODUCT_MASS}")
        print(f"\tMax. peptide length:\t{GEN_MAX_PEPTIDE_LENGTH}")

    def __len__(self):
        """
        Returns total number of PGN generated.

        Returns
        -------
        int
            Total number of PGN generated

        """
        return sum(len(self.PGN_dict[k]) for k in self.PGN_dict)

    def refresh_name(self, new_name, from_dir=False):
        """
        Updates name.

        Parameters
        ----------
        new_name : string
            New name or path.
        from_dir : boolean, optional
            If True, new_name is a path. The default is False.

        Returns
        -------
        None.

        """
        if from_dir:
            # convert dir to name
            path = Path(new_name)
            new_stem = path.stem
            new_name = re.sub("[0-9]{12}_","", new_stem) # remove date
        if self.name != new_name:
            self.name = new_name
            print(f"\nGenerator {self.name} is created.")
            parent_dir = Path(f"output/{self.name}")
            self.dir = Path(f"output/{self.name}/compounds")
            make_dir(parent_dir, self.dir)

    def refresh_parameters(self):
        """
        Compresses parameters into a dictionary for exporting.

        Returns
        -------
        None.

        """
        self.parameters = {
            "glycan_lmin": self.glycan_lmin,
            "glycan_lmax": self.glycan_lmax,
            "peptide_lmin": self.peptide_lmin,
            "peptide_lmax": self.peptide_lmax,
            "glycan_units": self.glycan_units,
            "peptide_residues": self.peptide_residues,
            "bridge_peptides": self.bridge_peptides,
            "num_polymerisations": self.num_polymerisations,
            "polymerisation_types": self.polymerisation_types,
            "modifications": self.modifications,
            "diffcalc_params": self.diffcalc_params,
            "diffcalc_units_min": self.diffcalc_units_min,
            "diffcalc_maxdiff_per_unit": self.diffcalc_maxdiff_per_unit
        }

    def update_from_parameters(self, loaded_parameters):
        """
        Updates parameters from a dictionary file.

        Parameters
        ----------
        loaded_parameters : dict
            Dictionary with parameters to be loaded.

        Returns
        -------
        None.

        """
        for key, value in loaded_parameters.items():
            setattr(self, key, value)

    def generate_PGN(self, start_step=1, exit_step=5):
        """
        Runs all generating functions to create a library of PGN.
        PGN is generated in 5 steps:
            1: Stem Peptide
            2: Glycan
            3: Monomer (requires 1 and 2)
            4: Modification (requires 3)
            5: Polymerisation (requires 3)

        Parameters
        ----------
        start_step : int, optional
            Indicates the integer of the first step.
            Use 1 to start from peptide generation etc. The default is 1.
        exit_step : int, optional
            Indicates the integer of the last step.
            Use 5 to exit after polymerisation. The default is 5.
        Returns
        -------
        None.

        """
        mods = [x for x in self.modifications if self.modifications[x]]
        mod_string = \
            "\n\t".join(
                "\n\t\t".join(
                    (mod, Peptidoglycan.modifications[mod])) for mod in mods)
        print(f"\nRunning generator with mods:\n\t{mod_string}")
        self.counter.disable_sleep()

        if start_step == 1 and exit_step >= start_step:
            self.PGN_dict[1] = {}  # 1 = monomers, 2 = dimers, 3 = trimers
            print("\nGenerating peptides...\t1/5")
            self.create_peptide_combs()
            start_step+=1
        if start_step == 2 and exit_step >= start_step:
            print("\nGenerating glycans...\t2/5")
            self.create_glycan_combs()
            start_step+=1
        if start_step == 3 and exit_step >= start_step:
            print("\nGenerating monomers...\t3/5")
            self.create_length_combs()
            start_step+=1
        if start_step == 4 and exit_step >= start_step:
            print("\nModifying monomers...\t4/5")
            if self.modifications["Braun LPP"]:
                self.form_LPP_dipeptide()
            if self.modifications["Alanine/Lactate Substitution"]:
                self.substitute_terminal_Ala_Lac()
            if self.modifications["Muramic Lactam"]:
                self.form_terminal_muramic_lactam()
            if self.modifications["Lactoyl Peptide"]:
                self.form_lactoyl_peptides()
            start_step+=1
        if start_step == 5 and exit_step >= start_step:
            if self.num_polymerisations > 0:
                print("\nCreating polymerisations...\t5/5")
                self.process_polymerisations()
                self.process_diffcalcs()
                for i in range(2, self.num_polymerisations+2):
                    self.create_polymerisations(i)
            print("\nCompleted...\t6/6")
            print(f"\tTotal:\t{len(self)}")
            start_step+=1

        self.counter.enable_sleep()
        msg = "\n".join(
            [f"\t{k}-mer: {len(self.PGN_dict[k])}" for k in self.PGN_dict])
        print(msg)
        self.counter.show_toast("PSN_MS2: Generator Job complete",
                                msg,
                                10)
        return msg

    def create_PGN_filter(self, PGN_components, PGN_bond, show_why=False):
        """
        Creates a function that filters PGNs by their constitutents. Used to
        filter PGN for polymerisation.

        Parameters
        ----------
        PGN_components : dict
            Dictionary describing composition restrictions imposed on a PGN unit
            in a polymerisation. Valid optional keys (except pep_len) include:
                1-5 :
                    allowed: List of amino acids allowed for this position in peptide
                    Accepts "None" and "any", "dicarboxy", "diamino".
                    bridge: Whether bridge peptides are required or forbidden.
                    bridges_allowed: List of bridge peptides allowed for this position.
                    Accepts "None" and "any".
                pep_len : List of lengths of stem peptide allowed.
                gly_len : list of lengths of glycan chain allowed.
                gly_allowed : List of glycans allowed to be in glycan chain.
                gly_rejected : List of glycans not allowed.
        PGN_bond : list
            Four part list [1,2,3,4] describing half of a bond:
                1) Integer, index of amino acid (1-5) involved. 0 indicates
                glycosidic bonds.
                2) String, describes part of amino acid or glycan involved.
                Can be "main", "bridge" or "glycan".
                3) String, describes type of bond involved. Can be "eupeptide"
                or "glycan".
                4) String, describes role in bond. "NH2"/"COOH" are used for
                peptide bonds and "Acc"/"Dnr" are used for glycosidic.
        show_why : boolean, optional
            If True. PGN_filter will print reason for rejections (for bugfixing).
            The default is False.
        Returns
        -------
        function
            PGN filtering function.

        """
        def filter_PGN(PGN):
            bond_idx = PGN_bond[0]
            # check peptide length
            if "pep_len" in PGN_components:
                if PGN.peptide_len(PGN.index) not in PGN_components["pep_len"]:
                    if show_why:
                        print("peptide length",
                              PGN_components["pep_len"])
                    return False
            # check glycan length
            if "gly_len" in PGN_components:
                if PGN.glycan_len(PGN.index) not in PGN_components["gly_len"]:
                    if show_why: print("glycan length",
                                       PGN_components["gly_len"])
                    return False
            # check valid glycan types
            if "gly_allowed" in PGN_components:
                #!!! KIV change this to PGN.check_glycan_identity
                gly_range = PGN_components["gly_allowed"]
                if gly_range != ["any"]:
                    gly_code = PGN.gly_code
                    if not all(gly in gly_range for gly in gly_code):
                        return False
            if "gly_rejected" in PGN_components:
                #!!! KIV change this to PGN.check_glycan_identity
                gly_range = PGN_components["gly_rejected"]
                gly_code = PGN.gly_code
                if any(gly in gly_range for gly in gly_code):
                    return False
            # check AAs
            for idx in range(1, self.peptide_lmax+1):
                # check if AA
                if idx == bond_idx:
                    if not PGN.check_AA_identity(idx, ["any"]):
                        return False
                if idx not in PGN_components:
                    continue
                # check bridge peptides
                if "bridge" in PGN_components[idx] and\
                        PGN_components[idx]["bridge"] != PGN.has_bridge(idx):
                    if show_why: print(idx,"bridge")
                    return False
                if "bridges_allowed" in PGN_components[idx] and\
                    not PGN.check_bridge_identity(
                        idx, PGN_components[idx]["bridges_allowed"]):
                    if show_why: print(idx,"bridges_allowed")
                    return False
                # check AA
                if "allowed" in PGN_components[idx] and\
                    not PGN.check_AA_identity(
                        idx, PGN_components[idx]["allowed"]):
                    if show_why: print(idx,"AA_allowed")
                    return False
            return True
        return filter_PGN

    # %%% Length

    def set_length(self, t, Lmin, Lmax):
        """
        Sets the minimum and maximum length of the peptide and glycan.

        Parameters
        ----------
        t : string
            t = peptide or glycan.
        Lmin : int
        Lmax : int

        Raises
        ------
        InputError
            Raises an error if lmin > lmax and if t is not recognised.

        Returns
        -------
        msg : string
            Description for GUI.

        """
        # sanitize
        if Lmin < 0:
            Lmin = 0
        if Lmax < 0:
            Lmax = 0
        if Lmin > Lmax:
            raise InputError("Length",
                             "Min. length greater than max. length.")
        if t not in ("peptide", "glycan"):
            raise InputError("Length",
                             f"Type must be 'peptide' or 'glycan' not {t}.")
        elif t == "peptide":
            self.peptide_lmin, self.peptide_lmax = Lmin, min(
                Lmax, GEN_MAX_PEPTIDE_LENGTH)
        else:
            self.glycan_lmin, self.glycan_lmax = Lmin, Lmax
        msg = f"Length range for {t} set to {Lmin}-{Lmax}"
        print(msg)
        return msg

    def create_length_combs(self):
        """
        Combines glycans and peptides to form monomers.

        Returns
        -------
        None.

        """
        length_combs = itertools.product(
            range(self.glycan_lmin, self.glycan_lmax+1),
            range(self.peptide_lmin, self.peptide_lmax+1))
        # Sort lists
        for glycan_len in range(self.glycan_lmin, self.glycan_lmax+1):
            if glycan_len > 0:
                self.glycan_combs[glycan_len].sort(key=lambda x: x.smiles)
        for peptide_len in range(self.peptide_lmin, self.peptide_lmax+1):
            if peptide_len > 0:
                self.peptide_combs[peptide_len].sort(key=lambda x: x.smiles)

        for glycan_len, peptide_len in length_combs:
            if glycan_len == 0 and peptide_len == 0:
                # NULL
                continue
            elif glycan_len > 0 and peptide_len == 0:
                # Glycans only
                to_add = [Peptidoglycan(
                    x.smiles, glycan=x,
                    modifications=["Amidase"])
                    for x in self.glycan_combs[glycan_len]]
            elif glycan_len == 0 and peptide_len > 0:
                # Peptides only
                to_add = [Peptidoglycan(
                    x.smiles, peptide=x)
                    for x in self.peptide_combs[peptide_len]]
                for PGN in to_add:
                    PGN.create_amidase_product()  # remove lactate
                    if self.modifications["EPase P1"] or \
                            self.modifications["EPase P2"]:
                        PGN.check_endopeptidase()
            else:
                if self.modifications["EPase P1"] or \
                        self.modifications["EPase P2"]:
                    valid_peptides = [peptide for peptide in
                                      self.peptide_combs[peptide_len]
                                      if peptide.len_null == 0]
                else:
                    valid_peptides = self.peptide_combs[peptide_len]
                # Typical PGNs, combine with __and__
                total = len(
                    self.glycan_combs[glycan_len]) * len(valid_peptides)
                if total == 0:
                    continue
                cache_key = create_mol_cache_key(
                    self.glycan_combs[glycan_len]+valid_peptides)
                to_add = cache_PGN_generation(self,
                                              total,
                                              cache_key,
                                              self.glycan_combs[glycan_len],
                                              valid_peptides)
            to_update = {
                PGN.InchiKey: PGN for PGN in to_add
                if PGN.mMass > GEN_MIN_PRODUCT_MASS}
            self.PGN_dict[1].update(to_update)
            print(
                f"\tG{glycan_len} & P{peptide_len}:\t{len(to_update)} monomers")
        print(f"\tTotal:\t{len(self.PGN_dict[1])} monomers")

    # %%% Peptides

    def check_if_valid_AA(self, AA_lst):
        """
        Checks if all given amino acid residues are valid.

        Parameters
        ----------
        AA_lst : list
            List of amino acids.

        Returns
        -------
        accepted : list
            List of valid amino acids.
        rejected : list
            List of invalid amino acids.

        """

        accepted = [x for x in AA_lst if x in Generator.valid_AAs]
        rejected = [x for x in AA_lst if x not in Generator.valid_AAs]
        return accepted, rejected

    def check_if_valid_peptide(self, chain):
        """
        Checks if all amino acids in peptide chain is valid.

        Parameters
        ----------
        chain : list
            List of amino acids

        Returns
        -------
        boolean
            True if all amino acids are valid.

        """
        return all(x in Generator.valid_AAs for x in chain)

    def create_peptide_filter(self, idx, AA_lst):
        """
        Creates a function that filters off peptides that do not have AAs from
        given list at given index.

        Parameters
        ----------
        idx : int
            Position to filter at.
        AA_lst : list
            List of amino acids to be retained.

        Returns
        -------
        function
            Peptide-filtering function.

        """
        def filter_peptide(peptide):
            return peptide.peptide_code[idx-1] in AA_lst
        return filter_peptide

    def get_peptide_residues(self):
        """
        Returns a dictionary which shows which amino acids can be used for at
        every position.

        Returns
        -------
        peptide_residues : dict
            Dictionary which shows list of amino acids for each position.

        """
        peptide_residues = {idx: value["AAs"]
                            for idx, value in self.peptide_residues.items()}
        return peptide_residues

    def set_peptide_residues(self, idx, AA_lst, precAA_lst=None):
        """
        Sets the amino acids that can be used at a certain position and with or
        without any conditions on the preceding amino acid.

        Parameters
        ----------
        idx : int
            Position of amino acids
        AA_lst : list
            List of amino acids
        precAA_lst : list, optional
            List of amino acids that must precede this list of amino acids.
            Only valid for positions 2-5. The default is None.

        Raises
        ------
        InputError
            Raises error if idx is greater than maximum length of peptide.

        Returns
        -------
        msg : string
            Description for GUI.

        """
        msg_parts = [] #parts of msg
        if type(idx) == int:
            if precAA_lst is None or idx == 0:
                condition = None
            else:
                accepted_pre, rejected_pre = self.check_if_valid_AA(precAA_lst)
                condition = accepted_pre
                if len(rejected_pre) > 0:
                    rejectedpre_msg = f"\t!Rejected precAAs at Pos[{idx}]:\t{rejected_pre}"
                    print(rejectedpre_msg)
                    msg_parts.append(rejectedpre_msg)
            accepted, rejected = self.check_if_valid_AA(AA_lst)
            self.peptide_residues[idx] = {"Condition": condition,
                                          "AAs": sorted(accepted)}
            num_added = len(accepted)
            added_msg = f"{num_added} AAs at Pos[{idx}]:\t{accepted}"
            if condition:
                append = f" if preceding AA = {accepted_pre}"
                added_msg += append
            print(added_msg)
            msg_parts.append(added_msg)
            if len(rejected) > 0:
                rejected_msg = f"\t!Rejected AAs at Pos[{idx}]:\t{rejected}"
                print(rejected_msg)
                msg_parts.append(rejected_msg)
            msg = "\n".join(msg_parts)
            return msg
        else:
            raise InputError("Stem Peptide",
                             f"Index {idx} not integer.")

    def create_peptide_combs(self):
        """
        Creates all possible combinations of stem peptide inclusive of the
        lactate moiety and dummy AAs for bridge peptides.

        Returns
        -------
        None.

        """
        self.peptide_combs = {}
        for idx in range(1, self.peptide_lmax+1):
            if idx == 1:
                self.peptide_combs[1] = [Lac+AMINO_ACID_DB[x]
                                         for x in self.peptide_residues[1]["AAs"]]
            else:
                condition = self.peptide_residues[idx]["Condition"]
                AAs = self.peptide_residues[idx]["AAs"]
                if condition:
                    peptide_filter = self.create_peptide_filter(
                        idx-1, condition)
                    valid_peptides = [peptide for peptide
                                      in self.peptide_combs[idx-1]
                                      if peptide_filter(peptide)]
                else:
                    valid_peptides = self.peptide_combs[idx-1]
                combs = itertools.product(valid_peptides, AAs)
                self.peptide_combs[idx] = [seq+AMINO_ACID_DB[ext]
                                           for seq, ext in combs]
            # Add dummy bridge peptides
            if idx in self.bridge_peptides:
                for grp in self.bridge_peptides[idx]:
                    # only combine peptides with valid side chain
                    valid_AAs = self.bridge_peptides[idx][grp]["valid_AAs"]
                    valid = self.filter_peptides_by_side_chain_smarts(
                        self.peptide_combs[idx], idx, grp, valid_AAs)
                    ext = self.dummy_AAs[idx]  # use dummy AA temporarily
                    if grp == "COOH":
                        self.peptide_combs[idx].extend(
                            [seq & ext for seq in valid])
                    else:
                        self.peptide_combs[idx].extend(
                            [ext-seq for seq in valid])
            # Endopeptidase modifications
            if self.modifications["EPase P1"] and idx == 1:
                print("\tAdding null monopeptide...")
                self.peptide_combs[idx].append(NULL_AA)
            if self.modifications["EPase P2"] and idx == 2:
                print("\tAdding null dipeptide...")
                self.peptide_combs[idx].append(NULL_AA+NULL_AA)
            num_peptides = len(self.peptide_combs[idx])
            print(f"\tP{idx}:\t{num_peptides} ")
        if self.bridge_peptides != {}:
            self.swap_dummyAA_for_bridges()

    # %%%% Bridge Peptides

    def set_bridge_peptides(self, idx, grp,
                           bridge_lst,  # each bridge peptide will be a list
                           valid_AAs="all"):
        """
        Sets the possible bridge peptides at a certain position attaching at
        either COOH or NH2.

        Parameters
        ----------
        idx : int
            Position of bridge peptide.
        grp : string
            Indicates attachment at either COOH or NH2.
        bridge_lst : list
            List of bridge peptides. Each peptide is a list, starting from the
            innermost amino acid to the outermost.
            i.e. ["Ala","Ala"], ["Ser", "Ala", "Thr", "Ala"]
        valid_AAs : list, optional
            List of amino acids that bridge peptide can connect to.
            The default is "all".

        Raises
        ------
        InputError
            Raises error if idx is not an integer or if group is neither NH2
            or COOH.

        Returns
        -------
        msg : string
            Description for GUI.

        """
        # Assemble bridge peptides, generate connecting rxn in advance
        valid = []
        rejected = []
        bridge_lst.sort()
        msg_parts = [] #parts of msg
        if type(idx) == int and grp in ("NH2", "COOH"):
            for pep in bridge_lst:
                if self.check_if_valid_peptide(pep):
                    valid.append(pep)
                else:
                    rejected.append(pep)
            if idx not in self.bridge_peptides:
                self.bridge_peptides[idx] = {}
            self.bridge_peptides[idx][grp] = {"valid_AAs": valid_AAs,
                                             "bridges": valid}
            num_added = len(valid)
            added_msg = f"{num_added} bridge peptides at AA[{idx}]-{grp} for {valid_AAs}:\
                \n\t{bridge_lst}"
            print(added_msg)
            msg_parts.append(added_msg)
            if len(rejected) > 0:
                rejected_msg = f"\t!Rejected bridge peptides at AA[{idx}]-{grp}:\t{rejected}"
                print(rejected_msg)
                msg_parts.append(rejected_msg)
            msg = "\n".join(msg_parts)
            return msg
        elif type(idx) != int:
            raise InputError("Bridge",
                             f"Index {idx} not integer.")
        else:
            raise InputError("Bridge",
                             f"Group must be 'NH2' or 'COOH'; not {grp}.")

    def process_bridge_peptides(self):
        """
        Assemble bridge peptides and create swapping reactions.
        A dummy AA is added during initial PGN creation. Dummy AAs are swapped
        for the bridge peptide later on.

        Returns
        -------
        None.

        """
        self.bridge_peptides_processed = {}
        for idx in self.bridge_peptides:
            dAA = self.dummy_AAs[idx]
            self.bridge_peptides_processed[idx] = {}
            for grp in self.bridge_peptides[idx]:
                self.bridge_peptides_processed[idx][grp] = []
                bridges = self.bridge_peptides[idx][grp]["bridges"]
                for bridge in bridges:
                    if type(bridge) != list:
                        bridge = list(bridge)  # force list datatype
                    # get connecting AA for this chain
                    cAA_code = bridge[0]
                    cAA = AMINO_ACID_DB[cAA_code]
                    # create swapping reaction
                    rxn = self.create_swapping_rxn(dAA, cAA, grp)
                    # assemble chain
                    if grp == "NH2":
                        bridge_mol = functools.reduce(
                            operator.add,
                            [AMINO_ACID_DB[x] for x in reversed(bridge)])
                    else:
                        bridge_mol = functools.reduce(
                            operator.add,
                            [AMINO_ACID_DB[x] for x in bridge])
                    self.bridge_peptides_processed[idx][grp].append(
                        (bridge_mol, rxn))

    def create_swapping_rxn(self, dAA, cAA, grp):
        """
        Creates a swapping reaction to replace the dummy amino acid with the
        bridge peptide.

        Parameters
        ----------
        dAA : AminoAcid
            Dummy amino acid.
        cAA : AminoAcid
            Connecting amino acid of bridge peptide - this is the innermost acid.
        grp : string
            Indicates attachment at either COOH or NH2.

        Returns
        -------
        rxn : ChemicalReaction
            Chemical reaction that will swap the dummy amino acid with the
            bridge peptide.

        """
        if grp == "COOH":
            # Use dAA C-terminus
            rxt1 = f"[NX3:51]-{dAA.eupeptide_smarts['COOH']}-[OH:52]"
            # Use cAA N-terminus
            rxt2 = f"[NH2:53]-{cAA.eupeptide_smarts['NH2']}-[OH,NX3:54]"
            pdt1 = f"[NX3:51]-{cAA.eupeptide_smarts['NH2']}-[OH,NX3:54]"
            pdt2 = f"[NH2:53]-{dAA.eupeptide_smarts['COOH']}-[OH:52]"
        elif grp == "NH2":
            # Use dAA N-terminus
            rxt1 = f"[NH2:51]-{dAA.eupeptide_smarts['NH2']}-[NX3:52]"
            # Use cAA C-terminus
            rxt2 = f"[NX3:53]-{cAA.eupeptide_smarts['COOH']}-[OH:54]"
            pdt1 = f"[NX3:53]-{cAA.eupeptide_smarts['COOH']}-[NX3:52]"
            pdt2 = f"[NH2:51]-{dAA.eupeptide_smarts['NH2']}-[OH:54]"
        rxn_smarts = f"{rxt1}.{rxt2}>>{pdt1}.{pdt2}"
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        return rxn

    def swap_dummy_AAs_with_bridge(self, rxn, peptide, bridge, lat_idx, grp):
        """
        Swaps the dummy amino acid with the bridge peptide and creates the
        corresponding Peptide.

        Parameters
        ----------
        rxn : ChemicalReaction
            Chemical reaction that will swap the dummy amino acid with the
            bridge peptide.
        peptide : Peptide
            Peptide molecule containing dummy amino acid.
        bridge : Peptide
            Bridge peptide.
        lat_idx : int
            Position of AA in which bridge peptide is attached to.
        grp : string
            Indicates attachment at either COOH or NH2.

        Raises
        ------
        NoReactionError
            No reaction occurred.
        TooManyPdtsError
            Too many possible products occurred.

        Returns
        -------
        pdt : Peptide
            Resulting peptide.

        """
        pdts1_smiles, _ = run_two_pdt_rxn(
            rxn, peptide.mol, bridge.mol)
        if len(pdts1_smiles) == 0:
            raise NoReactionError(peptide, bridge, f"swap {lat_idx}", rxn)
        elif len(pdts1_smiles) > 1:
            raise TooManyPdtsError(
                peptide, bridge, f"swap {lat_idx}", rxn, pdts1_smiles)
        pdt_smile = pdts1_smiles[0]
        # update side_chain at pos of side chain
        side_chain_dict = peptide.generate_updated_side_chain(
            bridge, grp, idx=lat_idx)
        # Create New Peptide
        pdt = Peptide(smiles=pdt_smile,
                      peptide_code=copy.deepcopy(
                          peptide.peptide_code),
                      eupeptide_smarts=copy.deepcopy(
                          peptide.eupeptide_smarts),
                      side_chain_smarts=copy.deepcopy(
                          peptide.side_chain_smarts),
                      side_chain_dict=side_chain_dict)
        return pdt

    def swap_dummyAA_for_bridges(self):
        """
        Swaps out all dummy acids for all bridges.

        Returns
        -------
        None.

        """
        print("Swapping dummy AAs...")
        self.process_bridge_peptides()
        # iterate through all peptides
        for idx in range(1, self.peptide_lmax+1):
            to_add = []  # swapped peptides
            to_remove = []  # remove original with dummy AAs
            for ori_peptide in self.peptide_combs[idx]:
                children = [[ori_peptide]]  # swapped peptides for this peptide
                # iterate through side chains
                for lat_idx in ori_peptide.side_chain_dict:
                    grp, _ = ori_peptide.side_chain_dict[lat_idx]
                    if grp == None:
                        continue
                    else:
                        children.append([])
                        for peptide in children[-2]:
                            children[-1].extend(
                                [self.swap_dummy_AAs_with_bridge(
                                    rxn, peptide, bridge, lat_idx, grp)
                                 for bridge, rxn in
                                 self.bridge_peptides_processed[lat_idx][grp]])
                if len(children) == 1:
                    continue  # no side chains substituted
                else:
                    to_remove.append(ori_peptide)
                    to_add.extend(children[-1])
            # after each length, clean up
            for r_peptide in to_remove:
                self.peptide_combs[idx].remove(r_peptide)
            self.peptide_combs[idx].extend(to_add)
            num_peptides = len(self.peptide_combs[idx])
            print(f"\tP{idx}:\t{num_peptides}")

    def filter_peptides_by_side_chain_smarts(self, peptides_lst,
                                             idx, grp, valid_AAs):
        """
        Filters a list of stem peptides to obtain stem peptides that contain
        specific amino acids (with correct connectivity COOH/NH2) at specific
        index.

        Parameters
        ----------
        peptides_lst : list
            List of stem peptides
        idx : int
            Position of amino acid in which bridge peptide is attached to.
        grp : string
            Indicates attachment at either COOH or NH2.
        valid_AAs : list
            List of valid amino acids.

        Returns
        -------
        list
            List of valid stem peptides

        """
        return [x for x in peptides_lst if
                type(x.side_chain_smarts[grp]) == str
                and x.check_AA_identity(idx, valid_AAs)]

    # %%% Glycans

    def check_if_valid_glycan(self, t, gly_lst):
        """
        Checks a list to see if glycans provided are valid.

        Parameters
        ----------
        t : string
            Indicates type of glycan; can be both -"All", glucosamine - "Glc"
            or muramic acid - "Mur".
        gly_lst : List
            List of glycans.

        Returns
        -------
        accepted : list
            List of valid glycans.
        rejected : list
            List of invalid glycans.

        """
        if t == "All":
            valid_glycans_lst = set(flatten(
                Generator.valid_glycans.values()))
        else:
            valid_glycans_lst = Generator.valid_glycans[t]
        accepted = [x for x in gly_lst if x in valid_glycans_lst]
        rejected = [x for x in gly_lst if x not in valid_glycans_lst]
        return accepted, rejected

    def set_glycan_units(self, t, gly_lst):
        """
        Set the glycans that can be used for Glc-type or Mur- type positions.

        Parameters
        ----------
        t : string
            Indicates type of glycan; glucosamine - "Glc" or muramic acid - "Mur".
        gly_lst : List
            List of glycans.

        Raises
        ------
        InputError
            Raises error if incorrect type given.

        Returns
        -------
        msg : string
            Description for GUI.

        """
        msg_parts = [] #parts of msg
        if t in ("Glc", "Mur"):
            accepted, rejected = self.check_if_valid_glycan(t, gly_lst)
            self.glycan_units[t] = sorted(accepted)
            num_added = len(self.glycan_units[t])
            added_msg = f"{num_added} {t}-type glycans:\t{self.glycan_units[t]}"
            print(added_msg)
            msg_parts.append(added_msg)
            if len(rejected) > 0:
                rejected_msg = f"\t Rejected {t}-type glycans:\t{rejected}"
                print(rejected_msg)
                msg_parts.append(rejected_msg)
            # assemble msg
            msg = "\n".join(msg_parts)
            return msg
        else:
            raise InputError("Glycan",
                             f"Type must be 'Glc' or 'Mur'; not {t}.")

    def create_glycan_combs(self):
        """
        Generates all possible glycan combinations. Glycans are built from the
        reducing terminus with "Mur" type glycans and alternate between "Glc"
        type and "Mur" type. Glycans are linked with 1,4 glycosidic bonds.

        Returns
        -------
        None.

        """
        for t in ["Glc", "Mur"]:
            assert t in self.glycan_units, f"No {t}-type glycans provided."
            assert len(self.glycan_units[t]
                       ) > 0, f"No {t}-type glycans provided."

        self.glycan_combs = {1: [GLYCAN_DB[x]
                                 for x in self.glycan_units["Mur"]]}
        for idx in range(2, self.glycan_lmax+1):
            t = "Glc" if idx % 2 == 0 else "Mur"
            combs = itertools.product(
                self.glycan_combs[idx-1], self.glycan_units[t])
            to_add = []
            for seq, ext in combs:
                try:
                    # fails with anhydro
                    to_add.append(GLYCAN_DB[ext]+seq)
                except NoReactionError:
                    continue  # skip
            self.glycan_combs[idx] = to_add
            print(f"\tG{idx}:\t{len(to_add)}")

    # %%% Difference Calc

    def get_diff_threshold(self, k):
        """
        Get the difference threshold.

        Parameters
        ----------
        k : integer
            The number of muropeptide units (aka PGN units) in a multimer.

        Returns
        -------
        diff : integer
            The difference threshold

        """
        # n = number of PGN units
        if (k-1) > len(self.diffcalc_maxdiff_per_unit):
            # exceeds list, return last value
            diff = self.diffcalc_maxdiff_per_unit[-1]
        else:
            diff = self.diffcalc_maxdiff_per_unit[k-2]
        return diff

    def set_diffcalc_units(self, min_units, maxdiffs):
        """
        Sets the values for minimum number.

        Parameters
        ----------
        min_units : integer
            Minimum k-mer for diff. calculation to be used, monomers excluded.
        maxdiffs : float
            Maximum number of diff. tolerated per PGN unit.

        Returns
        -------
        None.

        """
        self.diffcalc_units_min = min_units
        self.diffcalc_maxdiff_per_unit = maxdiffs
        msg = f'''Applying diff. calc. to >={self.diffcalc_units_min}-mers
        Max diff. per unit = {self.diffcalc_maxdiff_per_unit}'''
        print(msg)
        return msg

    def set_diffcalc_param(self, gly_len, gly_range,
                           pep_len, pep_range,
                           bridge_range=None,
                           polymerisation_range=None):
        """
        Sets the values used to describe a typical PGN unit.

        Parameters
        ----------
        gly_len : list
            List of glycan lengths.
        gly_range : list
            List of glycans (both types).
        pep_len : list
            List of peptide lengths.
        pep_range : list
            List of amino acids (stem peptide only).
        bridge_range : list, optional
            List of bridge peptides.
            If set to None, bridge peptides will not be checked.
        polymerisation_range : list, optional
            List of polymerisation names.
            If set to None, polymerisation types will not be checked.

        Returns
        -------
        None.

        """
        self.diffcalc_params = {
            "gly_len": gly_len,
            "gly_range": gly_range,
            "pep_len": pep_len,
            "pep_range": pep_range,
            "bridge_range": bridge_range,
            "polymerisation_range": polymerisation_range
        }
        msg = f'''Canonical muropeptide has...
        Glycan lengths:\t{gly_len}
        Glycans:\t{gly_range}
        Peptide lengths:\t{pep_len}
        Amino acids:\t{pep_range}
        Bridge peptides:\t{bridge_range}
        Polymerisation types:\t{polymerisation_range}
        '''
        print(msg)
        return msg

    def process_diffcalcs(self, show_breakdown=False):
        """
        From the values provided, generates a function that calculates the no.
        of differences between a Peptidoglycan Molecule and a typical PGN unit.

        Parameters
        ----------
        show_breakdown : boolean, optional
            If True, prints a breakdown of differences. The default is False.

        Returns
        -------
        calc_difference: function
            Calculates no. of differences.

        """
        if len(self.diffcalc_params) == 0 or\
                self.diffcalc_units_min > self.num_polymerisations+1:
            return
        gly_len = self.diffcalc_params["gly_len"]
        gly_range = self.diffcalc_params["gly_range"]
        pep_len = self.diffcalc_params["pep_len"]
        pep_range = self.diffcalc_params["pep_range"]
        bridge_range = self.diffcalc_params["bridge_range"]
        polymerisation_range = self.diffcalc_params["polymerisation_range"]

        def calc_difference(PGN):
            # difference in glycans
            gly_len_diff = min(
                [abs(PGN.glycan_len(PGN.index) - num) for num in gly_len])
            gly_range_diff = PGN.glycan_comp_diffcalc(gly_range)
            # difference in peptides
            pep_len_diff = min(
                [abs(PGN.peptide_len(PGN.index) - num)for num in pep_len])
            pep_range_diff = PGN.peptide_comp_diffcalc(pep_range)
            # difference in bridge peptides
            bridge_range_diff = PGN.bridge_comp_diffcalc(bridge_range)
            # difference in polymerisation types
            if polymerisation_range is None or PGN.polymerisation_type in polymerisation_range:
                polymerisation_type_diff = 0
            else:
                polymerisation_type_diff = 1
            # total
            diffs = (gly_len_diff, gly_range_diff,
                     pep_len_diff, pep_range_diff,
                     bridge_range_diff, polymerisation_type_diff
                     )
            total_diff = sum(diffs)
            if show_breakdown:
                diffs_name = ("GL", "GR", "PL", "PR", "SC", "CL")
                print(f"Diffcalc for {PGN.name} index-{PGN.index}:")
                print('\t'.join(diffs_name))
                print('\t'.join(map(str, diffs)))
            if PGN.index > 1:
                linked_diff = calc_difference(PGN.linked_PGN)
                total_diff += linked_diff
            return total_diff
        self.diffcalc = calc_difference
        return calc_difference

    # %%% Polymerisations

    def set_num_polymerisations(self, num):
        """
        Sets the maximum number of polymerisations used.

        Parameters
        ----------
        num : integer
            Maximum number of polymerisations

        Returns
        -------
        None.

        """
        self.num_polymerisations = num
        msg = f"Max polymerisations: {self.num_polymerisations}"
        print(msg)
        return msg

    def generate_polymerisation_name(self, P1_bond, P2_bond):
        """
        Generates name of polymerisation based on bonding.

        Parameters
        ----------
        P1_bond : list
            Four part list [1,2,3,4] describing half of a bond:
                1) Integer, index of amino acid (1-5) involved. 0 indicates
                glycosidic bonds.
                2) String, describes part of amino acid or glycan involved.
                Can be "main", "bridge" or "glycan".
                3) String, describes type of bond involved. Can be "eupeptide"
                or "glycan".
                4) String, describes role in bond. "NH2"/"COOH" are used for
                peptide bonds and "Acc"/"Dnr" are used for glycosidic.
        P2_bond : list
            Four part list [1,2,3,4] describing half of a bond. See above.

        Returns
        -------
        polymerisation_name : string
            Name of polymerisation.

        """

        def generate_polymerisation_name_part(Pn_bond):
            num, chain_type, bond, grp = Pn_bond
            chain_str = "" if chain_type in ("main", "glycan") else "br"
            type_str = "" if bond in ("eupeptide", "glycan") else "s"
            if num == 0:
                bond_str = f"G{chain_str}{type_str}"
            else:
                bond_str = f"{num}{chain_str}{type_str}"
            return grp, bond_str

        grp1, P1_bond_name = generate_polymerisation_name_part(P1_bond)
        grp2, P2_bond_name = generate_polymerisation_name_part(P2_bond)
        polymerisation_name = "-".join((P1_bond_name, P2_bond_name))
        return polymerisation_name

    def check_polymerisation_composition_restrictions(self, Pn):
        """
        Checks the polymerisation restrictions validity.

        Parameters
        ----------
        Pn : dict
            Dictionary describing composition restrictions imposed on a PGN unit
            in a polymerisation. Valid optional keys (except pep_len) include:
                1-5 :
                    allowed: List of amino acids allowed for this position in peptide
                    Accepts "None" and "any", "dicarboxy", "diamino".
                    bridge: Whether bridge peptides are required or forbidden.
                    bridges_allowed: List of bridge peptides allowed for this position.
                    Accepts "None" and "any".
                pep_len : List of lengths of stem peptide allowed.
                gly_len : list of lengths of glycan chain allowed.
                gly_allowed : List of glycans allowed to be in glycan chain.
                gly_rejected : List of glycans not allowed.

        Raises
        ------
        InputError
            Raises an error if peptide length is not provided.

        Returns
        -------
        no_errors : boolean
            If errors present: True; else False.

        """
        no_errors = True
        # compulsory condition
        if "pep_len" not in Pn:
            print(Pn)
            raise InputError("Polymerisation",
                             "Peptide length condition required.")
        # check glycan allowed list is valid
        gly_lst = Pn.get("gly_allowed", [])
        valid, invalid = self.check_if_valid_glycan("All", gly_lst)
        if len(invalid) > 0:
            invalid_msg = f"\t!Invalid accepted glycan:\t{invalid}"
            print(invalid_msg)
            no_errors = False
        # check glycan rejected list is valid
        gly_lst = Pn.get("gly_rejected", [])
        valid, invalid = self.check_if_valid_glycan("All", gly_lst)
        if len(invalid) > 0:
            invalid_msg = f"\t!Invalid rejected glycan:\t{invalid}"
            print(invalid_msg)
            no_errors = False
        # check peptide restrictions are valid
        max_Pn = max(Pn["pep_len"])
        aa_lst_Pn = flatten(
            Pn.get(i, {}).get("allowed", [])
            for i in range(1, max_Pn+1))
        valid, invalid = self.check_if_valid_AA(aa_lst_Pn)
        if len(invalid) > 0:
            invalid_msg = f"\t!Invalid restricted AA:\t{invalid}"
            print(invalid_msg)
            no_errors = False
        return no_errors

    def check_polymerisation_bond_restrictions(self, i, Pn_bond):
        """
        Checks whether the bond restrictions are valid.

        Parameters
        ----------
        i : integer
            Integer indicating Peptide 1 or Peptide 2.
        Pn_bond : list
            Four part list [1,2,3,4] describing half of a bond:
                1) Integer, index of amino acid (1-5) involved. 0 indicates
                glycosidic bonds.
                2) String, describes part of amino acid or glycan involved.
                Can be "main", "bridge" or "glycan".
                3) String, describes type of bond involved. Can be "eupeptide"
                or "glycan".
                4) String, describes role in bond. "NH2"/"COOH" are used for
                peptide bonds and "Acc"/"Dnr" are used for glycosidic.

        Raises
        ------
        InputError
            Raises an error if any of the list is invalid.

        Returns
        -------
        bool
            Returns True if valid.

        """
        num, chain_type, bond, grp = Pn_bond
        if num > self.peptide_lmax:
            raise InputError("Polymerisation",
                             "Index {num} greater than max length.")
        elif chain_type not in ("main", "bridge", "glycan"):
            raise InputError("Polymerisation",
                             f'''Type must be 'main' or 'bridge' or 'glycan'';
                             not {chain_type}.''')
        elif i == 2 and chain_type == "bridge":
            raise InputError("Polymerisation",
                             "Type must be 'main' on 2nd PGN.")
        elif bond not in ("eupeptide", "side", "glycan"):
            raise InputError("Polymerisation",
                             f'''Bond must be 'eupeptide' or 'side' or 'glycan';
                             not {bond}.''')
        elif grp not in ("NH2", "COOH", "Acc", "Dnr"):
            raise InputError("Polymerisation",
                             f'''Group must be 'NH2' or 'COOH' or "Acc" or "Dnr";
                             not {grp}.''')
        return True

    def check_polymerisation_bond_match(self, P1_bond, P2_bond):
        """
        Checks if the two bond descriptions match.

        Parameters
        ----------
        P1_bond : list
            Four part list [1,2,3,4] describing half of a bond:
                1) Integer, index of amino acid (1-5) involved. 0 indicates
                glycosidic bonds.
                2) String, describes part of amino acid or glycan involved.
                Can be "main", "bridge" or "glycan".
                3) String, describes type of bond involved. Can be "eupeptide"
                or "glycan".
                4) String, describes role in bond. "NH2"/"COOH" are used for
                peptide bonds and "Acc"/"Dnr" are used for glycosidic.
        P2_bond : list
            Four part list [1,2,3,4] describing half of a bond. See above.

        Raises
        ------
        InputError
            Raises an error if the two bonds do not form a valid pair.

        Returns
        -------
        bool
            Returns True if valid.

        """
        num1, chain_type1, bond1, grp1 = P1_bond
        num2, chain_type2, bond2, grp2 = P2_bond
        valid_pair1 = ("NH2", "COOH")
        valid_pair2 = ("Acc", "Dnr")
        for valid_pair in (valid_pair1, valid_pair2):
            if grp1 in valid_pair and grp2 in valid_pair:
                return True
        raise InputError("Polymerisation",
                         f'''{grp1} and {grp2} do not match.''')

    def clear_polymerisations(self):
        """
        Clears existing polymerisations.

        Returns
        -------
        None.

        """
        num_cleared = len(self.polymerisation_types)
        self.polymerisation_types = []  # criteria for polymerisation
        self.polymerisation_types_processed = None
        msg  = f"Cleared {num_cleared} polymerisation types"
        return msg

    def set_polymerisation_types(self, P1, P2, P1_bond, P2_bond, kmer_range=None):
        """
        Checks, names and sets the polymerisation described. P1 and P2 refer to
        polymerised peptide chains, and the bridge peptide (if any) is assumed to be on P1.
        P1/2 describes the requirements of the PGN unit.
        P1/2_bond describes the requirements of the polymerisation bond.

        Parameters
        ----------
        P1 : dict
            Dictionary describing composition restrictions imposed on a PGN unit
            in a polymerisation. Valid optional keys (except pep_len) include:
                1-5 :
                    allowed: List of amino acids allowed for this position in peptide
                    Accepts "None" and "any", "dicarboxy", "diamino".
                    bridge: Whether bridge peptides are required or forbidden.
                    bridges_allowed: List of bridge peptides allowed for this position.
                    Accepts "None" and "any".
                pep_len : List of lengths of stem peptide allowed.
                gly_len : list of lengths of glycan chain allowed.
                gly_allowed : List of glycans allowed to be in glycan chain.
                gly_rejected : List of glycans not allowed.
        P2 : dict
            Dictionary describing composition restrictions imposed on a PGN unit
            in a polymerisation. See above.
        P1_bond : list
            Four part list [1,2,3,4] describing half of a bond:
                1) Integer, index of amino acid (1-5) involved. 0 indicates
                glycosidic bonds.
                2) String, describes part of amino acid or glycan involved.
                Can be "main", "bridge" or "glycan".
                3) String, describes type of bond involved. Can be "eupeptide"
                or "glycan".
                4) String, describes role in bond. "NH2"/"COOH" are used for
                peptide bonds and "Acc"/"Dnr" are used for glycosidic.
        P2_bond : list
            Four part list [1,2,3,4] describing half of a bond. See above.
        kmer_range : list, optional.
            List of integers describing the valid kmer range of polymerisation.
            Defaults to [1,99].

        Raises
        ------
        InputError
            Raises an error if either PGN composition or bond is incorrect.

        Returns
        -------
        name : string
            Name of polymerisation type.

        """
        if kmer_range is None:
            kmer_range = [1, 99]
        # check composition; raise InputError
        b1 = self.check_polymerisation_composition_restrictions(P1)
        b2 = self.check_polymerisation_composition_restrictions(P2)
        # check bond; raise InputError
        self.check_polymerisation_bond_restrictions(1, P1_bond)
        self.check_polymerisation_bond_restrictions(2, P2_bond)
        name = self.generate_polymerisation_name(P1_bond, P2_bond)
        # check if groups match; raise InputError
        self.check_polymerisation_bond_match(P1_bond, P2_bond)
        self.polymerisation_types.append(
            [name, P1, P2, P1_bond, P2_bond, kmer_range])
        if b1 and b2:
            msg = f"Added {name} polymerisation with range:{kmer_range}"
        else:
            msg = f"\t!Added {name} polymerisation with errors."
        print(msg)
        return msg

    def process_polymerisations(self, show_why=False):
        """
        Creates and saves a function that filters PGN that can participate as
        either P1 or P2 in each polymerisation.

        Parameters
        ----------
        show_why : boolean, optional
            If True. PGN_filter will print reason for rejections (for bugfixing).
            The default is False.

        Returns
        -------
        None.

        """
        self.polymerisation_types_processed = []
        for name, P1, P2, P1_bond, P2_bond, kmer_range in self.polymerisation_types:
            P1_filter = self.create_PGN_filter(P1, P1_bond, show_why)
            P2_filter = self.create_PGN_filter(P2, P2_bond, show_why)
            self.polymerisation_types_processed.append(
                [name, P1_filter, P2_filter, P1_bond, P2_bond, kmer_range])
            print(f"\tProcessed {name} polymerisation")

    def create_polymerisations(self, k):
        """
        Creates and saves polymerised PGN units for k-mers;
        where k = no. of total PGN units.

        Parameters
        ----------
        k : integer
            No. of total PGN units in new multimer.

        Returns
        -------
        None.

        """
        self.PGN_dict[k] = {}
        to_add = []
        print(f"\nPolymerising to create {k}-mers")
        for line in self.polymerisation_types_processed:
            name, P1_filter, P2_filter, P1_bond, P2_bond, kmer_range = line
            if min(kmer_range) > k or k > max(kmer_range):
                print(f"\t{name} polymerisation out of range for {k}-mers")
                continue  # skip out of range
            # P1 always monomers
            valid_P1 = [PGN for PGN in self.PGN_dict[1].values()
                        if P1_filter(PGN)]
            valid_P2 = [PGN for PGN in self.PGN_dict[k-1].values()
                        if P2_filter(PGN)]
            len_P1 = len(valid_P1)
            len_P2 = len(valid_P2)
            total = len_P1*len_P2
            if total == 0:
                if len_P1 == 0:
                    print(f"\t!!No valid P1 PGN with 1 units for {name}!!")
                else:
                    print(f"\t!!No valid P2 PGN with {k-1} units for {name}!!")
                continue
            # reduce variety with diffcalc
            if k >= self.diffcalc_units_min:
                diff_threshold = self.get_diff_threshold(k)
                P1_threshold = diff_threshold
                valid_P1 = [PGN for PGN in valid_P1
                            if self.diffcalc(PGN) <= P1_threshold]
                P2_threshold = diff_threshold*(k-1)
                valid_P2 = [PGN for PGN in valid_P2
                            if self.diffcalc(PGN) <= P2_threshold]
                mod_len_P1 = len(valid_P1)
                mod_len_P2 = len(valid_P2)
                mod_total = mod_len_P1*mod_len_P2
                reduction = int(100*(total-mod_total)/total)
                print(f"\tRestricted {name} {k}-mers:")
                print(f"\t\tP1: {P1_threshold} diff(s); P2: {P2_threshold} diff(s)")
                print(f"\t\tReduced generated {k}-mers by {reduction}%")
                total = mod_total
                len_P1 = mod_len_P1
                len_P2 = mod_len_P2
            if len_P1 == 0 or len_P2 == 0:
                continue
            # check amount
            print(f"\t{name} polymerised {k}-mers: {len_P1}x{len_P2}")
            if total > self.warning_threshold:
                prompt = f'''Continue with {total} {name}? Y/N/X'''
                print(prompt)
                answer = "N"
                if answer == "N":
                    continue
                elif answer == "X":
                    return
            # sort
            valid_P1.sort(key=lambda x: x.smiles)
            valid_P2.sort(key=lambda x: x.smiles)

            cache_key = create_mol_cache_key(valid_P1+valid_P2)
            to_add.extend(
                cache_polymer_generation(
                    self, total, name, P1_bond, P2_bond,
                    cache_key, valid_P1, valid_P2))
        to_update = {PGN.InchiKey: PGN for PGN in to_add}
        self.PGN_dict[k].update(to_update)
        print(f"\tPolymerised {k}-mers:\t{len(to_update)}")

    # %%% Modifications

    def add_modifications(self, listmods):
        """
        Adds modifications to Generator. Modifications not listed are set to False.

        Parameters
        ----------
        listmods : list
            List of applicable modifications.

        Returns
        -------
        None.

        """
        self.modifications = {
            mod: mod in listmods for mod in Peptidoglycan.modifications}
        total_mods = len(self.modifications)
        mods_added = [
            mod for mod in self.modifications if self.modifications[mod]]
        msg = f"{len(mods_added)}/{total_mods} modifications applied:\n{mods_added}"
        print(msg)
        return msg

    def update_PGN_dict(self, k):
        """
        Updates all the keys for all k-mers for in-place modifications.

        Parameters
        ----------
        k : integer
           The number of muropeptide units (aka PGN units) in a multimer.

        Returns
        -------
        None.

        """
        self.PGN_dict[k] = {
            PGN.InchiKey: PGN for PGN in self.PGN_dict[k].values()}

    def reduce_terminal_muramic_acid(self):
        """
        Reduces terminal "Mur"-type glycans in monomers.

        Returns
        -------
        None.

        """
        # pick all monomers
        num_cpds = sum(PGN.reduce_terminal_muramic_acid()
                       for PGN in self.PGN_dict[1].values())
        self.update_PGN_dict(k=1)
        print(f"\tReduced terminal Mur:\t{num_cpds} monomers")

    def substitute_terminal_Ala_Lac(self):
        """
        Creates a copy of all pentapeptide PGN with terminal Ala and substitutes
        the amide bond with an ester bond.

        Returns
        -------
        None.

        """
        # pick peptides with terminal Ala (5th residue)
        valid = [copy.deepcopy(PGN) for PGN in self.PGN_dict[1].values()
                 if PGN.peptide_len(PGN.index) == 5 and PGN.peptide_code[-1] == "Ala"]
        valid_rxt = [PGN for PGN in valid if PGN.substitute_terminal_Ala_Lac()]
        num_cpds = len(valid_rxt)
        self.PGN_dict[1].update({PGN.InchiKey: PGN for PGN in valid_rxt})
        print(f"\tSubstituted terminal Ala with Lac:\t{num_cpds} monomers")

    def form_terminal_muramic_lactam(self):
        """
        Forms lactams in terminal "Mur"-type glycans in monomers.

        Returns
        -------
        None.

        """
        # pick all monomers absent peptides
        valid = [copy.deepcopy(PGN) for PGN in self.PGN_dict[1].values()
                 if PGN.peptide is None]
        valid_rxt = [
            PGN for PGN in valid if PGN.form_terminal_muramic_lactam()]
        num_cpds = len(valid_rxt)
        self.PGN_dict[1].update({PGN.InchiKey: PGN for PGN in valid_rxt})
        print(f"\tFormed muramic lactams:\t{num_cpds} monomers")

    def form_lactoyl_peptides(self):
        """
        Forms new lactoyl peptides for each peptide.

        Returns
        -------
        None.

        """
        valid = itertools.chain.from_iterable(
            self.peptide_combs[pep_len] for pep_len in range(
                max(self.peptide_lmin, 1),
                self.peptide_lmax+1))
        lactoyl_peptides = [Peptidoglycan(x.smiles, glycan=gly_Lac, peptide=x)
                            for x in valid if x.len_null == 0]
        for x in lactoyl_peptides:
            x.modifications["Lactoyl Peptide"] = True
        num_cpds = len(lactoyl_peptides)
        self.PGN_dict[1].update(
            {PGN.InchiKey: PGN for PGN in lactoyl_peptides})
        print(f"\tFormed lactoyl peptides:\t{num_cpds} monomers")

    def form_LPP_dipeptide(self):
        """
        Forms Brauns lipoprotein attachment points for monomers.

        Returns
        -------
        None.

        """
        # check if Lys, Arg are in residues
        r4_isolysine = "ε-Lys" in self.peptide_residues[4]["AAs"]
        r5_arginine = "Arg" in self.peptide_residues[5]["AAs"]
        if r4_isolysine and r5_arginine:  # search and label
            def has_LPP(PGN):
                boolean = (PGN.peptide_len(PGN.index) == 5
                           and PGN.peptide_code[3] == "ε-Lys"
                           and PGN.peptide_code[4] == "Arg")
                return boolean
            valid = [PGN for PGN in self.PGN_dict[1].values() if has_LPP(PGN)]
            for PGN in valid:
                PGN.modifications["Braun LPP"] = True
            print(f"Labeled {len(valid)} Braun LPP monomers")
        else:  # add dipeptide to valid length 3
            total = 0
            isolys = AMINO_ACID_DB["ε-Lys"]
            arg = AMINO_ACID_DB["Arg"]
            LPP = isolys + arg
            for glycan_len in range(self.glycan_lmin, self.glycan_lmax+1):
                pep_LPP = [x + LPP for x in self.peptide_combs[3]
                           if x.len_null == 0]
                if glycan_len == 0:
                    # Peptides only
                    to_add = [Peptidoglycan(x.smiles, peptide=x)
                              for x in pep_LPP]
                    for PGN in to_add:
                        PGN.create_amidase_product()  # remove lactate
                        PGN.modifications["Braun LPP"] = True
                else:
                    # Typical PGNs, combine with __and__
                    total = len(self.glycan_combs[glycan_len]) * len(pep_LPP)
                    if total == 0:
                        continue
                    cache_key = create_mol_cache_key(
                        self.glycan_combs[glycan_len]+pep_LPP)
                    to_add = cache_PGN_generation(self,
                                                  total,
                                                  cache_key,
                                                  self.glycan_combs[glycan_len],
                                                  pep_LPP)
                    for PGN in to_add:
                        PGN.modifications["Braun LPP"] = True
                total += len(to_add)
                to_update = {
                    PGN.InchiKey: PGN for PGN in to_add
                    if PGN.mMass > GEN_MIN_PRODUCT_MASS}
                self.PGN_dict[1].update(to_update)
            print(f"\tBraun LPP: \t {total} monomers")

    # %%% Importing/Exporting

    def get_kmers_range(self, kmers):
        """
        Processes input to a valid k-mers range.

        Parameters
        ----------
        kmers : None, list or integer
            If None, gives the maximum range. If integer, converts it to a list.
            If list, no changes.

        Returns
        -------
        list
            Valid k-mer ranges for the exporting functions.

        """
        if kmers is None:
            return range(1, self.num_polymerisations+2)
        elif type(kmers) == list:
            return kmers
        else:
            return [kmers]

    def export_dataframes(self, formats, overwrite=False, kmers=None):
        """
        Creates and exports dataframes in various formats for k=mers.

        Parameters
        ----------
        formats : list
            List of formats, currently supports xlsx, sdf, csv, txt.
        overwrite : boolean, optional
            If True, create and overwrite existing dataframe. If False,
            use existing dataframe (unless empty). The default is False.
        kmers : None, list or integer
            If None, gives the maximum range. If integer, converts it to a list.
            If list, no changes.

        Returns
        -------
        None.

        """
        kmers = self.get_kmers_range(kmers)
        if overwrite or self.master_dataframe is None:
            self.create_dataframe()
        for ext in formats:
            if ext == "xlsx":
                self.export_dataframe_as_xlsx(kmers)
            elif ext == "sdf":
                self.export_dataframe_as_sdf(kmers)
            elif ext == "csv":
                self.export_dataframe_as_csv(kmers)
            elif ext == "txt":
                self.export_MS1_lib(kmers=kmers)
            else:
                print(f"Unknown ext. {ext}")

    def create_dataframe(self):
        """
        Creates dataframe from PGN_dict.

        Returns
        -------
        None.

        """
        print("\nCreating dataframes...")

        def create_df_part(k):
            selected = self.PGN_dict[k].values()
            k_df = pd.DataFrame([PGN.convert_to_dict() for PGN in selected])
            return k_df
        self.master_dataframe = pd.concat(
            create_df_part(k) for k in self.PGN_dict)

    def export_dataframe_as_xlsx(self, kmers):
        """
        Exports the dataframe in spreadsheet format.

        Parameters
        ----------
        kmers : list
            Valid k-mer range to export.

        Returns
        -------
        None.

        """
        print("\nExporting dataframes as xlsx...")
        # Export xlsx
        xlsx_filename = self.dir/f"{TIME_STRING}_{self.name}.xlsx"
        with pd.ExcelWriter(xlsx_filename) as writer:
            for k in range(1, self.num_polymerisations+2):
                if k not in kmers:
                    continue
                sheet_name = f"{k}-mers"
                df = self.master_dataframe[
                    self.master_dataframe["PGN Units"] == k]
                df = df.drop("Mol", axis="columns")
                df.to_excel(writer, sheet_name=sheet_name, index=False)
        print(f"\tDataframe saved as:\t{xlsx_filename}")

    def export_dataframe_as_csv(self, kmers, export_structure=True):
        """
        Exports the dataframe in csv format. To be depreciated.

        Parameters
        ----------
        kmers : list
            Valid k-mer range to export.
        export_structure : boolean, optional
            If True, exports structures as MolBlocks. The default is True.

        Returns
        -------
        None.

        """
        print("\nExporting dataframes as csv...")
        csv_filename = self.dir/f"{TIME_STRING}_{self.name}.csv"
        df = self.master_dataframe.drop("Mol", axis="columns")
        if export_structure:
            df["Structure"] = self.master_dataframe["Mol"].apply(
                Chem.MolToMolBlock)
            csv_filename = self.dir/"{TIME_STRING}_{self.name}_withstruct.csv"
        else:
            csv_filename = self.dir/f"{TIME_STRING}_{self.name}_withstruct.csv"
        df.to_csv(csv_filename, index=False)
        print(f"\tDataframe saved as:\t{csv_filename}")

    def export_dataframe_as_sdf(self, kmers):
        """
        Exports the dataframe in sdf format. To be depreciated.

        Parameters
        ----------
        kmers : list
            Valid k-mer range to export.

        Returns
        -------
        None.

        """
        print("\nExporting dataframes as sdf...")
        sdf_filename = self.dir/f"{TIME_STRING}_{self.name}.sdf"
        PandasTools.WriteSDF(self.master_dataframe, sdf_filename,
                             molColName="Mol", idName="Name",
                             properties=["SMILES", "Formula", "Monoisotopic Mass"])
        print(f"\tDataframe saved as:\t{sdf_filename}")

    def export_MS1_lib(self, kmers):
        """
        Exports the dataframe in txt format. This file can be used with MS-DIAL.

        Parameters
        ----------
        kmers : list
            Valid k-mer range to export.

        Returns
        -------
        None.

        """
        headers_lst = ["NAME", "MZ", "RT", "Adduct", "InChIKey",
                       "Formula", "SMILES", "Ontology",  # end of defaults
                       "PGN units", "Polymerised PGN", "Modifications",
                       "Glycan units", "Glycan",
                       "Stem Peptide units", "Stem Peptide",
                       "Monoisotopic Mass", "clogP"]
        txt_filename = self.dir/f"{TIME_STRING}_{self.name}.txt"
        with open(txt_filename, "w", encoding="ASCII") as f:
            headers = "\t".join(headers_lst)
            f.write(headers)
            f.write("\n")
            for k in self.PGN_dict:
                if k not in kmers:
                    continue
                for InchiKey in self.PGN_dict[k]:
                    PGN = self.PGN_dict[k][InchiKey]
                    info_lst = [PGN.formula, PGN.smiles,
                                PGN.ontology,
                                PGN.index, PGN.linked_PGN, PGN.mods,
                                PGN.glycan_len(), PGN.glycan_name(),
                                PGN.peptide_len(), PGN.peptide_name(),
                                PGN.mMass, PGN.clogP]
                    info = "\t".join(map(str, info_lst))
                    for adduct, mz in PGN.mz.items():
                        if adduct in self.output_adducts:
                            mz_info_lst = [PGN.name, mz, "0", adduct, InchiKey]
                            mz_info = "\t".join(map(str, mz_info_lst))
                            PGN_info = "\t".join([mz_info, info])
                            f.write(PGN_info)
                            f.write("\n")
            f.write("\n")
        print(f"\nMS1 library saved as:\t{txt_filename}")
        print(f"\tAdducts:{self.output_adducts}")

    def export_reference_dict(self):
        """
        Exports the keys for PGN_dict (InChiKeys)

        Returns
        -------
        None.

        """
        # is this useful?
        df = pd.DataFrame.from_dict(self.reference_dict, orient="index")
        df.index.rename("Batch No.")
        xlsx_filename = self.dir/f"{TIME_STRING}_{self.name}_InchiKeys.xlsx"
        df.to_excel(xlsx_filename)
        print(f"\nReference dict saved as:\t{xlsx_filename}")

    def export_PGN_as_pickle(self, kmers=None):
        """
        Exports PGN_dict in pickle format as batches.

        Parameters
        ----------
        kmers : None, list or integer
            If None, gives the maximum range. If integer, converts it to a list.
            If list, no changes.

        Returns
        -------
        pickle_filename : Path
            Path of pickle file.

        """
        kmers = self.get_kmers_range(kmers)
        pickle_filename = self.dir/f"{TIME_STRING}_{self.name}.pickle"
        self.reference_dict = {}
        print(f"\nSaving as {pickle_filename}")
        with open(pickle_filename, "wb") as f:
            idx = 0
            for k in self.PGN_dict:
                if k not in kmers:
                    continue
                InchiKeys = list(self.PGN_dict[k].keys())
                num_cpds = len(InchiKeys)
                batches = int(np.ceil(num_cpds/BATCH_SIZE))
                for i in range(0, batches):
                    nmin = i*BATCH_SIZE
                    nmax = (i+1)*BATCH_SIZE
                    selected = {InchiKey: self.PGN_dict[k][InchiKey]
                                for InchiKey in InchiKeys[nmin:nmax]}
                    pickle.dump(selected, f)
                    print(f"\tSaved batch {idx+1}:\t{k}-mers")
                    self.reference_dict[f"Batch{idx+1}"] = InchiKeys[nmin:nmax]
                    idx += 1
        self.export_reference_dict()
        print(f"Compound data saved as:\t{pickle_filename}")
        return pickle_filename

    def import_PGN_as_pickle(self, pickle_filename,
                             inherit_name=True, overwrite=True,
                             kmers=None):
        """
        Imports pickle file to use as PGN_dict.

        Parameters
        ----------
        pickle_filename : Path
            Path of pickle file.
        inherit_name : bool, optional
            Determines if name is inherited from pickle file.
            The default is True.
        overwrite : boolean, optional
            If True, clears existing data. The default is True.
        kmers : None, list or integer
            If None, gives the maximum range. If integer, converts it to a list.
            If list, no changes.

        Returns
        -------
        None.

        """
        print(f"\nLoading {pickle_filename}")
        if overwrite:
            print("\tOverwriting current data.")
            self.PGN_dict = {}
        total_cpds = 0
        subtotal_cpds = 0
        with open(pickle_filename, "rb") as f:
            k = 1
            idx = 0
            while True:
                try:
                    d = pickle.load(f)
                    num_cpds = len(d)
                    print(f"\tLoaded batch {idx+1}:\t{num_cpds} {k}-mers")
                    if k not in self.PGN_dict:
                        self.PGN_dict[k] = {}
                    self.PGN_dict[k].update(d)
                    idx += 1
                    total_cpds += num_cpds
                    subtotal_cpds += num_cpds
                    if num_cpds < BATCH_SIZE:
                        print(f"\t\tSubtotal:\t{subtotal_cpds} {k}-mers")
                        k += 1
                        subtotal_cpds = 0
                except EOFError:  # end of file
                    break
        print(f"Loaded {total_cpds} cpds.")
        if inherit_name:
            self.refresh_name(pickle_filename, from_dir=True)

    def export_settings_as_pickle(self):
        """
        Export Generator settings in a pickle format.

        Returns
        -------
        settings_filename : Path
            Path of settings file.

        """
        settings_filename = self.dir/f"{TIME_STRING}_{self.name}.settings"
        self.refresh_parameters()
        with open(settings_filename, "wb") as f:
            pickle.dump(self.parameters, f)
        print(f"\nGenerator settings saved as:\t{settings_filename}")
        return settings_filename

    def import_settings_as_pickle(self, settings_filename, inherit_name=True):
        """
        Import Generator settings in a pickle format.

        Parameters
        ----------
        settings_filename : Path
            Path of settings file.
        inherit_name : bool, optional
            Determines if name is inherited from pickle file.
            The default is True.

        Returns
        -------
        None.


        """
        with open(settings_filename, "rb") as f:
            loaded_parameters = pickle.load(f)
        self.update_from_parameters(loaded_parameters)
        print(f"\nGenerator settings loaded from:\t{settings_filename}")
        if inherit_name:
            self.refresh_name(settings_filename, from_dir=True)

    def export_settings_as_yaml(self):
        """
        Export Generator settings in a yaml format.

        Returns
        -------
        yaml_filename :  Path
            Path of settings file.

        """
        yaml_filename = self.dir/f"{TIME_STRING}_{self.name}.yaml"
        self.refresh_parameters()
        with open(yaml_filename, "w") as file:
            yaml.dump(self.parameters, file)
        print(f"\nGenerator settings saved as:\t{yaml_filename}")
        return yaml_filename

    def import_settings_as_yaml(self, yaml_filename, inherit_name=True):
        """
        Import Generator settings in a yaml format.

        Parameters
        ----------
        yaml_filename :  Path
            Path of settings file.
        inherit_name : bool, optional
            Determines if name is inherited from yaml file.
            The default is True.

        Returns
        -------
        None.

        """
        with open(yaml_filename, "r") as f:
            loaded_parameters = yaml.safe_load(f)
        self.update_from_parameters(loaded_parameters)
        print(f"\nGenerator settings loaded from:\t{yaml_filename}")
        if inherit_name:
            self.refresh_name(yaml_filename, from_dir=True)

    def export_settings_as_image(self):
        """
        Export Generator settings as an graphical summary. Not all details
        are captured, however.

        Returns
        -------
        None.

        """
        print("\nCreating graphical summary...")
        illustrator = Illustrator(ILL_SIZE, ILL_FONTSIZE, BOND_CHAR)
        illustrator.plot_all(self)
        svg_filename = self.dir/f"{TIME_STRING}_{self.name}_graphical_summary.svg"
        illustrator.save_figure(svg_filename)
