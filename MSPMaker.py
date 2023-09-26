'''
Fragmenter and MSPMaker Classes
env = chem
'''

import os
import pandas as pd
import numpy as np
import pickle
import yaml
import re

from pathlib import Path

from Reactions import PEPTIDE_FRAG_TYPES, GLYCAN_FRAG_TYPES
from Molecules import Molecule
from Molecules import run_fragment_rxn

from Common import BATCH_SIZE, PRECISION_MASS, PRECISION_INTENSITY,\
    TIME_STRING, MEMORY, OUTPUT_ADDUCTS
from Common import Counter
from Common import make_dir, sigmoid

# %% Globals

# Default values are False
# Terminal = Not valid precursor; Sum = unique structures not req.
# Type = valid for precursors of adduct type

FRAG_TYPES_DICT = {
    # "M+[Na]": {"int_factor": 1, "Terminal": True, "Type": "Na"},
    # "Lac+[Na]": {"int_factor": 0.75, "Type": "Na"},
    # "M+[K]": {"int_factor": 1, "Terminal": True, "Type": "K"},
    # "Lac+[K]": {"int_factor": 0.75, "Type": "K"},
    "Lac": {"int_factor": 0.9},
    "Gly-Lac": {"int_factor": 0.1},
    "Lac[r]": {"int_factor": 0.6},
    "Gly-Lac[r]": {"int_factor": 0.1},
    "Gly-B": {"int_factor": 0.6},
    "Gly-Y": {"int_factor": 0.25},
    "Gly-C": {"int_factor": 0.05},
    "Gly-Z": {"int_factor": 0.3},
    "Gly-Y[r]": {"int_factor": 0.4},
    "Gly-Z[r]": {"int_factor": 0.5},
    "Pep-b": {"int_factor": 0.9},
    "Pep-y": {"int_factor": 0.9},
    "Pep-q1": {"int_factor": 0.4},
    "Pep-q2": {"int_factor": 0.6},
    "Pep-e1": {"int_factor": 0.15},
    "Pep-e2": {"int_factor": 0.2},
    "Pep-K1": {"int_factor": 0.8, "Terminal": True},
    "Mur-1": {"int_factor": 0.4, "Terminal": True},  # 138
    "Glc-1": {"int_factor": 0.3, "Terminal": True},  # 186
    "Glc-2": {"int_factor": 0.3},  # 168
    "Glc-3": {"int_factor": 0.9, "Terminal": True},  # 126
    "Acyl-H2O": {"int_factor": 0.12, "Terminal": True, "Sum": True},
    "Acyl-NH3": {"int_factor": 0.18, "Terminal": True, "Sum": True}
}

FRAG_RXN_MAPPING = {
    # "Lac [Na]": {"moiety": "glycan", "frag_types": ["Lac+[Na]", None],
    #              "precursor_frag_type": ["Parent", "Lac+[Na]"]},
    # "Lac [K]": {"moiety": "glycan", "frag_types": ["Lac+[K]", None],
    #             "precursor_frag_type": ["Parent", "Lac+[K]"]},
    # "Adduct [Na]": {"moiety": "glycan", "frag_types": ["M+[Na]", None],
    #                 "precursor_frag_type": ["Parent"]},
    # "Adduct [K]": {"moiety": "glycan", "frag_types": ["M+[K]", None],
    #                "precursor_frag_type": ["Parent"]},
    "Gly. B/Y": {"moiety": "glycan", "frag_types": ["Gly-B", "Gly-Y"]},
    "Gly. C/Z": {"moiety": "glycan", "frag_types": ["Gly-C", "Gly-Z"]},
    "Lac": {"moiety": "glycan", "frag_types": ["Lac", "Gly-Lac"]},
    "Gly. B/Y[r]": {"moiety": "glycan", "frag_types": ["Gly-B", "Gly-Y[r]"]},
    "Gly. C/Z[r]": {"moiety": "glycan", "frag_types": ["Gly-C", "Gly-Z[r]"]},
    "Lac[r]": {"moiety": "glycan", "frag_types": ["Lac[r]", "Gly-Lac[r]"]},
    "Glc 1": {"moiety": "glycan", "frag_types": ["Glc-1", "Glc-2"]},
    "Glc 2": {"moiety": "glycan", "frag_types": ["Glc-3", None]},
    "Mur 1": {"moiety": "glycan", "frag_types": ["Mur-1", None]},
    "Peptide y": {"moiety": "peptide", "frag_types": ["Pep-y", None]},
    "Peptide b": {"moiety": "peptide", "frag_types": ["Pep-b", None]},
    "Peptide q": {"moiety": "peptide", "frag_types": ["Pep-q1", "Pep-q2"]},
    "Peptide e": {"moiety": "peptide", "frag_types": ["Pep-e1", "Pep-e2"]},
    "Peptide K": {"moiety": "peptide", "frag_types": ["Pep-K1", None]},
    # "-H2O(r5)": {"moiety": "peptide", "frag_types": ["Acyl-H2O", None]},
    # "-NH3(r5)": {"moiety": "peptide", "frag_types": ["Acyl-NH3", None]},
    "-H2O": {"moiety": "peptide", "frag_types": ["Acyl-H2O", None]},
    "-NH3": {"moiety": "peptide", "frag_types": ["Acyl-NH3", None]}
}

FRAG_RXN_ORDER = [
    #"Lac [Na]", "Lac [K]", "Adduct [Na]", "Adduct [K]",
    "Gly. B/Y", "Gly. C/Z", "Lac",
    "Gly. B/Y[r]", "Gly. C/Z[r]", "Lac[r]",
    "Glc 1", "Glc 2", "Mur 1",
    "Peptide y", "Peptide b", "Peptide q", "Peptide e", "Peptide K",
    #"-H2O(r5)", "-NH3(r5)",
    "-H2O", "-NH3"]

FRAG_MODEL = {
    "max_times": 3,
    "min_precursor_percent": 0.03,
    "min_product_mz": 100,
    "min_product_int": 0.01,
    "max_num_ions_export": 50,  # kept in simple msp
    "max_num_ions_save": 200,  # kept in full peaklist
    "iontypes_dict": FRAG_TYPES_DICT,
    "mapping": FRAG_RXN_MAPPING}

MSP_OMIT_IONS_GREATER_THAN_PRECURSOR = True
MSP_OMIT_MZ_THRESHOLD = 5
MSP_PRECURSOR_MZ_MAX = 2500
MSP_PRECURSOR_MZ_MIN = 250


def mztxt(num):
    '''Formats m/z (num) to specified precision.'''
    return f"{num:.{PRECISION_MASS}f}"


def inttxt(num):
    '''Formats intensity (num) to specified precision.'''
    return f"{num:.{PRECISION_INTENSITY}f}"


def get_adduct_type(adduct):
    '''Returns type of adduct.'''
    if adduct == "[M+Na]+":
        return "Na"
    elif adduct == "[M+K]+":
        return "K"
    else:
        return "H"

# %% Cached


@MEMORY.cache(ignore=["frag", "PGN"])
def cache_spectral_data(frag, PGN, InchiKey, model=FRAG_MODEL):
    '''Caches simulated PGN spectra for current fragmentation model
    using InchiKey as hash.'''
    frag.load_PGN(PGN, InchiKey)
    frag.fragment_all()
    return frag.create_full_peaklist_dataframe()

# %% Fragmenter


class Fragmenter:
    '''Handles and logs fragmentation of PGN'''
    PARAM_FILEPATH = Path("settings/frag/settings.yaml")
    ADDUCT_FORMS = ["[M+H]+", "[M+2H]2+", "[M+3H]3+", "[M+Na]+"]

    def __init__(self, name="test",
                 frag_mode="formula",
                 network=None,
                 quiet=False,
                 aggregate=True,
                 adjust_int_each_rd=True,
                 min_precursor_int=None,
                 min_precursor_percent=None,
                 min_product_mz=None,
                 min_product_int=None):
        """
        Handles and logs fragmentation of PGN.

        Parameters
        ----------
        name : string, optional
            Name of Fragmenter (used in filenames). The default is "test".
        frag_mode : string, optional
            Type of calculations employed to calculate fragment intensity.
            The default is "formula".
        network : Network object, optional
            Network used to calculate fragment intensity (if required).
            The default is None.
        quiet : bool, optional
            If True, most messages will be muted. The default is False.
        aggregate : bool, optional
            If True, aggregation of fragments will be done during df creation.
        adjust_int_each_rd : bool, optional
            If True, maximum precursor intensity will be adjusted up each round.
        min_precursor_int: float, optional
            Min. intensity for precursor to be fragmented.
            If None, defaults to parameter defined in FRAG_MODEL.
        min_precursor_percent: float, optional
            If adjust_int_each_rd = True, min_precursor_int will be increased
            after each rd to min_precursor_percent * maximum fragment intensity
            If None, defaults to parameter defined in FRAG_MODEL.
        min_product_mz: float, optional
            Min. mz for precursor ion to be saved.
            If None, defaults to parameter defined in FRAG_MODEL.
        min_product_int: float, optional
            Min. intensity for product ion to be saved.
            If None, defaults to parameter defined in FRAG_MODEL.

        Returns
        -------
        None.

        """

        self.quiet = quiet
        self.min_precursor_percent = min_precursor_percent
        self.min_precursor_int = min_precursor_int
        self.min_product_int = min_product_int
        self.min_product_mz = min_product_mz
        self.update_null_variables()  # update variables from FRAG_MODEL
        self.name = None
        self.dir = None
        self.refresh_name(name)
        self.adduct_type = None
        self.frag_mode = frag_mode
        self.aggregate = aggregate
        self.adjust_int_each_rd = adjust_int_each_rd
        self.add_debugging()  # initialise debug parameters

        if self.frag_mode == "formula":
            # trial and error formula
            self.int_adj = FRAG_TYPES_DICT
            self.frag_network = None
        elif self.frag_mode == "show_all":
            # show all possible peaks
            self.int_adj = None
            self.frag_network = network
            self.adjust_int_each_rd = False
        else:
            # neural network with intermediates
            self.frag_network = network
            self.int_adj = None
        self.counter = Counter()
        self.peaklist_df_full = None
        self.peaklist_df_simple = None
        self.PGN_fragments = None
        self.log_settings()
        if not self.quiet:
            print(
                f'''\nFragmenting up to {FRAG_MODEL["max_times"]} times with
                min. precursor % int. {self.min_precursor_int} and
                min. product m/z {self.min_product_mz} and
                int. {self.min_product_int}''')

    def add_debugging(self,
                      debug_fragments=None,
                      debug_frag_types=None,
                      debug_reactions=None,
                      debug_mz=None):
        """
        Updates debugging parameters. Called during __init__.

        Parameters
        ----------
        debug_fragments : List, optional
            A list of SMILES. If provided, a debugging message will be printed
            whenever the fragment is formed. The default is None.
        debug_frag_types : List, optional
            A list of fragment_types. If provided, a debugging message will be
            printed whenever the fragment_type is generated. The default is None.
        debug_reactions: list, optional
            A list of reactions. If provided, a debugging message will be
            printed whenever the reaction produces fragments. The default is None.
        debug_mz : List, optional
            A list of mz. If provided, a debugging message will be
            printed whenever a fragment with said mz is formed. The default is None.

        Returns
        -------
        None.

        """

        if debug_fragments is not None:
            print("Accepted debugging fragments as SMILES")
            self.debug_fragments = debug_fragments
        else:
            self.debug_fragments = []
        if debug_frag_types is not None:
            print("Accepted debugging frag_types as strings")
            self.debug_frag_types = debug_frag_types
        else:
            self.debug_frag_types = []
        if debug_reactions is not None:
            print("Accepted debugging reaction as strings")
            self.debug_reactions = debug_reactions
        else:
            self.debug_reactions = []
        if debug_mz is not None:
            print("Accepted debugging mz as floats")
            self.debug_mz = debug_mz
        else:
            self.debug_mz = []

    def update_null_variables(self):
        """
        Update variables from FRAG_MODEL if variable == None.

        Returns
        -------
        None.

        """
        if self.min_precursor_percent is None:
            self.min_precursor_percent = FRAG_MODEL["min_precursor_percent"]
        if self.min_precursor_int is None:
            self.min_precursor_int = FRAG_MODEL["min_precursor_percent"]
        if self.min_product_int is None:
            self.min_product_int = FRAG_MODEL["min_product_int"]
        if self.min_product_mz is None:
            self.min_product_mz = FRAG_MODEL["min_product_mz"]

    def log_settings(self):
        """
        Compares current fragmentation model with previous.
        Saves a copy of current model if modified.

        Returns
        -------
        None.

        """
        with open(Fragmenter.PARAM_FILEPATH, "r") as f:
            prev_settings = yaml.safe_load(f)
        if prev_settings != FRAG_MODEL:
            yaml_filename = Path(f"settings/frag/{TIME_STRING}_settings.yaml")
            with open(yaml_filename, "w") as f:
                yaml.dump(FRAG_MODEL, f)
            print(f"Saved new fragmentation parameters as:\t{yaml_filename}")
            with open(Fragmenter.PARAM_FILEPATH, "w") as f:
                yaml.dump(FRAG_MODEL, f)

    def refresh_name(self, new_name, from_dir=False):
        """
        Refreshes the name of Fragmenter.

        Parameters
        ----------
        new_name : String or Path
            New name/path.
        from_dir : bool, optional
            Indicates if name is from path. The default is False.

        Returns
        -------
        None.

        """
        if from_dir:
            # convert dir to name
            path = Path(new_name)
            new_stem = path.stem
            new_name = re.sub("[0-9]{12}_", "", new_stem)  # remove date
        if self.name != new_name:
            self.name = new_name
            if not self.quiet:
                print(f"\nFragmenter {self.name} is created.")
            self.dir = Path(f"output/{self.name}/peaklists")
            make_dir(f"output/{self.name}", self.dir)

    def load_PGN(self, PGN, PGN_name="test", adduct_type="[M+H]+",
                 auto_process=False, **kwargs):
        """
        Loads a PGN molecule for fragmentation.

        Parameters
        ----------
        PGN : Molecule object, optional
            Molecule object that represents PGN. Can be loaded later.
            The default is None.
        PGN_name : String, optional
            Name of PGN (used in filenames). The default is "test".
        adduct_type: string, optional
            Type of adduct. Only used when frag_mode = "network"
            The default is "[M+H]+".
        auto_process : bool, optional
            If True, loaded PGN will be fragmented on initiation. Simple and
            full dataframes will be saved. The default is False.
        **kwargs :
            Passed to export_peaklist.

        Returns
        -------
        None.

        """
        self.PGN_name = PGN_name
        self.PGN = PGN
        if self.PGN.mz == {}:
            self.PGN.update_mz()
        self.count = 1
        self.adduct_type = adduct_type
        self.PGN_fragments = {self.PGN.smiles:
                              {"idx": 0,
                               "gen": 0,
                               "mz": 0,
                               "frag_type": "Parent",
                               "p_idx": None,
                               "adduct_type": "H",
                               "terminal": False,
                               "molecule": self.PGN,
                               "formula": self.PGN.formula,
                               "neutral": self.PGN.neutral_form,
                               "neutral_key": self.PGN.neutral_InchiKey,
                               "intensity": 1}}
        self.PGN_fragment_history = {}  # record precursors for each frag
        self.PGN_fragment_int = {}  # record intensities for each parent
        self.peaklist_df_simple = None
        self.peaklist_df_full = None
        self.create_full_peaklist_dataframe(
            False, False, False, False)  # initialise df
        self.PGN_params = self.PGN.property_array
        if auto_process:
            self.auto_process(adduct_type=adduct_type, **kwargs)

    def auto_process(self, printout=True, **kwargs):
        """
        Automatically processes a molecule, fragmenting it and optionally, saves
        detailed and simple peaklists.

        Parameters
        ----------
        printout : Boolean, optional
            Determines if peaklists will be saved. The default is True.
        **kwargs :
            Passed to export_peaklist.

        Returns
        -------
        None.

        """
        self.fragment_all()
        if printout:
            self.export_peaklist(overwrite=True, simple=False, **kwargs)
            self.export_peaklist(simple=True, **kwargs)


# %%% Fragmentation

    def derive_fragment_intensity(self, frag_type, f_smiles,
                                  num_f, attributes):
        """
        Returns relative fragment intensity.

        Parameters
        ----------
        frag_type : String
            Name of fragment type.
        f_smiles : String
            Fragment SMILES
        num_f : Integer
            Total number of fragments generated with this fragment.
        attributes : Dict
            A dictionary storing attributes of fragment

        Returns
        -------
        Float
            Relative fragment intensity

        """

        if self.frag_mode == "formula":
            int_factor = self.int_adj[frag_type].get("int_factor", 1)
            mass_ratio = attributes["f_molecule"].mMass / \
                attributes["p_molecule"].mMass
            peptide_ratio = attributes["f_molecule"].num_peptide_bond/mass_ratio
            ratio = sigmoid(max(mass_ratio, peptide_ratio))
            f_int = int_factor*attributes["p_intensity"]*ratio/num_f
        elif self.frag_mode == "show_all":
            return 1
        else:
            # Prepare parameters
            f_params = attributes["f_molecule"].morgan_fingerprint(
                3, 128)  # 128
            # p_params = attributes["p_molecule"].morgan_fingerprint(3, 128) # 128
            # nonf_params = p_params - f_params  # 9
            # extra_params = np.array([1/num_f, attributes["p_intensity"]])  # 2
            params = f_params  # 256
            cat = {"f": frag_type, "p": attributes["p_frag_type"],
                   "adduct": self.adduct_type}
            int_factor = self.frag_network.calculate_value(
                params, cat=cat)
            f_int = int_factor
        return np.round(f_int, PRECISION_INTENSITY)

    def eligible_precursors(self, frag_rxn, precursor_frag_type):
        """
        Filters all stored fragments to find valid precursors.

        Parameters
        ----------
        frag_rxn : String
            Name of fragmentation reaction.
        precursor_frag_type : List or None
            List of valid frag_types which precursor must be from. If None,
            restricts to proton-only frag_types.

        Returns
        -------
        Set
            A set containing all eligible precursors.

        """
        df = self.peaklist_df_full
        precursor_done = ~df["SMILES"].isin(
            self.PGN_fragment_history[frag_rxn])
        gen_valid = df["gen"] < FRAG_MODEL["max_times"]
        int_valid = df["intensity"] > self.min_precursor_int
        not_terminal_fragment = ~df["terminal"]

        if precursor_frag_type is None:  # filter for H+ adducts
            eligible_fragment = df["adduct_type"] == "H"
        else:
            eligible_fragment = df["frag_type"].isin(precursor_frag_type)

        precursor_df = df[precursor_done &
                          gen_valid &
                          int_valid &
                          not_terminal_fragment &
                          eligible_fragment]
        return precursor_df

    def fragment_once(self, frag_rxn, moiety, frag_types,
                      precursor_frag_type=None):
        """
        Carries out fragmentation reaction on all eligible precursors.

        Parameters
        ----------
        moiety : string
            'Peptide' or 'Glycan'.
        frag_rxn : string
            Name of fragmentation reaction.
        type1 : string
            Name of first fragment.
        type2 : string
            Name of second fragment.

        Raises
        ------
        ValueError
            Moiety is not valid.

        Returns
        -------
        bool
            Returns True if fragmentation proceeded.
            Returns False if no fragmentation occurred.

        """
        if frag_rxn not in self.PGN_fragment_history:
            self.PGN_fragment_history[frag_rxn] = set()

        precursor_df = self.eligible_precursors(frag_rxn, precursor_frag_type)
        if len(precursor_df) == 0:
            return False
        else:
            self.PGN_fragment_history[frag_rxn] |= set(precursor_df["SMILES"])

        if moiety == "peptide":
            rxn = PEPTIDE_FRAG_TYPES[frag_rxn]
        elif moiety == "glycan":
            rxn = GLYCAN_FRAG_TYPES[frag_rxn]
        else:
            raise ValueError

        def fragment_all_precursors(row):
            smiles = row["SMILES"]
            precursor = row["molecule"]
            cpd = row["neutral"]
            idx = row["idx"]

            if not cpd:
                # print(self.PGN_fragments[smiles]["idx"],
                #       self.PGN_fragments[smiles]["frag_type"])
                return

            mMass = precursor.mMass
            pdts1_smiles, pdts2_smiles = run_fragment_rxn(
                frag_rxn, mMass, rxn, cpd)
            num_f = len(pdts1_smiles)
            if frag_rxn in self.debug_reactions and num_f > 0:
                print(
                    f'''\n{num_f} fragments formed from {frag_rxn} rxn from ion {idx}''')
            for frag_type, pdt_lst in zip(frag_types,
                                          [pdts1_smiles, pdts2_smiles]):
                if frag_type is None or len(pdt_lst) == 0:
                    if frag_type in self.debug_frag_types:
                        print(
                            f'''\nNo {frag_type} pdts from ion {idx}''')
                    continue  # ignore this fragment
                elif FRAG_TYPES_DICT[frag_type].get("Sum", False):
                    # Sum fragments, no need to record unique structures
                    num_f = 1
                    f_smiles = pdt_lst[0]
                    self.record_fragment(
                        smiles, frag_type, f_smiles, num_f)
                else:
                    if frag_type in self.debug_frag_types:
                        print(
                            f'''\n{num_f} {frag_type} pdts from ion {idx}''')
                    for f_smiles in pdt_lst:
                        self.record_fragment(
                            smiles, frag_type, f_smiles, num_f)

        precursor_df.apply(fragment_all_precursors, axis="columns")

        return True

    def record_fragment(self, p_smiles, frag_type, f_smiles, num_f):
        """
        Saves fragment to self.PGN_fragments (if valid).

        Parameters
        ----------
        p_smiles : String
            Precursor SMILES
        frag_type : String
            Name of fragment_Type
        f_smiles : String
            Fragment SMILES
        num_f : Integer
            Total number of fragments generated with this fragment.

        Raises
        ------
        Exception
            Precursor not in self.PGN_fragments.
        ValueError
            Reaction / fragment is not valid.

        Returns
        -------
        bool
            Returns True if fragment is saved.
            Returns False if not.

        """
        # check if fragment already exists
        existing_fragment = f_smiles in self.PGN_fragments
        if existing_fragment:
            molecule = self.PGN_fragments[f_smiles]["molecule"]
        else:
            molecule = Molecule(f_smiles)

        # calculate attributes of fragment
        attributes = {}
        attributes["f_molecule"] = molecule

        # abort if charge incorrect or mMass too low
        if molecule.charge != 1:
            print(p_smiles, frag_type, f_smiles, num_f)
            raise ValueError(f"Created molecule with charge {molecule.charge}")
            # ignore for now
        if molecule.mMass < self.min_product_mz:
            return False

        # grab parent fragment attributes
        if p_smiles in self.PGN_fragments:
            attributes.update(
                ("p_"+k, v) for k, v in
                self.PGN_fragments[p_smiles].items()
                if k[:2] != "p_")
            # duplicate some keys
            attributes["f_gen"] = attributes["p_gen"]+1
            attributes["f_p_idx"] = attributes["p_idx"]
            attributes["f_p_molecule"] = attributes["p_molecule"]
        else:
            raise Exception("No record of precursor")

        # calculate intensity
        f_intensity = self.derive_fragment_intensity(
            frag_type, f_smiles, num_f, attributes)

        # debugging
        if frag_type in self.debug_frag_types:
            print(
                f'''{frag_type} intensity: {f_intensity} for {molecule.formula},
                {molecule.mMass}''')

        # To depreciate
        if p_smiles in self.PGN_fragment_int:
            self.PGN_fragment_int[p_smiles] += f_intensity
        else:
            self.PGN_fragment_int[p_smiles] = f_intensity

        # debugging
        if f_smiles in self.debug_fragments:
            print("\n"+f_smiles)
            parent_frag_type = attributes['p_frag_type']
            parent_idx = attributes['p_idx']
            parent_int = attributes['p_intensity']
            parent_gen = attributes["p_gen"]
            print(
                f"{frag_type} fragment, int= {f_intensity}"
                f"\nfrom {parent_frag_type} gen-{parent_gen} [{parent_idx}], int= {parent_int}")

        # discard if intensity too low
        if f_intensity < self.min_product_int:
            return False

        # accumulate if precursor already present
        if existing_fragment:
            self.PGN_fragments[f_smiles]["intensity"] += f_intensity
            if self.PGN_fragments[f_smiles]["gen"] > attributes["f_gen"]:
                # update parent
                self.PGN_fragments[f_smiles]["gen"] = attributes["f_gen"]
                self.PGN_fragments[f_smiles]["p_idx"] = attributes["f_p_idx"]

        # add new entry
        else:
            adduct_type = FRAG_TYPES_DICT[frag_type].get("Type", "H")
            terminal = FRAG_TYPES_DICT[frag_type].get("Terminal", False)
            self.PGN_fragments[f_smiles] = {
                "idx": self.count,
                "frag_type": frag_type,
                "adduct_type": adduct_type,
                "terminal": terminal,
                "intensity": f_intensity,
                "formula": molecule.formula,
                "mz": molecule.mMass,
                "neutral": molecule.neutral_form,
                "neutral_key": molecule.neutral_InchiKey}
            # rename attributes
            self.PGN_fragments[f_smiles].update(
                (k[2:], v) for k, v in attributes.items() if k[:2] == "f_")
            self.count += 1

        return True

    def fragment_all(self):
        """
        Fragments PGN molecule according to model.

        Returns
        -------
        None.

        """
        # moiety, frag rxn name, fragment 1, fragment 2, parent_only
        frag_again = True
        frag_counts = 0
        fragments = 1
        if not self.quiet:
            print("\n###Fragmentation###\n")
        while frag_again:
            frag_bools = []
            frag_counts += 1
            if not self.quiet:
                print(f"\nFragmenting... round {frag_counts}")
            # iterate through each fragmentation rxn
            for frag_rxn in FRAG_RXN_ORDER:  # dict keys order not preserved
                # retrieve parameters for rxn
                frag_parameters = FRAG_RXN_MAPPING[frag_rxn]
                frag_bools.append(self.fragment_once(
                    frag_rxn, **frag_parameters))
            # clean up after each round
            self.create_full_peaklist_dataframe(True, False, False, False)
            # adjust max int
            if self.adjust_int_each_rd:
                max_intensity = self.peaklist_df_full["intensity"].max()
                self.min_precursor_int = round(
                    self.min_precursor_percent * max_intensity, 3)
                if not self.quiet:
                    print(
                        f"\tAdjusted min precursor int. to {self.min_precursor_int}")
            if not self.quiet:
                added = len(self.peaklist_df_full) - fragments
                print(f"\tAdded {added} fragment(s)")
                fragments += added
            frag_again = any(frag_bools)
        if not self.quiet:
            print(f"\tFinished at round {frag_counts}")
        return frag_counts

    # %%% Tabulation

    def set_dataframe(self, df, simple=False):
        """
        Sets appropriate dataframe variable to 'df'.

        Parameters
        ----------
        df : DataFrame
            Dataframe itself
        simple : Boolean, optional
            If True, use simple dataframe, else full dataframe.
            The default is False.

        Returns
        -------
        df : DataFrame
            Dataframe itself

        """
        if simple:
            self.peaklist_df_simple = df
        else:
            self.peaklist_df_full = df
        return df

    def pick_dataframe(self, overwrite, simple=False):
        """
        Picks appropriate dataframe. If it does not exist, create it.

        Parameters
        ----------
        overwrite : Boolean, optional
            If True, create and overwrite existing dataframe. If False,
            use existing dataframe (unless empty). The default is False.
        simple : Boolean, optional
            If True, use simple dataframe, else full dataframe.
            The default is False.

        Returns
        -------
        df : DataFrame
            Dataframe itself

        """
        if simple:
            if overwrite or not isinstance(self.peaklist_df_simple, pd.DataFrame):
                df = self.create_simple_peaklist_dataframe(overwrite=overwrite)
            else:
                if not self.quiet:
                    print("\tUsing existing simple dataframe...")
                df = self.peaklist_df_simple
        else:
            if overwrite or not isinstance(self.peaklist_df_full, pd.DataFrame):
                df = self.create_full_peaklist_dataframe()
            else:
                if not self.quiet:
                    print("\tUsing existing full dataframe...")
                df = self.peaklist_df_full
        return df

    def create_initial_dataframe(self):
        """
        Creates initial dataframe from self.PGN_fragments. Highest level of detail.

        Returns
        -------
        peaklist_df : DataFrame
            Unprocessed, unsorted dataframe with all columns.

        """
        if not self.quiet:
            print("\t\tCreating dataframe from fragments...")
        peaklist_df = pd.DataFrame.from_dict(
            self.PGN_fragments, orient='index')
        # peaklist_df.dropna(inplace=True)
        peaklist_df.index.rename("SMILES", inplace=True)
        peaklist_df.reset_index(inplace=True)
        peaklist_df.sort_values(by=["mz", "intensity"], inplace=True)
        return peaklist_df

    def create_full_peaklist_dataframe(self,
                                       aggregate_ions=True,
                                       trim_columns=True,
                                       trim_rows=True,
                                       sort_rows=True):
        """
        Creates a full, detailed dataframe from self.PGN_fragments.

        Parameters
        ----------
        aggregate_ions : Boolean, optional
            Determines if ions are aggregated based on mz and neutral form.
            The default is True.
        trim_columns : Boolean, optional
            Determines if unnecessary columns are trimmed.
            The default is True.
        trim_rows : Boolean, optional
            Determines if excess least intense rows are removed.
            The default is True.
        sort_rows : Boolean, optional
            Determines if rows are sorted by mz and intensity.
            The default is True.

        Returns
        -------
        full_df : DataFrame
            Processed, sorted dataframe with relevant columns.

        """
        if not self.quiet:
            creation_string = "".join(
                map(str,
                    (map(int,
                         [aggregate_ions, trim_columns, trim_rows, sort_rows]))))
            print(f"\tCreating full dataframe [{creation_string}]...")
        full_df = self.create_initial_dataframe()
        if aggregate_ions and self.aggregate:
            bef_agg = len(full_df)
            groupby = ["mz", "neutral_key"]
            agg_fn = {col: "first" for col in full_df.columns}
            agg_fn["intensity"] = "sum"
            full_df = full_df.groupby(
                groupby, as_index=False).aggregate(agg_fn)
            # full_df.reset_index(inplace=True)
            after_agg = len(full_df)
            if not self.quiet:
                print(
                    f"\t\tAggregating by mz and neutral form: {bef_agg} --> {after_agg}")
        # rearrange
        if trim_columns:
            columns = ["mz", "intensity", "frag_type", "idx", "gen",
                       "p_idx", "formula", "SMILES"]
            try:
                full_df = full_df[columns]
            except Exception as exc:
                print(full_df.head())
                raise exc
        if trim_rows and len(full_df) > FRAG_MODEL["max_num_ions_save"]:
            if not self.quiet:
                print("\tTruncating least abundant ions to save...")
            full_df = full_df.nlargest(FRAG_MODEL["max_num_ions_save"],
                                       "intensity", keep="all")
        # sort by mz (ascending), intensity (descending)
        if sort_rows:
            full_df.sort_values(by=["mz", "intensity"],
                                ascending=[True, False],
                                inplace=True)
        # debug
        for mz in self.debug_mz:
            selected = full_df["mz"].between(mz-0.05, mz+0.05)
            intensity = full_df[selected]["intensity"].sum()
            if intensity > 0 and not self.quiet:
                print(f"\tDebug mz={mz}; intensity = {intensity}")
        return self.set_dataframe(full_df, simple=False)

    def create_simple_peaklist_dataframe(self, overwrite=False):
        """
        Creates a simple dataframe by aggregating fragments with the same m/z
        from the 'full' dataframe.

        Parameters
        ----------
        overwrite : Boolean, optional
            If True, create and overwrite existing dataframe. If False,
            use existing dataframe (unless empty). The default is False.

        Returns
        -------
        simple_df : DataFrame
            Simplified dataframe with fewer columns and aggregated fragments.

        """
        if not self.quiet:
            print("\tCreating simple dataframe...")
        # Get full dataframe
        peaklist_df = self.pick_dataframe(overwrite, simple=False)
        columns = ["mz", "intensity", "frag_type", "SMILES"]
        try:
            peaklist_df = peaklist_df[columns]
        except Exception as exc:
            print(peaklist_df.columns)
            raise exc
        # agg. by first as peaklist_df is sorted in decreasing intensity
        agg_fn = {"intensity": "sum", "frag_type": "first",
                  "SMILES": "first"}
        bef_agg = len(peaklist_df)
        simple_df = peaklist_df.groupby("mz", as_index=False).aggregate(agg_fn)
        simple_df.sort_values(by="mz", inplace=True)
        after_agg = len(simple_df)
        if not self.quiet:
            print(f"\t\tAggregating by mz: {bef_agg} --> {after_agg}")
        return self.set_dataframe(simple_df, simple=True)

    def filter_and_sum_peaks(self, precursor_mz, adduct_type,
                             omit_greater=MSP_OMIT_IONS_GREATER_THAN_PRECURSOR,
                             omit_mz_threshold=MSP_OMIT_MZ_THRESHOLD,
                             overwrite=False):
        """
        Creates a copy of the simple peaklist.
        Filters peaks outside m/z range.
        Removes peaks that do not suit adduct type.

        Parameters
        ----------
        precursor_mz : float
            m/z of precursor ion.
        adduct_type : string
            Adduct type of precursor ion.
        omit_greater : boolean, optional
            Determines if fragment ions with mz greater than the precursor mz
            are omitted. The default is MSP_OMIT_IONS_GREATER_THAN_PRECURSOR.
        omit_mz_threshold : boolean, optional
            mz threshold for omission. The default is MSP_OMIT_MZ_THRESHOLD.
        overwrite : Boolean, optional
            If True, create and overwrite existing dataframe. If False,
            use existing dataframe (unless empty). The default is False.

        Returns
        -------
        simple_df : DataFrame
            Simplified dataframe with fewer columns and aggregated fragments.

        """
        simple_df = self.pick_dataframe(overwrite, simple=True).copy()
        # omits ions greater than precursor mz
        if omit_greater:
            simple_df = simple_df[
                simple_df["mz"] < precursor_mz+omit_mz_threshold]
        # omits ions not fitting adduct type
        adduct_type = get_adduct_type(adduct_type)
        valid_frag_types = [
            f for f in FRAG_TYPES_DICT
            if FRAG_TYPES_DICT[f].get("Type", "H") == adduct_type]
        simple_df = simple_df[
            simple_df["frag_type"].isin(valid_frag_types)]
        # remove least intense peaks
        if len(simple_df) > FRAG_MODEL["max_num_ions_export"]:
            if not self.quiet:
                print("\tTruncating least abundant ions to export...")
            simple_df = simple_df.nlargest(
                FRAG_MODEL["max_num_ions_export"],
                "intensity",
                keep="all")
        simple_df.sort_values(by="mz", inplace=True)
        return simple_df

    def export_peaklist(self, overwrite=False, simple=False,
                        filename=None, **kwargs):
        """
        Prints peaklist of all fragments except parent.

        Parameters
        ----------
        simple : Boolean, optional
            If True, use simple dataframe, else full dataframe.
            The default is False.
        overwrite : Boolean, optional
            If True, create and overwrite existing dataframe. If False,
            use existing dataframe (unless empty). The default is False.

        Returns
        -------
        Path
            Filename of csv.

        """
        if not self.quiet:
            print(
                f"\nExporting peaklist; overwrite={overwrite}, simple={simple}")
        if simple:
            # export simplified peaklist
            return self.export_simple_peaklist(overwrite=overwrite,
                                               filename=filename,
                                               **kwargs)
        # export full peaklist
        peaklist_df = self.pick_dataframe(overwrite, simple=False)
        if filename is None:
            csv_filename = f"{self.dir}\{TIME_STRING}_{self.PGN_name}_full.csv"
        else:
            csv_filename = f"{self.dir}\{filename}_full.csv"
        try:
            peaklist_df.to_csv(csv_filename, index=False)
        except PermissionError:
            return False
        else:
            return csv_filename

    def export_simple_peaklist(self, overwrite=False,
                               filename=None,
                               adduct_type=None,
                               omit_greater=MSP_OMIT_IONS_GREATER_THAN_PRECURSOR,
                               omit_mz_threshold=MSP_OMIT_MZ_THRESHOLD):
        """
        Prints simple peaklist of all fragments produced by a parent of
        'adduct type' (except parent itself).

        Parameters
        ----------
        adduct_type : string, optional
            Type of adduct: [M+Na]+, [M+2H]2+, etc.
            Defaults to self.adduct_type or "[M+H]+" if None provided.
        omit_greater : boolean, optional
            Determines if fragment ions with mz greater than the precursor mz
            are omitted. The default is MSP_OMIT_IONS_GREATER_THAN_PRECURSOR.
        omit_mz_threshold : boolean, optional
            mz threshold for omission. The default is MSP_OMIT_MZ_THRESHOLD.
        overwrite : Boolean, optional
            If True, create and overwrite existing dataframe. If False,
            use existing dataframe (unless empty). The default is False.

        Returns
        -------
        Path
            Filename of csv.

        """
        if adduct_type is None:
            if self.adduct_type is not None:
                adduct_type = self.adduct_type
            else:
                adduct_type = "[M+H]+"
        simple_df = self.filter_and_sum_peaks(
            precursor_mz=self.PGN.mz[adduct_type],
            adduct_type=adduct_type,
            omit_greater=MSP_OMIT_IONS_GREATER_THAN_PRECURSOR,
            omit_mz_threshold=MSP_OMIT_MZ_THRESHOLD,
            overwrite=overwrite)
        if filename is None:
            csv_filename = \
                f'''{self.dir}\{TIME_STRING}_{self.PGN_name}_{adduct_type}_simple.csv'''
        else:
            csv_filename = \
                f'''{self.dir}\{filename}_{adduct_type}_simple.csv'''
        try:
            simple_df.to_csv(csv_filename, index=False)
        except PermissionError:
            return False
        else:
            return csv_filename

# %%% Testing

    def test_fragmenter(self, molecule=None, **kwargs):
        """
        Tests fragmenter with molecule provided or with a default molecule.

        Parameters
        ----------
        molecule : Molecule, optional
            Molecule to be fragmented. The default is None. If None, uses a
            default molecule.
        **kwargs :
            Passed to auto_process() and export_peaklist()

        Returns
        -------
        None.

        """
        if molecule is None:
            from Molecules import Tan_str_a
            self.load_PGN(Tan_str_a, auto_process=True, **kwargs)
        else:
            self.load_PGN(molecule, auto_process=True, **kwargs)

# %%% ML

    def modify_frag_intensity(self):

        def get_new_intensity(fragment_dict):
            frag_type = fragment_dict["frag_type"]
            frag_occurences = fragment_dict["intensity"]/100
            f_molecule = fragment_dict["molecule"]
            p_molecule = fragment_dict["p_molecule"]
            f_mass = f_molecule.mMass/5000
            p_mass = p_molecule.mMass/5000
            f_peptide = f_molecule.num_peptide_bond/100
            p_peptide = p_molecule.num_peptide_bond/100
            f_dipeptide = f_molecule.num_dipeptide_bond/100
            p_dipeptide = p_molecule.num_dipeptide_bond/100
            f_glycosidic = f_molecule.num_glycosidic_bond/100
            p_glycosidic = p_molecule.num_glycosidic_bond/100
            f_acidic = f_molecule.num_acidic/100
            p_acidic = p_molecule.num_acidic/100
            f_basic = f_molecule.num_acidic/100
            p_basic = p_molecule.num_basic/100
            f_lactoyl = f_molecule.num_lactoyl_bond/100
            p_lactoyl = p_molecule.num_lactoyl_bond/100
            f_lactate = f_molecule.num_lactate/100
            p_lactate = p_molecule.num_lactate/100
            params = np.array([
                f_mass, p_mass-f_mass,
                f_acidic, p_acidic-f_acidic,
                f_basic, p_basic-f_basic,
                f_peptide, p_peptide-f_peptide,
                f_dipeptide, p_dipeptide-f_dipeptide,
                f_glycosidic, p_glycosidic-f_glycosidic,
                f_lactoyl, p_lactoyl-f_lactoyl,
                f_lactate, p_lactate-f_lactate,
                frag_occurences])  # 17
            cat = {"f": frag_type, "adduct": self.adduct_type}
            f_int = self.frag_network.calculate_value(params, cat=cat)
            return f_int

        for frag, frag_dict in self.PGN_fragments.items():
            if frag_dict["frag_type"] != "Parent":
                new_int = get_new_intensity(frag_dict)
                frag_dict["intensity"] = new_int

        self.create_full_peaklist_dataframe()

# %% Fragmenter Testing


Ec940 = Molecule(
    "CC(=O)NC1C(OC2C(CO)OC(O)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O")
Sa1251 = Molecule(
    "CC(=O)NC1C(OC2C(CO)OC(O)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCCNC(=O)CNC(=O)CNC(=O)CNC(=O)CNC(=O)CN)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(N)=O)OC(CO)C(O)C1O")
Efs1108 = Molecule(
    "CC(=O)NC1C(OC2C(CO)OC(O)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCCNC(=O)C(C)NC(=O)C(C)N)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(N)=O)OC(CO)C(O)C1O")

# %% MSPMaker


class MSPMaker():

    def __init__(self, name="test", ionmode="Positive",
                 output_adducts=OUTPUT_ADDUCTS,
                 omit_greater=MSP_OMIT_IONS_GREATER_THAN_PRECURSOR,
                 omit_mz_threshold=MSP_OMIT_MZ_THRESHOLD,
                 precursor_mz_min=MSP_PRECURSOR_MZ_MIN,
                 precursor_mz_max=MSP_PRECURSOR_MZ_MAX
                 ):
        """


        Parameters
        ----------
        name : string, optional
            Name of MSPMaker (used in filenames). The default is "test".
        ionmode : string, optional
            Ion mode. Currently serves no major function.
            The default is "Positive".
        output_adducts : list, optional
            List of output adducts. The default is OUTPUT_ADDUCTS.
        omit_greater : boolean, optional
            Determines if fragment ions with mz greater than the precursor mz
            are omitted. The default is MSP_OMIT_IONS_GREATER_THAN_PRECURSOR.
        omit_mz_threshold : boolean, optional
            mz threshold for omission. The default is MSP_OMIT_MZ_THRESHOLD.
        precursor_mz_min : float, optional
            Min. mz of precursor. The default is MSP_PRECURSOR_MZ_MIN.
        precursor_mz_max : float, optional
            Max. mz of precursor. The default is MSP_PRECURSOR_MZ_MAX.

        Returns
        -------
        None.

        """
        self.name = None
        self.dir = None
        self.refresh_name(name)
        self.ionmode = ionmode
        self.counter = Counter()
        self.peaklist_df_data = {}  # InchiKey: [mz:{int/smiles}]
        self.fragmenter = Fragmenter(self.name, frag_mode="formula",
                                     quiet=True)
        self.omit_greater = omit_greater
        self.omit_mz_threshold = omit_mz_threshold
        self.precursor_mz_min = precursor_mz_min
        self.precursor_mz_max = precursor_mz_max
        self.output_adducts = output_adducts
        print(
            f'''\nCreating msp for adducts: \n{self.output_adducts}
            Omitting ions > precursor: {self.omit_greater}
            with mz threshold: {self.omit_mz_threshold}
            Precursor mz: {self.precursor_mz_min} - {self.precursor_mz_max}''')
        self.msp_filename = None

    def refresh_name(self, new_name, from_dir=False):
        """
        Refreshes the name of MSPMaker.

        Parameters
        ----------
        new_name : String or Path
            New name/path.
        from_dir : bool, optional
            Indicates if name is from path. The default is False.

        Returns
        -------
        None.

        """
        if from_dir:
            # convert dir to name
            path = Path(new_name)
            new_stem = path.stem
            new_name = re.sub("[0-9]{12}_", "", new_stem)  # remove date
        if self.name != new_name:
            self.name = new_name
            print(f"\nMSPMaker {self.name} is created.")
            self.dir = f"output/{self.name}/msp"
            make_dir(f"output/{self.name}", self.dir)

    # %%% Processing

    def create_mspfile(self, overwrite=False):
        """
        Creates an empty mspfile.

        Parameters
        ----------
        overwrite : boolean, optional
            If True, overwrites existing msp file. The default is False.

        Returns
        -------
        self.msp_filename: Path
            msp path.

        """
        if self.msp_filename is None or overwrite:
            self.msp_filename = f"{self.dir}/{TIME_STRING}_{self.name}.msp"
        file_exists = os.path.exists(self.msp_filename)
        if file_exists and not overwrite:
            print(f"\nAppending to existing mspfile:\t{self.msp_filename}")
        else:
            with open(self.msp_filename, "w", encoding="ASCII") as f:
                f.close()
            print(f"\nCreated mspfile:\t{self.msp_filename}")
        return self.msp_filename

    def process_bool(self, k, idx, k_allowed, idx_allowed):
        """
        Checks if compound data of index 'idx' containing 'k-mers' should
        be processed.

        Parameters
        ----------
        k : integer
            Type of -mer (1-mer, 2-mer) found in compound/peaklist data.
        idx : integer
            Index of compound/peaklist data.
        k_allowed : list
            List of allowed 'k's
        idx_allowed : list
            List of allowed 'idx's

        Returns
        -------
        boolean
            Returns True if compound data should be processed.

        """
        k_bool = (k_allowed == "all" or k in k_allowed)
        idx_bool = (idx_allowed == "all" or idx in idx_allowed)
        return k_bool and idx_bool

    def process_compound_data(self, k, idx):
        """
        Processes current compound data and generates peaklist data from it.

        Parameters
        ----------
        k : integer
            Type of -mer (1-mer, 2-mer) found in compound/peaklist data.
        idx : integer
            Index of compound/peaklist data.

        Returns
        -------
        num_cpds : integer
            Number of compounds processed.

        """
        num_cpds = len(self.compound_data)
        print(f"\tFragmenting batch {idx}:\t{num_cpds} {k}-mers")
        self.counter.set_counter(num_cpds, f"{k}-mers")
        frag = self.fragmenter
        # possible parallelization
        for InchiKey in self.compound_data:
            PGN = self.compound_data[InchiKey]
            fragment_data = cache_spectral_data(frag, PGN, InchiKey)
            self.peaklist_df_data[InchiKey] = fragment_data
            self.counter.update_counter()
        return num_cpds

    def import_PGN_as_pickle(self, pickle_filename, inherit_name=True,
                             k_allowed="all", idx_allowed="all",
                             overwrite=False, msp_filename=None,
                             export="all"):
        """
        Imports a pickle file generated by Generator to create a mspfile.

        Parameters
        ----------
        pickle_filename : Path
            Path of pickle file.
        inherit_name : bool, optional
            Determines if name is inherited from pickle file.
            The default is True.
        k_allowed : List of integers or "all", optional
            Generate msp for k-mers (k>0).
            If "all", generates msp for all compounds. The default is "all".
        idx_allowed : List of integers or "all", optional
            Generate msp for batches. (idx>0).
            If "all", generates msp for all batches. The default is "all".
        overwrite : bool, optional
            Determines if current msp (if it exists) is overwritten.
            The default is False.
        msp_filename : Path, optional
            msp_file to be extended (if overwrite is False).
        export : List of strings, optional
            "pickle" : Pickled compound data (fragments)
            "msp" : mspfile
            "peaklist" : Detailed compound data in xlsx.
            "ionlist" : Consolidated list of ions for all compounds in xlsx.
            "all" : All outputs.
            The default is "all".

        Returns
        -------
        None.

        """
        # Check export logic
        valid_exports = ["pickle", "msp", "peaklist", "ionlist"]
        if export == "all":
            export = valid_exports
        if type(export) == str:
            export == [export]
        export = [x for x in export if x in valid_exports]
        if len(export) == 0:
            print("\n!!!Nothing exported!!!\n")
        elif "ionlist" in export or "peaklist" in export:
            if "pickle" not in export:
                print("'Pickle' files must be exported for 'ionlist'/'peaklist'.")
                export.append("pickle")
        # Create and name files
        self.counter.disable_sleep()
        if inherit_name:
            self.refresh_name(pickle_filename, from_dir=True)
            self.fragmenter.refresh_name(pickle_filename, from_dir=True)
            overwrite = True
        if msp_filename is not None:
            self.msp_filename = msp_filename
        if "msp" in export:
            self.create_mspfile(overwrite=overwrite)
            self.write_pkl_path(pickle_filename)
        peaklists_pkl_list = []

        print(f"\nOpening\t{pickle_filename}")
        with open(pickle_filename, "rb") as f:
            k = 1
            idx = 1
            while True:
                try:
                    self.compound_data = pickle.load(f)
                    curr_batch_size = len(self.compound_data)
                    print(
                        f"\nLoaded batch {idx}:\t {curr_batch_size} {k}-mers")
                    if self.process_bool(k, idx, k_allowed, idx_allowed):
                        # fragment
                        self.process_compound_data(k, idx)
                        # export
                        if "pickle" in export:
                            peaklists_pickle = self.export_peaklists_as_pickle(
                                idx)
                            peaklists_pkl_list.append(peaklists_pickle)
                        if "msp" in export:
                            self.write_msp()
                        # clear
                        self.peaklist_df_data = {}
                        self.compound_data = {}
                        print(f"\tBatch {idx} completed.")
                    else:
                        print(f"\tBatch {idx} skipped.")
                        if "pickle" in export:
                            # Insert null filename into pickle list when skipped
                            peaklists_pkl_list.append(None)
                    if curr_batch_size < BATCH_SIZE:
                        k += 1
                    idx += 1
                except EOFError:  # end of file
                    break
        print(f"\nClosed {pickle_filename}")
        if "ionlist" in export:
            self.export_consolidated_ion_data(peaklists_pkl_list,
                                              inherit_name=False)
        if "peaklist" in export:
            self.export_consolidated_peaklists(pickle_filename,
                                               peaklists_pkl_list,
                                               inherit_name=False)
        self.counter.show_toast("PSN_MS2: MSPMaker Job complete",
                                f"{self.name} completed.",
                                10)
        self.counter.enable_sleep()
        return pickle_filename, peaklists_pkl_list

    # %%% Exporting

    def write_pkl_path(self, cpd_pkl_filename):
        """
        Writes a short text file to reference the path of the pickle file in
        which compound data is taken from.

        Parameters
        ----------
        cpd_pkl_filename : Path
            Path of pickle file containing compound data.

        Returns
        -------
        file_path : Path
            Reference file path.

        """
        file_path = f"{self.dir}/{TIME_STRING}_{self.name}.txt"
        with open(file_path, "w", encoding="ASCII") as f:
            f.write(str(cpd_pkl_filename))
            f.close()
        print(f"Source pickle path saved at:\t{file_path}")
        return file_path

    def write_msp(self):
        """
        Writes mspfile.

        Returns
        -------
        bool
            True when completed without errors.

        """
        frag = self.fragmenter
        with open(self.msp_filename, "a", encoding="ASCII") as f:
            def write_peak(row):
                # format info for each peak
                mz = mztxt(row["mz"])
                rel_int = inttxt(row["intensity"])
                frag_type = row["frag_type"]
                f.write(f"{mz}\t{rel_int}\t{frag_type}\n")
            for InchiKey in self.compound_data:
                fragment_data = self.peaklist_df_data[InchiKey]
                PGN = self.compound_data[InchiKey]
                frag.load_PGN(PGN, InchiKey)
                frag.set_dataframe(fragment_data, simple=False)
                frag.create_simple_peaklist_dataframe()
                for adduct in self.output_adducts:
                    # exit loop if adduct absent
                    if adduct not in PGN.mz:
                        continue
                    precursor_mz = PGN.mz[adduct]
                    # exit loop if precursor mz above/below mz range
                    if self.precursor_mz_min > precursor_mz or\
                            self.precursor_mz_max < precursor_mz:
                        continue
                    # omit product ions > precursor
                    peaklist_data = frag.filter_and_sum_peaks(
                        precursor_mz, adduct,
                        omit_greater=self.omit_greater,
                        omit_mz_threshold=self.omit_mz_threshold)
                    # exit loop if no product ions
                    if len(peaklist_data) == 0:
                        continue
                    # convert characters to ASCII
                    f.write(f"NAME: {PGN.name}\n")
                    f.write(f"PRECURSORMZ: {mztxt(precursor_mz)}\n")
                    f.write(f"PRECURSORTYPE: {adduct}\n")
                    f.write(f"SMILES: {PGN.smiles}\n")
                    f.write(f"INCHIKEY: {PGN.InchiKey}\n")
                    f.write(f"FORMULA: {PGN.formula}\n")
                    f.write("RETENTIONTIME: 0\n")
                    f.write(f"IONMODE: {self.ionmode}\n")
                    f.write(f"COMPOUNDCLASS: {PGN.ontology}\n")
                    f.write("Comment: theoretical MS2\n")
                    f.write(f"Num Peaks: {len(peaklist_data)}\n")
                    peaklist_data.apply(write_peak, axis='columns')
                    f.write("\n")
        print("\tWrote to mspfile.")
        return True

    def write_peaklist_batched(self, writer, sheet_name):
        """
        Writes peaklist with 'writer' in 'sheet_name'. Helper function for
        'export_consolidated_peaklists'.

        Parameters
        ----------
        writer : ExcelWriter
            Writer.
        sheet_name : string
            Name of sheet.

        Returns
        -------
        None

        """
        startrow = 2
        for InchiKey, PGN in self.compound_data.items():
            fragment_data = self.peaklist_df_data[InchiKey]
            fragment_data.to_excel(
                writer, sheet_name=sheet_name,
                startrow=startrow, index=False)
            # write identifiers, cell_index starts from 1,1
            writer.sheets[sheet_name].cell(startrow-1, 1, "Parent Name")
            writer.sheets[sheet_name].cell(startrow-1, 2, "Parent Synonym")
            writer.sheets[sheet_name].cell(startrow-1, 3, "Parent Mass")
            writer.sheets[sheet_name].cell(startrow-1, 4, "Parent Formula")
            writer.sheets[sheet_name].cell(startrow-1, 5, "Parent InChiKey")
            writer.sheets[sheet_name].cell(startrow-1, 6, "Parent SMILES")
            writer.sheets[sheet_name].cell(startrow, 1, PGN.name)
            writer.sheets[sheet_name].cell(startrow, 2, PGN.synonym)
            writer.sheets[sheet_name].cell(startrow, 3, PGN.mMass)
            writer.sheets[sheet_name].cell(startrow, 4, PGN.formula)
            writer.sheets[sheet_name].cell(startrow, 5, InchiKey)
            writer.sheets[sheet_name].cell(startrow, 6, PGN.smiles)
            startrow += fragment_data.shape[0]+4

    def export_consolidated_peaklists(self, cpd_pkl_filename,
                                      peaklists_pkl_list,
                                      inherit_name=True, inherit_index=0):
        """
        Consolidates peaklists from a list of peaklist pickle files.

        Parameters
        ----------
        cpd_pkl_filename : Path
            Path of pickle file containing compound data.
        peaklists_pkl_list : List
            Paths of saved peaklist data (pickles) in a list.
        inherit_name : boolean, optional
            If True, MSPMaker's name will be inherited from pickle file.
            The default is True.
        inherit_index : integer, optional
            If inherit_name, name willbe inherited from pickle file with this
            index. The default is 0.

        Returns
        -------
        None.

        """
        print("\nConsolidating peaklists...")
        if inherit_name:
            self.refresh_name(
                peaklists_pkl_list[inherit_index], from_dir=True)
        xlsx_filename = f"output/{self.name}/peaklists/{TIME_STRING}_{self.name}_spectraldata.xlsx"
        print(f"\tCreated peaklist:\t{xlsx_filename}")
        writer = pd.ExcelWriter(xlsx_filename)
        # Open compound data as pickle
        with open(cpd_pkl_filename, "rb") as f:
            k = 1
            idx = 1
            while True:
                try:
                    # Load compound data
                    self.compound_data = pickle.load(f)
                    curr_batch_size = len(self.compound_data)
                    print(
                        f"\nLoaded batch {idx}:\t {curr_batch_size} {k}-mers")
                    sheet_name = f"Batch{idx}"
                    # check for peaklist data
                    peaklists_pkl_file = peaklists_pkl_list[idx-1]
                    if peaklists_pkl_file is None:
                        # skip placeholder (None)
                        print(f"\t\tBatch {idx} skipped")
                    else:
                        self.import_peaklist_data_as_pickle(peaklists_pkl_file)
                        self.write_peaklist_batched(writer, sheet_name)
                        print(f"\tExported peaklists for batch {idx}")
                        # clear
                        self.peaklist_df_data = {}
                        self.compound_data = {}
                    if curr_batch_size < BATCH_SIZE:
                        k += 1
                    idx += 1
                except EOFError:  # end of file
                    break
        writer.close()
        print(f"\nClosed {cpd_pkl_filename}")

    def export_consolidated_ion_data(self, peaklists_pkl_list,
                                     inherit_name=True, inherit_index=0):
        """
        Consolidates fragments from a list of peaklist pickle files. Fragments
        with the same SMILES are removed.

        Parameters
        ----------
        peaklists_pkl_list : List
            Paths of saved peaklist data (pickles) in a list.
        inherit_name : boolean, optional
            If True, MSPMaker's name will be inherited from pickle file.
            The default is True.
        inherit_index : integer, optional
            If inherit_name, name willbe inherited from pickle file with this
            index. The default is 0.

        Returns
        -------
        xlsx_filename: Path
            Consolidated ion data path.

        """
        print("\nConsolidating ion data...")
        if inherit_name:
            self.refresh_name(
                peaklists_pkl_list[inherit_index], from_dir=True)
        column_names = ["mz", "formula", "SMILES"]
        consolidated_df = pd.DataFrame(columns=column_names)
        # go through each file and consolidate
        for idx, file in enumerate(peaklists_pkl_list):
            if file == None:
                # skip placeholder (None)
                print(f"\t\tBatch {idx+1} skipped")
            else:
                self.import_peaklist_data_as_pickle(file)
                batch_df = pd.concat(self.peaklist_df_data.values())
                batch_df = batch_df[column_names]
                batch_df.drop_duplicates(subset="SMILES", inplace=True)
                consolidated_df = pd.concat(
                    [consolidated_df, batch_df], ignore_index=True)
                print(f"\t\tBatch {idx+1} completed")
        xlsx_filename = f"output/{self.name}/peaklists/{TIME_STRING}_{self.name}_iondata.xlsx"
        consolidated_df.drop_duplicates(subset="SMILES", inplace=True)
        # remove parents
        consolidated_df = consolidated_df[consolidated_df["mz"] != 0]
        consolidated_df.sort_values(by=["mz"], inplace=True)
        consolidated_df.to_excel(xlsx_filename, index=False,
                                 sheet_name="IonData")
        print(f"\tIon data saved as:\t{xlsx_filename}")
        return xlsx_filename

    def export_peaklists_as_pickle(self, idx):
        """
        Exports peaklist data as pickle file.

        Parameters
        ----------
        idx : integer
            Index of compound/peaklist data.

        Returns
        -------
        pickle_filename : Path
            Pickle file path.

        """
        pickle_filename = f"{self.dir}/{TIME_STRING}_{self.name}_{idx}.pickle"
        with open(pickle_filename, "wb") as f:
            pickle.dump(self.peaklist_df_data, f)
            f.close()
        print(f"\tPeaklist data saved as:\t{pickle_filename}")
        return pickle_filename

    def import_peaklist_data_as_pickle(self, pickle_filename):
        """
        Imports peaklist data as pickle file.

        Parameters
        ----------
        pickle_filename : Path
            Pickle file path.

        Returns
        -------
        None.

        """
        with open(pickle_filename, "rb") as f:
            self.peaklist_df_data = pickle.load(f)
        print(f"\tPeaklist data imported from: \t{pickle_filename}")
