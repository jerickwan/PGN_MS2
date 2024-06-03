'''
AA, Peptide, Glycan Molecular Classes
env = chem
'''

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdAbbreviations
from rdkit.Chem import Draw  # drawing functions
from rdkit.Chem.Descriptors import MolLogP
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdmolops import GetFormalCharge
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt, CalcTPSA, CalcNumHeavyAtoms

import itertools
import numpy as np
import pandas as pd

import copy

from base.Reactions import PEPTIDE_RXN_TYPES, GLYCAN_RXN_TYPES
from base.Reactions import PEPTIDE_PATT_TYPES, MISC_PATT_TYPES
from base.Reactions import MISC_RXN_UNCHARGE_CO_A, MISC_RXN_UNCHARGE_M
from base.Exceptions import NoReactionError, TooManyPdtsError
from base.Common import BOND_CHAR, BRANCH_CHAR, shBRANCH_CHAR,\
    PRECISION_MASS, MEMORY, OUTPUT_ADDUCTS
from base.Common import flatten


# %% Globals

DATABASE = "data/PGN.xlsx"
AMINO_ACID_DF = None
AMINO_ACID_DB = {}
AMINO_ACID_AMIDATED_LST = ["mDAP(NH2)", "γ-isoGln", "β-isoAsn", "Asn"]
GLYCAN_DF = None
GLYCAN_DB = {}
UNCHARGE_RXNS = (MISC_RXN_UNCHARGE_CO_A, MISC_RXN_UNCHARGE_M)
if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    AUTOCALC_MZ = True
else:
    AUTOCALC_MZ = False

# %% Handling Mol Objects


def create_mol_cache_key(molecule_lst):
    """
    Creates a string from a list of molecules.

    Parameters
    ----------
    molecule_lst : list
        List of Molecule objects.

    Returns
    -------
    string
        Unique string that corresponds to the list of molecules.

    """
    return ",".join(x.InchiKey for x in molecule_lst)


def get_cation_mz(mMass, atom, num):
    """
    Calculates m/z for cationic adduct [M+(num)(atom)] for given monoisotopic mass (mMass).

    Parameters
    ----------
    mMass : float
        Monoisotopic molecular mass.
    atom : letter
        Atom added as adduct.
    num : int
        Number of atoms added.

    Returns
    -------
    float
        m/z for cationic adduct.

    """
    atoms = {"H": 1.00783, "Na": 22.98977, "K": 38.96371}
    e_mass = 0.00054858
    mz = (mMass + num*atoms[atom] - num*e_mass)/num
    return round(mz, PRECISION_MASS)


def create_peptide_rxn(self, other,
                       COOH_side_chain, NH2_side_chain):
    """
    Creates custom peptide rxn for self (COOH terminus) and other (NH2 terminus),
    taking into account side chains (COOH or NH2).

    Parameters
    ----------
    other : AminoAcid, Peptide
    COOH_side_chain : boolean
    NH2_side_chain : boolean

    Returns
    -------
    Chem.Reaction
        Reaction for the two amino acids.

    """
    # Select Template
    COOH_amino_acid = isinstance(self, AminoAcid)
    NH2_amino_acid = isinstance(other, AminoAcid)
    if COOH_amino_acid and NH2_amino_acid:
        smarts = PEPTIDE_RXN_TYPES["Template Eupeptide"]
    elif COOH_amino_acid:
        smarts = PEPTIDE_RXN_TYPES["Template Eupeptide COOH"]
    else:
        smarts = PEPTIDE_RXN_TYPES["Template Eupeptide NH2"]
    # Modify COOH
    if COOH_amino_acid:
        if COOH_side_chain and self.is_smarts_valid("COOH",side_chain=True):
            COOH_smarts = self.side_chain_smarts["COOH"]
        else:
            COOH_smarts = self.eupeptide_smarts["COOH"]
        if type(COOH_smarts) != str:
            raise NoReactionError(self, other, "peptide")
        else:
            smarts = smarts.replace("*COOH*", COOH_smarts)
    # Modify NH2
    if NH2_amino_acid:
        if NH2_side_chain and other.is_smarts_valid("NH2",side_chain=True):
            NH2_smarts = other.side_chain_smarts["NH2"]
        else:
            NH2_smarts = other.eupeptide_smarts["NH2"]
        if type(NH2_smarts) != str:
            raise NoReactionError(self, other, "peptide")
        else:
            smarts = smarts.replace("*NH2*", NH2_smarts)
    return AllChem.ReactionFromSmarts(smarts)


def count_peptide_bonds(mol, bond_type="Peptide"):
    """
    Returns number of peptide bonds for peptide (mol).

    Parameters
    ----------
    mol : Chem.Mol
        Mol object.
    bond_type : string, optional
        Type of peptide bond. The default is "peptide".

    Returns
    -------
    num_peptide_bonds : int
        Number of peptide bonds that match pattern.

    """
    substruct_patt = PEPTIDE_PATT_TYPES[bond_type]
    num_peptide_bonds = len(mol.GetSubstructMatches(substruct_patt))
    return num_peptide_bonds


def count_peptide_linearity(mol):
    """
    Returns linearity of peptide (mol); i.e. how similar is this peptide to a
    typical peptide chain of the same length.

    Parameters
    ----------
    mol : Chem.Mol
        Mol object.

    Returns
    -------
    int
        Measure of peptide linearity.

    """

    return sum([count_peptide_bonds(mol, bond_type="Alpha C Eupeptide"),
                count_peptide_bonds(mol, bond_type="Alpha N Eupeptide"),
                count_peptide_bonds(mol, bond_type="Dipeptide")*2,
                count_peptide_bonds(mol, bond_type="isoGlx Dipeptide")*2,
                count_peptide_bonds(mol, bond_type="Tripeptide")*3,
                count_peptide_bonds(mol, bond_type="Lax Tripeptide")*3])


def uncharge_molecule(mol):
    """
    Uncharges molecule by removing H+ (with re-arrangement if needed)

    Parameters
    ----------
    mol : Chem.Mol
        Mol object.

    Returns
    -------
    mol : Chem.Mol
        Returns uncharged Mol object if successful else False.
    """
    charge = GetFormalCharge(mol)
    if charge != 0:
        uncharger = rdMolStandardize.Uncharger()
        # remove H+
        mol = uncharger.uncharge(mol)
        # check for success
        charge = GetFormalCharge(mol)
        if charge != 0:
            for uncharge_rxn in UNCHARGE_RXNS:
                if uncharge_rxn.RunReactantInPlace(mol):
                    Chem.SanitizeMol(mol)
                    break
            charge = GetFormalCharge(mol)
            # all uncharge_rxns failed
            if charge != 0:
                return False
    return mol


def pick_correct_pdt(pdts, pick_min):
    """
    Picks the maximum or minimum linearity product from a list of products (pdts).

    Parameters
    ----------
    pdts : list
        List of Chem.Mol objects
    pick_min : boolean
        If True, picks the product with the minimum linearity else maximum.

    Returns
    -------
    mol : Chem.Mol
        Returns Mol object if only one product is picked else returns False.

    """
    pl = [count_peptide_linearity(Pep.mol) for Pep in pdts]
    max_pl = max(pl)
    min_pl = min(pl)

    if pick_min:
        min_pl_pdts = [
            x for x in pdts if count_peptide_linearity(x.mol) == min_pl]
        if len(min_pl_pdts) > 1:
            return False
        else:
            return min_pl_pdts[0]

    else:
        max_pl_pdts = [
            x for x in pdts if count_peptide_linearity(x.mol) == max_pl]
        if len(max_pl_pdts) > 1:
            return False
        else:
            return max_pl_pdts[0]

# %%% Cached


def run_mono_pdt_rxn(rxn, cpd1, cpd2):
    """
    Runs one pdt reaction (rxn) between cpd1 and cpd2 (as mol),
    returns list of sanitized smiles.

    Parameters
    ----------
    rxn : Chem.Reaction
        Reaction object
    cpd1 : Molecule
        Compound 1
    cpd2 : Molecule
        Compound 2

    Returns
    -------
    pdts_smiles : List
        List of sanitized SMILES corresponding to product.
    """
    try:
        pdts_mols = list(itertools.chain.from_iterable(
            rxn.RunReactants((cpd1, cpd2))))
    except ValueError as exc: #ChemicalParserException
        #reaction has duplicate labels
        smarts = AllChem.ReactionToSmarts(rxn)
        smiles1 = Chem.MolToSmiles(cpd1)
        smiles2 = Chem.MolToSmiles(cpd2)
        print(smiles1, smiles2)
        print(smarts)
        raise exc

    for i, pdt in enumerate(pdts_mols):
        Chem.SanitizeMol(pdt)  # sanitize molecule
        # img = Draw.MolToImage(pdt, size=(500, 500))
        # img.save(f"img/last_rxn_{i}.png")

    # sanity check on molecule
    pdts_smiles = [Chem.MolToSmiles(mol, canonical=True)
                   for mol in pdts_mols]
    return pdts_smiles


def run_two_pdt_rxn(rxn, cpd1, cpd2):
    """
    Runs two pdt reaction (rxn) between cpd1 and cpd2 (as mol),
    returns list of sanitized smiles.

    Parameters
    ----------
    rxn : Chem.Reaction
        Reaction object
    cpd1 : Molecule
        Compound 1
    cpd2 : Molecule
        Compound 2

    Returns
    -------
    pdts1_smiles : List
        List of sanitized SMILES corresponding to product 1.
    pdts2_smiles : List
        List of sanitized SMILES corresponding to product 2.

    """
    try:
        pdts_mols = rxn.RunReactants((cpd1, cpd2))
    except ValueError as exc: #ChemicalParserException
        #reaction has duplicate labels
        smarts = AllChem.ReactionToSmarts(rxn)
        smiles1 = Chem.MolToSmiles(cpd1)
        smiles2 = Chem.MolToSmiles(cpd2)
        print(smiles1, smiles2)
        print(smarts)
        raise exc
    pdts1_smiles = []
    pdts2_smiles = []

    for i, (pdt1, pdt2) in enumerate(pdts_mols):
        Chem.SanitizeMol(pdt1)
        Chem.SanitizeMol(pdt2)
        # img = Draw.MolsToImage(mols=[pdt1, pdt2], subImgSize=(500, 500))
        # img.save(f"img/last_rxn_{i}.png")
        # sanity check on molecule
        pdts1_smiles.append(Chem.MolToSmiles(pdt1, canonical=True))
        # sanity check on molecule
        pdts2_smiles.append(Chem.MolToSmiles(pdt2, canonical=True))

    return pdts1_smiles, pdts2_smiles


def run_fragment_rxn(frag_rxn, mMass, rxn, cpd):
    """
    Runs fragmentation reaction of parent molecule and returns list of
    sanitized smiles [pdt1,pdt2].

    Parameters
    ----------
    frag_rxn : string
        Type of fragmentation
    mMass : float
        Monoisotopic mass of parent molecule.
    rxn : Chem.Reaction
        Reaction object for fragmentation event
    cpd : Molecule
        Parent molecule undergoing fragmentation

    Returns
    -------
    pdts1_smiles : List
        List of sanitized SMILES corresponding to product 1.
    pdts2_smiles : List
        List of sanitized SMILES corresponding to product 2.

    """
    pdts_mols = rxn.RunReactants((cpd,))
    pdts1_smiles = []
    pdts2_smiles = []
    for i, (pdt1, pdt2) in enumerate(pdts_mols):
        Chem.SanitizeMol(pdt1)
        Chem.SanitizeMol(pdt2)
        pdts1_smiles.append(Chem.MolToSmiles(pdt1, canonical=True))
        pdts2_smiles.append(Chem.MolToSmiles(pdt2, canonical=True))
    return pdts1_smiles, pdts2_smiles


@MEMORY.cache
def run_polymerisation_rxn(polymerisation_name, P1_bond, P2_bond, P1, P2):
    """
    Runs polymerisation reaction between two peptidoglycan molecules (P1,P2)

    Parameters
    ----------
    polymerisation_name : string
        Name of polymerisation.
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
    P1 : Peptidoglycan
        Peptidoglycan molecule 1
    P2 : Peptidoglycan
        Peptidoglycan molecule 2

    Raises
    ------
    ValueError
        Incorrect bonding provided.
    NoReactionError
        No reaction occurs.
    TooManyPdtsError
        Too many similar products created.

    Returns
    -------
    pdt : Peptidoglycan
         Peptidoglycan polymer of P1 and P2

    """

    P1_AApos, P1_chain, P1_rxntype, P1_grp = P1_bond
    P2_AApos, _, P2_rxntype, P2_grp = P2_bond

    if P1_AApos == 0 and P2_AApos == 0:
        # Glycan polymerisation
        results = run_glycan_polymerisation(
            polymerisation_name, P1_bond, P2_bond, P1, P2)
    else:
        # Peptide polymerisation
        results = run_peptide_polymerisation(
            polymerisation_name, P1_bond, P2_bond, P1, P2)

    pdts_smiles, name_ext, rxn, NH2_side_chain, COOH_side_chain = results
    if len(pdts_smiles) == 0:
        raise NoReactionError(P1, P2, name_ext, rxn)
    # combine modifications
    combined_mods = {
        mod: P1.modifications[mod] or P2.modifications[mod]
        for mod in Peptidoglycan.modifications}
    pdts = [Peptidoglycan(x,
                          glycan=P1.glycan,
                          peptide=P1.peptide,
                          modifications=combined_mods,
                          polymerisation_type=polymerisation_name,
                          linked_PGN=P2) for x in set(pdts_smiles)]
    if len(pdts) == 1:
        pdt = pdts[0]
    else:
        pick_min = NH2_side_chain and COOH_side_chain
        pdt = pick_correct_pdt(pdts, pick_min)
        if pdt == False:
            raise TooManyPdtsError(
                P1, P2, polymerisation_name, rxn, pdts_smiles)
    return pdt


def run_peptide_polymerisation(polymerisation_name, P1_bond, P2_bond, P1, P2):
    """
    Runs polymerisation reaction between two peptidoglycan molecules (P1,P2)
    through the peptide moiety.

    Parameters
    ----------
    polymerisation_name : string
        Name of polymerisation.
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
    P1 : Peptidoglycan
        Peptidoglycan molecule 1
    P2 : Peptidoglycan
        Peptidoglycan molecule 2

    Raises
    ------
    ValueError
        Incorrect bonding provided.

    Returns
    -------
    pdts_smiles : list
        List of sanitized SMILES corresponding to product.
    name_ext : string
        String with an extended description of the polymerisation.
        Used when raising NoReactionError.
    rxn : Chem.Reaction
        Reaction object for polymerisation.
    NH2_side_chain : boolean
    COOH_side_chain : boolean

    """
    # Peptide polymerisation
    P1_AApos, P1_chain, P1_rxntype, P1_grp = P1_bond
    P2_AApos, _, P2_rxntype, P2_grp = P2_bond
    # Pick correct amino acids on P1/P2 for bonding
    if P1_chain == "main":
        P1_AA = P1.peptide_code[P1_AApos-1]
    elif P1_chain == "bridge":
        _, peptide_code = P1.peptide.side_chain_dict[P1_AApos]
        P1_AA = peptide_code[-1]  # outermost AA
    # P2 is always on main chain
    P2_AA = P2.peptide_code[P2_AApos-1]
    # Pick correct orientation of bonding
    if P1_grp == "COOH" and P2_grp == "NH2":
        AA_C = AMINO_ACID_DB[P1_AA]
        AA_N = AMINO_ACID_DB[P2_AA]
        COOH_side_chain = P1_rxntype == "side"
        NH2_side_chain = P2_rxntype == "side"
        P_C = P1
        P_N = P2
    elif P2_grp == "COOH" and P1_grp == "NH2":
        AA_C = AMINO_ACID_DB[P2_AA]
        AA_N = AMINO_ACID_DB[P1_AA]
        COOH_side_chain = P2_rxntype == "side"
        NH2_side_chain = P1_rxntype == "side"
        P_C = P2
        P_N = P1
    else:
        raise ValueError
    rxn = create_peptide_rxn(
        AA_C, AA_N, COOH_side_chain, NH2_side_chain)
    pdts_smiles = run_mono_pdt_rxn(rxn, P_C.mol, P_N.mol)
    name_ext = f"{polymerisation_name}_{AA_C}_{COOH_side_chain}_{AA_N}_{NH2_side_chain}"
    return pdts_smiles, name_ext, rxn, NH2_side_chain, COOH_side_chain


def run_glycan_polymerisation(polymerisation_name, P1_bond, P2_bond, P1, P2):
    """
    Runs polymerisation reaction between two peptidoglycan molecules (P1,P2)
    through the glycan moiety.

    Parameters
    ----------
    polymerisation_name : string
        Name of polymerisation.
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
    P1 : Peptidoglycan
        Peptidoglycan molecule 1
    P2 : Peptidoglycan
        Peptidoglycan molecule 2

    Raises
    ------
    ValueError
        Incorrect bonding provided.

    Returns
    -------
    pdts_smiles : list
        List of sanitized SMILES corresponding to product.
    name_ext : string
        String with an extended description of the polymerisation.
        Used when raising NoReactionError.
    rxn : Chem.Reaction
        Reaction object for polymerisation.
    NH2_side_chain : boolean
    COOH_side_chain : boolean

    """
    P1_AApos, P1_chain, P1_rxntype, P1_grp = P1_bond
    P2_AApos, _, P2_rxntype, P2_grp = P2_bond
    if P1_grp == "Acc":
        Gly_Acc = P1
        Gly_Dnr = P2
    elif P1_grp == "Dnr":
        Gly_Acc = P2
        Gly_Dnr = P1
    else:
        raise ValueError

    rxn = GLYCAN_RXN_TYPES["Glycosidic"]
    pdts_smiles = run_mono_pdt_rxn(rxn, Gly_Acc.mol, Gly_Dnr.mol)
    name_ext = f"{polymerisation_name}"
    return pdts_smiles, name_ext, rxn, False, False


# %% Molecules

def state_molecular_comb_logic():
    print(
        '''Molecular Combination Logic
    __add__ (+): \nGlycosidic bond (1,4) or \nPeptide bond with no side chain peptidation (self = COOH, other = NH2)
    __and__ (&): \nGlycopeptide bond (requires overlapping Lac moiety) or \nPeptide bond with side chain peptidation on self COOH. Does not extend code, retains code of main chain (self).
    __sub__ (-): \nPeptide bond with side chain peptidation on other NH2. Does not extend code, retains code of main chain (other).\n''')


class Molecule:
    '''Represents a molecule. Wraps around a Chem.Mol object.'''
    abbrevs = rdAbbreviations.GetDefaultAbbreviations()

    def __init__(self, smiles):
        """

        Parameters
        ----------
        smiles : string
            Molecule SMILES.

        Returns
        -------
        None.

        """
        self.mol = Chem.MolFromSmiles(smiles)
        self.add_coords = False
        self.smiles = Chem.MolToSmiles(self.mol, canonical=True)
        self.retention_time = 0
        # Hidden Properties
        self.reset_hidden_var()
        # mz
        self.mz = {}
        if AUTOCALC_MZ:
            self.update_mz()

    def reset_hidden_var(self):
        # Hidden Properties
        self.__property_array = None
        self.__morgan_fingerprint = {}
        self.__name = None
        self.__synonym = None
        self.__formula = None
        self.__mMass = None
        self.__charge = None
        self.__neutral_form = None
        self.__InchiKey = None
        self.__neutral_InchiKey = None
        self.__clogP = None

    # %%%% Magic Methods

    def __repr__(self):
        return f"{self.smiles}"

    def __eq__(self, other):
        '''Compares SMILES'''
        return isinstance(other, Molecule) and self.smiles == other.smiles

    def __add__(self, other):
        '''no side chain peptides'''
        return self.peptide_bond(other)

    def __str__(self):
        '''shows longer name'''
        return self.smiles

    # %%%% Properties

    @property
    def name(self):
        '''Returns shorter name'''
        if self.__name is None:
            self.__name = shorten_molecule_name(self.synonym)
        return self.__name

    @property
    def synonym(self):
        '''Returns longer name'''
        if self.__synonym is None:
            self.__synonym = str(self)
        return self.__synonym

    @property
    def formula(self):
        '''Returns chemical formula'''
        if self.__formula is None:
            self.__formula = CalcMolFormula(self.mol)
        return self.__formula

    @property
    def mMass(self):
        '''Returns exact molecular mass'''
        if self.__mMass is None:
            self.__mMass = CalcExactMolWt(self.mol)
        return round(self.__mMass, PRECISION_MASS)

    @property
    def charge(self):
        '''Returns num of positive charge'''
        if self.__charge is None:
            self.__charge = GetFormalCharge(self.mol)
        return self.__charge

    @property
    def neutral_form(self):
        '''Returns neutral form of ion (if possible)'''
        if self.__neutral_form is None:
            if self.charge == 0:
                self.__neutral_form = self.mol
            else:
                self.__neutral_form = uncharge_molecule(self.mol)
        return self.__neutral_form

    @property
    def InchiKey(self):
        '''Returns InchiKey'''
        if self.__InchiKey is None:
            self.__InchiKey = MolToInchiKey(self.mol)
        return self.__InchiKey

    @property
    def neutral_InchiKey(self):
        '''Returns InchiKey of neutral form'''
        if self.__neutral_InchiKey is None:
            if self.neutral_form:
                self.__neutral_InchiKey = MolToInchiKey(self.neutral_form)
            else:
                self.__neutral_InchiKey = False
        return self.__neutral_InchiKey

    @property
    def clogP(self):
        '''Returns computed logP (Crippen method)'''
        if self.__clogP is None:
            self.__clogP = round(MolLogP(self.mol), 5)
        return self.__clogP

    @property
    def num_heavy_atoms(self):
        '''Returns number of non-H atoms'''
        return CalcNumHeavyAtoms(self.mol)

    @property
    def num_acidic(self):
        '''Returns number of acidic groups'''
        return len(self.mol.GetSubstructMatches(
            MISC_PATT_TYPES["Acidic"]))

    @property
    def num_basic(self):
        '''Returns number of basic groups'''
        return len(self.mol.GetSubstructMatches(
            MISC_PATT_TYPES["Basic"]))

    @property
    def num_hydroxyl(self):
        '''Returns number of hydroxyl groups'''
        return len(self.mol.GetSubstructMatches(
            MISC_PATT_TYPES["Hydroxyl"]))

    @property
    def num_amino(self):
        '''Returns number of amino groups'''
        return len(self.mol.GetSubstructMatches(
            MISC_PATT_TYPES["Amino"]))

    @property
    def num_peptide_bond(self):
        '''Returns number of peptide bonds'''
        return len(self.mol.GetSubstructMatches(
            MISC_PATT_TYPES["Peptide"]))

    @property
    def num_eupeptide_bond(self):
        '''Returns number of eupeptide bonds'''
        return len(self.mol.GetSubstructMatches(
            PEPTIDE_PATT_TYPES["Eupeptide"]))

    @property
    def num_dipeptide_bond(self):
        '''Returns number of dipeptide (eupeptide) bonds'''
        return len(self.mol.GetSubstructMatches(
            PEPTIDE_PATT_TYPES["Dipeptide"]))

    @property
    def num_isoglx_peptide_bond(self):
        '''Returns number of iso-Glx type bonds'''
        return len(self.mol.GetSubstructMatches(
            PEPTIDE_PATT_TYPES["isoGlx Dipeptide"]))

    @property
    def num_glycosidic_bond(self):
        '''Returns number of glycosidic bonds'''
        return len(self.mol.GetSubstructMatches(
            MISC_PATT_TYPES["Gly Bond"]))

    @property
    def num_diglycosidic_bond(self):
        '''Returns number of glycosidic bonds'''
        return len(self.mol.GetSubstructMatches(
            MISC_PATT_TYPES["DiGly Bond"]))

    @property
    def num_glycosidic_donor(self):
        '''Returns number of glycosidic donors'''
        return len(self.mol.GetSubstructMatches(
            MISC_PATT_TYPES["Gly Donor"]))

    @property
    def num_glycosidic_acceptor(self):
        '''Returns number of glycosidic acceptors'''
        return len(self.mol.GetSubstructMatches(
            MISC_PATT_TYPES["Gly Acceptor"]))

    @property
    def num_glycan(self):
        '''Returns number of glycan units'''
        return self.num_glycosidic_bond + self.num_glycosidic_acceptor

    @property
    def num_lactoyl_bond(self):
        '''Returns number of lactoyl bonds'''
        return len(self.mol.GetSubstructMatches(
            MISC_PATT_TYPES["Lactoyl Bond"]))

    @property
    def num_lactate(self):
        '''Returns number of free lactate groups'''
        return len(self.mol.GetSubstructMatches(
            MISC_PATT_TYPES["Lactate"]))

    @property
    def TPSA(self):
        '''Returns calculated total polar surface area'''
        return round(CalcTPSA(self.mol), 3)

    @property
    def property_array(self):
        '''Returns an array of properties for ML'''
        try:
            self.__property_array
        except AttributeError:
            self.__property_array = None
        if self.__property_array is None or len(self.__property_array) != 9:
            # save in hidden variable
            self.__property_array = np.array([
                self.num_heavy_atoms/1000,
                self.num_peptide_bond/100,
                self.num_dipeptide_bond/100,
                self.num_isoglx_peptide_bond/100,
                self.num_glycosidic_bond/100,
                self.num_diglycosidic_bond/100,
                self.num_lactoyl_bond/100,
                self.num_lactate/100,
                self.TPSA/5000
            ])
        return self.__property_array

    # %%%% Methods

    def morgan_fingerprint(self, radius, nbits):
        """
        Retrieve Morgan fingerprint of molecule

        Parameters
        ----------
        radius : integer
            Radius of Morgan fingerprint
        nbits : integer
            Bits to be used for Morgan fngerprint

        Returns
        -------
        bitvect : rdKit.BitVect
            BitVect object from rdKit

        """
        bitvect = self.__morgan_fingerprint.get((radius,nbits),None)
        if bitvect is None:
            bitvect = AllChem.GetMorganFingerprintAsBitVect(self.mol, radius, nbits)
        return bitvect

    def add_chem_coords(self):
        """
        Adds chemical co-ordinates to mol object.

        Returns
        -------
        None.

        """
        Chem.rdCoordGen.AddCoords(self.mol)
        self.add_coords = True

    def refresh(self):
        """
        Sanitizes self.mol (Chem.Mol object) and refreshes SMILES after
        a reaction has occurred in-place.

        Returns
        -------
        None.

        """
        Chem.SanitizeMol(self.mol)
        self.smiles = Chem.MolToSmiles(self.mol, canonical=True)
        self.add_coords = False
        self.reset_hidden_var() # reset

    def set_atom_map_numbers(self):
        """
        Numbers atoms.

        Returns
        -------
        None.

        """
        for atom in self.mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())

    def clear_atom_map_numbers(self):
        """
        Clear atom numbers.

        Returns
        -------
        None.

        """
        for atom in self.mol.GetAtoms():
            atom.SetAtomMapNum(0)

    def update_mz(self):
        """
        Updates m/z of hard-coded adducts.

        Returns
        -------
        None.

        """
        self.mz = {"[M+H]+": get_cation_mz(self.mMass, "H", 1),
                   "[M+2H]2+": get_cation_mz(self.mMass, "H", 2),
                   "[M+3H]3+": get_cation_mz(self.mMass, "H", 3),
                   "[M+4H]4+": get_cation_mz(self.mMass, "H", 4),
                   "[M+5H]5+": get_cation_mz(self.mMass, "H", 5),
                   "[M+Na]+": get_cation_mz(self.mMass, "Na", 1),
                   "[M+2Na]2+": get_cation_mz(self.mMass, "Na", 2),
                   "[M+K]+": get_cation_mz(self.mMass, "K", 1),
                   "[M+2K]2+": get_cation_mz(self.mMass, "K", 2)}

    def draw(self):
        """
        Display image of molecule.

        Returns
        -------
        PIL Image
            Image of molecule.

        """
        if not self.add_coords:
            self.add_chem_coords()
        mol = rdAbbreviations.CondenseMolAbbreviations(
            self.mol, Molecule.abbrevs)
        return Draw.MolToImage(mol, size=(800, 400))

    def uncharge(self):
        """
        Uncharges self and refreshes if successful.

        Returns
        -------
        uncharged : Chem.Mol object
            Returns uncharged Chem.Mol object if successful else False.

        """
        uncharged = uncharge_molecule(self.mol)
        if uncharged:
            self.mol = uncharged
            self.refresh()
        return uncharged

    def peptide_bond(self, other,
                     COOH_side_chain=False, NH2_side_chain=False):
        """
        Forms a peptide bond with (other).

        Parameters
        ----------
        other : AminoAcid, Peptide


        Raises
        ------
        NoReactionError
            No products.
        TooManyPdtsError
            More than one valid product.

        Returns
        -------
        Peptide
            Peptide molecule.

        """
        if not isinstance(other, (AminoAcid, Peptide)):
            raise NoReactionError(self, other, "invalid")

        rxn = create_peptide_rxn(self, other, COOH_side_chain, NH2_side_chain)

        pdt_peptide_code = other.peptide_code
        pdt_eupeptide_smarts = copy.deepcopy(other.eupeptide_smarts)
        pdt_eupeptide_smarts.pop("NH2")  # prevent further reaction at N
        pdt_side_chain_smarts = copy.deepcopy(other.side_chain_smarts)
        pdts_smiles = run_mono_pdt_rxn(rxn, self.mol, other.mol)
        pdts = [AminoAcid(smiles=smiles,
                          peptide_code=pdt_peptide_code,
                          eupeptide_smarts=pdt_eupeptide_smarts,
                          side_chain_smarts=pdt_side_chain_smarts) for smiles in set(pdts_smiles)]

        if len(pdts) == 1:
            return pdts[0]
        elif len(pdts) == 0:
            raise NoReactionError(self, other, "peptide", rxn)
        else:
            raise TooManyPdtsError(self, other, "peptide", rxn, pdts_smiles)


# %%% Amino Acids


class AminoAcid(Molecule):
    '''Represents an Amino Acid.'''

    def __init__(self, smiles="", peptide_code="",
                 eupeptide_smarts={"NH2": None, "COOH": None},
                 side_chain_smarts={}):
        """
        Initialization with no parameters provided will create a "Null" AA

        Parameters
        ----------
        smiles : string
            Molecule SMILES.
        peptide_code : string
            String representing Amino Acid.
        eupeptide_smarts : dict
            Dict with eupeptide reactions.
        side_chain_smarts : dict
            Dict with side chain reactions.

        Returns
        -------
        None.

        """
        super().__init__(smiles)
        if type(peptide_code) == str:
            # tracks main chain; N terminus --> C terminus
            self.peptide_code = [peptide_code]
        else:
            self.peptide_code = peptide_code

        # {idx: [side_chain_code] where code[idx] is connected to side chain}
        self.side_chain_dict = {1: (None, None)}
        # idx == 1 denotes the first AA
        self.eupeptide_smarts = eupeptide_smarts  # {"COOH","NH2" :smarts}
        self.side_chain_smarts = side_chain_smarts  # {"COOH","NH2": smarts}

    # %%%% Magic Methods

    def __add__(self, other):
        '''Peptide bond with no side chain reactivity.'''
        return self.peptide_bond(other)

    def __and__(self, other):
        '''Peptide bond with side chain reactivity for self (if present).'''
        return self.peptide_bond(other, COOH_side_chain=True)

    def __sub__(self, other):
        '''Peptide bond with side chain reactivity for other (if present).'''
        return self.peptide_bond(other, NH2_side_chain=True)

    def __len__(self):
        '''Number of residues (minus side chains)'''
        return len(self.peptide_code)

    def __str__(self):
        return self.peptide_code[0]

    # %%%% Properties

    @property
    def all_amino_acids(self):
        '''Returns all amino acids present in peptide'''
        all_sc_AAs = flatten(AAs for grp, AAs in
                             self.side_chain_dict.values() if AAs is not None)
        all_AAs = copy.deepcopy(self.peptide_code)
        all_AAs.extend(all_sc_AAs)
        return all_AAs

    @property
    def all_side_chains(self):
        '''Returns all side chains present in peptide'''
        side_chains = [side_chain for _, side_chain in
                       self.side_chain_dict.values() if side_chain is not None]
        return side_chains

    @property
    def degree_amidation(self):
        '''Returns number of amidated (CONH2) amino acids in peptide'''
        return sum(AA in AMINO_ACID_AMIDATED_LST for AA in self.all_amino_acids)

    @property
    def len_main(self):
        '''Number of residues (minus side chains)'''
        return len(self)

    @property
    def len_side(self):
        '''Number of residues (only side chains)'''
        len_side = sum(len(side_chain) for side_chain
                       in self.all_side_chains)
        return len_side

    @property
    def len_null(self):
        '''Number of null residues in main chain'''
        return self.peptide_code.count("")

    @property
    def peptide_units(self):
        '''Number of non-null residues in main chain'''
        return self.len_main - self.len_null

    @property
    def ontology(self):
        '''Peptide ontology'''
        if self.len_side == 0:
            # should null residues be counted?
            return f"P{self.len_main-self.len_null}"
        else:
            return BRANCH_CHAR.join([f"P{self.len_main}", f"S{self.len_side}"])

    # %%%% Methods

    def is_smarts_valid(self, grp, side_chain=False):
        """
        Checks if amino acid has SMARTS code for 'NH2'/'COOH'.

        Parameters
        ----------
        grp : string
            'NH2' or 'COOH'; returns False by default otherwise.
        side_chain : boolean, optional
            If True, checks the SMARTS for the side chain, else checks SMARTS
            for the eupeptide bond. The default is False.

        Returns
        -------
        boolean

        """
        if grp not in ["COOH", "NH2"]:
            return False
        if side_chain:
            smarts = self.side_chain_smarts.get(grp,None)
        else:
            smarts = self.eupeptide_smarts.get(grp,None)
        return type(smarts) == str

    def has_side_chain(self, pos):
        """
        Checks if peptide has side chain.

        Parameters
        ----------
        pos : integer
            Integer representing position in peptide, 1 = N-terminus.

        Returns
        -------
        boolean

        """
        if pos not in self.side_chain_dict:
            return False
        else:
            return self.side_chain_dict[pos] != (None, None)

    def generate_updated_side_chain(self, other, side_chain_grp=None,
                                    inplace=False, idx=None):
        """
        Creates copy of side_chain_dict, updates and returns it.

        Parameters
        ----------
        other : AminoAcid, Peptide, None
            Object to be added as either main or side chain.
        side_chain_grp : string, optional
            Describes connection between side chain and stem peptide.
            Can be either "NH2" or "COOH". The default is None. If None,
            connection will be done on main chain.
        inplace : boolean, optional
            If True, modifies side_chain_dict inplace. The default is False.
        idx : integer, optional
            Integer representing position in peptide, 1 = N-terminus.
            Defaults to C-terminus if None.

        Raises
        ------
        ValueError
            Incorrect bonding provided

        Returns
        -------
        scd : dict
            Updated copy of side_chain_dict

        """
        if idx is None:
            idx = len(self)
        if inplace:
            scd = self.side_chain_dict
        else:
            scd = copy.deepcopy(self.side_chain_dict)

        # not possible to have 2 side chains?
        if isinstance(other, (AminoAcid, Peptide)):
            peptide_code = other.peptide_code
            if side_chain_grp == "NH2":
                # tracks side chain; C terminus --> N terminus
                scd[idx] = side_chain_grp, list(reversed(peptide_code))
            elif side_chain_grp == "COOH":
                # tracks side chain; N terminus --> C terminus
                scd[idx] = side_chain_grp, peptide_code
            elif side_chain_grp is None:
                # addition as main chain, check if self inherits other's scd
                other_scd = copy.deepcopy(other.side_chain_dict)
                for other_idx,sc in other_scd.items():
                    scd[other_idx+idx]=sc
            else:
                raise ValueError(f"Unknown side_chain_grp: {side_chain_grp}")

        elif idx not in scd:
            # indicate Null, don't overwrite
            scd[idx] = (None, None)

        return scd

    def peptide_bond(self, other, COOH_side_chain=False, NH2_side_chain=False):
        """
        Forms a peptide bond with (other). Joins COOH of self to NH2 of other.
        X_side_chain >> use alternative COOH/NH2.

        Parameters
        ----------
        other : AminoAcid, Peptide
        COOH_side_chain : boolean, optional
            If True, use alternative COOH. The default is False.
        NH2_side_chain : boolean, optional
            If True, use alternative NH2. The default is False.

        Raises
        ------
        NoReactionError
            No products.
        TooManyPdtsError
            More than one valid product.

        Returns
        -------
        Peptide
            Peptide molecule.

        """

        if self.peptide_code[-1] == "" and not COOH_side_chain:
            return self.null_peptide_bond(other, null="C")
        elif other.peptide_code[0] == "" and not NH2_side_chain:
            return self.null_peptide_bond(other, null="N")

        # check for validity of side chains
        if COOH_side_chain:
            assert self.is_smarts_valid("COOH", side_chain=True),\
                f"C-terminus of {self.name} has no side chain COOH"
        if NH2_side_chain:
            assert other.is_smarts_valid("NH2", side_chain=True),\
                f"N-terminus of {other.name} has no side chain NH2"

        rxn = create_peptide_rxn(
            self, other, COOH_side_chain, NH2_side_chain)

        # Tracking of side chains and main chains...
        if COOH_side_chain:  # addition of other as side chain
            pdt_peptide_code = self.peptide_code  # not main chain
            side_chain_dict = self.generate_updated_side_chain(other, "COOH")

        elif NH2_side_chain:  # addition of self as side chain
            pdt_peptide_code = other.peptide_code
            side_chain_dict = other.generate_updated_side_chain(self, "NH2")

        else:  # normal behaviour, addition to main chain
            pdt_peptide_code = self.peptide_code+other.peptide_code
            # add placeholder
            side_chain_dict = self.generate_updated_side_chain(other, None)

        # Tracking of future peptide bond reactions
        if COOH_side_chain:  # addition of other as side chain
            # retain current eupeptide reactivity
            pdt_eupeptide_smarts = copy.deepcopy(self.eupeptide_smarts)
            # retain self side chain reactivity minus COOH
            pdt_side_chain_smarts = copy.deepcopy(self.side_chain_smarts)
            pdt_side_chain_smarts.pop("COOH")  # prevent further rxn at side C

        elif NH2_side_chain:  # addition of self as side chain
            # inherit other's eupeptide reactivity
            pdt_eupeptide_smarts = other.eupeptide_smarts
            # inherit other's side chain reactivity minus NH2
            pdt_side_chain_smarts = copy.deepcopy(other.side_chain_smarts)
            pdt_side_chain_smarts.pop("NH2")  # prevent further rxn at side C

        else:  # normal behaviour, addition to main chain
            # retain self eupeptide NH2 reactivity
            pdt_eupeptide_smarts = copy.deepcopy(self.eupeptide_smarts)
            # inherit other eupeptide COOH and side_chain reactivity
            pdt_eupeptide_smarts["COOH"] = other.eupeptide_smarts["COOH"]
            pdt_side_chain_smarts = copy.deepcopy(other.side_chain_smarts)

        pdts_smiles = run_mono_pdt_rxn(rxn, self.mol, other.mol)
        pdts = [Peptide(smiles=smiles,
                        peptide_code=pdt_peptide_code,
                        eupeptide_smarts=pdt_eupeptide_smarts,
                        side_chain_smarts=pdt_side_chain_smarts,
                        side_chain_dict=side_chain_dict) for smiles in set(pdts_smiles)]

        for x in pdts:
            # add placeholder
            x.generate_updated_side_chain(None, inplace=True)

        if len(pdts) == 1:
            return pdts[0]
        elif len(pdts) == 0:
            raise NoReactionError(self, other, "peptide", rxn)

        else:
            pick_min = NH2_side_chain and COOH_side_chain
            pdt = pick_correct_pdt(pdts, pick_min)
            if pdt:
                return pdt
            else:
                raise TooManyPdtsError(
                    self, other, "peptide", rxn, pdts_smiles)

    def null_peptide_bond(self, other, null="C"):
        """
        Forms a peptide bond with (other) when other is a null amino acid.

        Parameters
        ----------
        other : AminoAcid, Peptide
        null : string, optional
            Indicates terminus in which null is added. The default is "C".

        Raises
        ------
        ValueError
            Incorrect bonding provided

        Returns
        -------
        Peptide
            Peptide molecule.

        """
        if null == "C":
            smiles = other.smiles
            pdt_peptide_code = self.peptide_code+other.peptide_code
            pdt_eupeptide_smarts = {}
            # inherit other eupeptide COOH and side chain
            pdt_eupeptide_smarts["COOH"] = other.eupeptide_smarts["COOH"]
            pdt_eupeptide_smarts["NH2"] = self.eupeptide_smarts["NH2"]
            pdt_side_chain_smarts = copy.deepcopy(other.side_chain_smarts)
            side_chain_dict = self.generate_updated_side_chain(other,None)
            pdt = Peptide(smiles=smiles,
                          peptide_code=pdt_peptide_code,
                          eupeptide_smarts=pdt_eupeptide_smarts,
                          side_chain_smarts=pdt_side_chain_smarts,
                          side_chain_dict=side_chain_dict)
            return pdt
        elif null == "N":
            smiles = self.smiles
            pdt_peptide_code = self.peptide_code+other.peptide_code
            pdt_eupeptide_smarts = {}
            # inherit self eupeptide NH2 and side chain
            pdt_eupeptide_smarts["COOH"] = other.eupeptide_smarts["COOH"]
            pdt_eupeptide_smarts["NH2"] = self.eupeptide_smarts["NH2"]
            pdt_side_chain_smarts = copy.deepcopy(self.side_chain_smarts)
            side_chain_dict = other.generate_updated_side_chain(self,None)
            pdt = Peptide(smiles=smiles,
                          peptide_code=pdt_peptide_code,
                          eupeptide_smarts=pdt_eupeptide_smarts,
                          side_chain_smarts=pdt_side_chain_smarts,
                          side_chain_dict=side_chain_dict)
            return pdt
        else:
            raise ValueError

    def check_AA_identity(self, pos, identity):
        """
        Checks if amino acid at position (pos) has a specific identity.'

        Parameters
        ----------
        pos : integer
            Integer representing position in peptide, 1 = N-terminus.
        identity : list
            List of identities to check against. Can be names of amino acids
            or groups or None. Valid entries for groups include: "diamino",
            "dicarboxy" and "all". "all" does not include null amino acids. None
            indicates shorter peptides are accepted; i.e. length < pos

        Returns
        -------
        boolean

        """
        valid_AAs = list(copy.deepcopy(identity))
        # Account for special names
        if "diamino" in identity:
            valid_AAs.extend(AMINO_ACID_DIAMINO)
            valid_AAs.remove("diamino")
        if "dicarboxy" in identity:
            valid_AAs.extend(AMINO_ACID_DICARBOXY)
            valid_AAs.remove("dicarboxy")
        if "any" in identity:
            valid_AAs.extend(ALL_AMINO_ACID_LST)
            valid_AAs.remove("any")
        if len(self) < pos:
            # check if None is valid
            return None in valid_AAs
        else:
            return self.peptide_code[pos-1] in valid_AAs

    def check_side_chain_identity(self, pos, side_chains):
        """
        Checks if peptide has specific side chain.

        Parameters
        ----------
        pos : integer
            Integer representing position in peptide, 1 = N-terminus.
        side_chains : list
            List of strings that represent valid side chains.

        Returns
        -------
        boolean

        """
        if pos not in self.side_chain_dict:
            return False
        # Account for special names
        elif "any" in side_chains:
            return True
        else:
            g, c = self.side_chain_dict[pos]
            return c in side_chains

    def peptide_comp_diffcalc(self, pep_range):
        """
        Calculates differences in peptide composition.

        Parameters
        ----------
        pep_range : list
            List of strings which describes valid amino acids.

        Returns
        -------
        integer

        """
        if pep_range is None:
            return 0
        else:
            diff_AAs = [AA for AA in self.peptide_code
                        if AA not in pep_range and AA != ""]
            return len(diff_AAs)

    def side_chain_comp_diffcalc(self, sc_range):
        """
        Calculates differences in side chain composition.

        Parameters
        ----------
        sc_range : list
            List of lists which describes valid side chains.

        Returns
        -------
        integer

        """
        if sc_range is None:
            return 0
        else:
            diff_side_chains = [
                sc for sc in self.all_side_chains if sc not in sc_range]
            return len(diff_side_chains)

# %%% Peptides


class Peptide(AminoAcid):

    def __init__(self, smiles, peptide_code, eupeptide_smarts,
                 side_chain_smarts, side_chain_dict):
        '''

        Parameters
        ----------
        smiles : string
            Molecule SMILES.
        peptide_code : list
            Name of peptide.
        eupeptide_smarts : dict
            Dict with eupeptide reactions.
        side_chain_smarts : dict
            Dict with side chain peptide reactions.
        side_chain_dict : dict
            Dict describing side chains.

        Returns
        -------
        None.

        '''
        super().__init__(smiles, peptide_code, eupeptide_smarts, side_chain_smarts)
        # supersede null dictionary created by AA
        self.side_chain_dict = side_chain_dict

    # %%%% Magic Methods

    def __str__(self):
        main_string = BOND_CHAR.join(r for r in self.peptide_code if r != "")
        side_chain_lst = []

        for pos in range(1, len(self)+1):
            grp, side_chain = self.side_chain_dict.get(pos, (None, None))
            if grp is None or side_chain is None:
                continue
            else:
                side_chain_lst.append(
                    f"{pos}{self.get_side_chain_string(grp,side_chain)}")

        if len(side_chain_lst) == 0:
            return main_string
        else:
            return f"{main_string}[{','.join(side_chain_lst)}]"

    # %%%% Methods

    def state_peptide_str(self):
        """
        Prints the structure of the peptide.

        Returns
        -------
        None.

        """
        for i, peptide_code in enumerate(self.peptide_code):
            pos = i+1
            grp, side_chain = self.side_chain_dict.get(pos, (None, None))

            if grp:
                print(
                    f"AA{pos}: \t{peptide_code}{self.get_side_chain_string(grp,side_chain)}")
            else:
                print(f"AA{pos}: \t{peptide_code}")

    def get_side_chain_string(self, grp, side_chain):
        """
        Formats the side chain at grp as a string.

        Parameters
        ----------
        grp : string
            Describes connection between stem peptide and side chain.
            Can be 'COOH' or 'NH2'.
        side_chain : list
            List of strings that represent the side chain.

        Returns
        -------
        str
            Formatted string that describes side chain.

        """
        # !!! Known issue: side chains present in side_chain are omitted
        return f"{BRANCH_CHAR}{grp}{BRANCH_CHAR}{BOND_CHAR.join(side_chain)}"


# %%% Glycans


class Glycan(Molecule):
    code_exceptions = ["Lac"]

    def __init__(self, smiles, gly_code, init_gly_type, fin_gly_type):
        """


        Parameters
        ----------
        smiles : string
            Molecule SMILES.
        gly_code : str or list
            string(s) describing the glycan
        init_gly_type : string
            Designates the first glycan in the (1-->4) direction
            as "Glc" or "Mur".
        fin_gly_type : string
            Designates the last glycan in the (1-->4) direction
            as "Glc" or "Mur".

        Returns
        -------
        None.

        """
        super().__init__(smiles)
        if type(gly_code) == str:  # Tracks main chain of glycan (1-->4) direction
            self.gly_code = [gly_code]
        else:
            self.gly_code = gly_code
        # type of first glycan (1-->4) direction
        self.init_gly_type = init_gly_type
        # type of last glycan (1-->4) direction
        self.fin_gly_type = fin_gly_type

    # %%%% Magic Methods

    def __add__(self, other):
        if isinstance(other, Glycan):
            # Glycosidic bond
            return self.glycosidic_bond(other)
        else:
            return super().__add__(other)

    def __and__(self, other):
        if isinstance(other, (AminoAcid, Peptide)):
            return self.glycopeptide_bond(other)
        else:
            raise NoReactionError(self, other, "Invalid reaction.", None)

    def __len__(self):
        return len([gly for gly in self.gly_code
                    if gly not in Glycan.code_exceptions])

    def __str__(self):
        return BOND_CHAR.join(self.gly_code)

    # %%%% Properties

    @property
    def degree_acetylation(self):
        return sum(gly.count("Ac") for gly in self.gly_code)

    @property
    def ontology(self):
        return f"G{len(self)}"

    # %%%% Methods

    def glycosidic_bond(self, other):
        """
        Forms a 1,4 glycosidic bond between self and other.

        Parameters
        ----------
        other : Glycan

        Raises
        ------
        NoReactionError
            No products.
        TooManyPdtsError
            More than one valid product.

        Returns
        -------
        Glycan
            Returns the resulting glycan.

        """
        rxn = GLYCAN_RXN_TYPES["Glycosidic"]
        pdt_gly_code = self.gly_code+other.gly_code
        pdts_smiles = run_mono_pdt_rxn(rxn, self.mol, other.mol)
        pdts = [Glycan(smiles=smiles,
                       gly_code=pdt_gly_code,
                       init_gly_type=self.init_gly_type,
                       fin_gly_type=other.fin_gly_type) for smiles in set(pdts_smiles)]
        if len(pdts) > 1:
            raise TooManyPdtsError(self, other, "glycosidic", rxn, pdts_smiles)
        elif len(pdts) == 0:
            raise NoReactionError(self, other, "glycosidic", rxn)
        else:
            return pdts[0]

    def glycopeptide_bond(self, other):
        """
        Generates a peptidoglycan by combining glycan with peptide by
        overlapping the Lac moiety.

        Parameters
        ----------
        other : Peptide

        Raises
        ------
        NoReactionError
            No products.
        TooManyPdtsError
            More than one valid product.

        Returns
        -------
        Peptidoglycan

        """
        rxn = GLYCAN_RXN_TYPES["Lac Merger"]
        pdts_smiles = run_mono_pdt_rxn(rxn, self.mol, other.mol)
        pdts = [Peptidoglycan(smiles=smiles, glycan=self, peptide=other)
                for smiles in set(pdts_smiles)]
        if len(pdts) > 1:
            raise TooManyPdtsError(self, other, "Lac Merger", rxn, pdts_smiles)
        elif len(pdts) == 0:
            raise NoReactionError(self, other, "Lac Merger", rxn)
        else:
            return pdts[0]

    def check_glycan_identity(self, identity):
        #!!! TODO
        '''
        To implement and support special groups of glycans:
            Glc-type
            Mur-type
            anhydro
            reduced
        '''
        pass

    def glycan_comp_diffcalc(self, glycan_range):
        """
        Calculates differences in glycan composition

        Parameters
        ----------
        glycan_range : List
            List of strings which describes valid glycans.


        Returns
        -------
        integer

        """
        if glycan_range is None or len(self.gly_code) == 0:
            return 0
        valid_glycans = list(copy.deepcopy(glycan_range))
        # Account for special names
        #!!! KIV move this to check_glycan_identity
        if "any" in glycan_range:
            valid_glycans.extend(ALL_GLC_LST)
            valid_glycans.extend(ALL_MUR_LST)
            valid_glycans.remove("any")
            valid_glycans.append(None)
        if "Glc-type" in glycan_range:  # All Glucosamines
            valid_glycans.extend(ALL_GLC_LST)
            valid_glycans.remove("Glc-type")
        if "Mur-type" in glycan_range:  # All Muramic Acids
            valid_glycans.extend(ALL_MUR_LST)
            valid_glycans.remove("Mur-type")
        if "anhydroMur-type" in glycan_range:  # All anhydro-muramic acids
            valid_glycans.extend(ALL_ANHYDROMUR_LST)
            valid_glycans.remove("anhydroMur-type")
        if "Mur[r]-type" in glycan_range:  # All reduced muramic acids
            valid_glycans.extend(ALL_REDMUR_LST)
            valid_glycans.remove("Mur[r]-type")
        else:
            return len(self.gly_code) - sum(g in valid_glycans for g in self.gly_code)

# %%% Peptidoglycans


class Peptidoglycan(Molecule):

    modifications = {
        "Alanine/Lactate Substitution":
            "Last Ala in pentapeptides substituted with lactate",
        "Muramic Lactam":
            "Lactam formed between acetyl and lactoyl groups",
        "Lactoyl Peptide":
            "Glycans removed with lactoyl moiety remaining",
        "Amidase":
            "Amidase acting on lactoyl peptide bond",
        "EPase P1":
            "Endopeptidase acting on 1st stem peptide bond",
        "EPase P2":
            "Endopeptidase acting on 2nd stem peptide bond",
        "Braun LPP":
            "Addition of ε-Lys-Arg dipeptide for residues 4-5"
    }

    def __init__(self, smiles, glycan=None, peptide=None,
                 modifications=None, polymerisation_type=None, linked_PGN=None):
        """

        Parameters
        ----------
        smiles : string
            Molecule SMILES.
        glycan : Glycan, optional
            The default is None.
        peptide : Peptide, optional
            The default is None.
        modifications : dict or list, optional
            Dictionary of modifications present. Use dict to inherit
            modifications. Use list for modifications that are present. The
            default is None.
        polymerisation_type : string, optional
            The type of polymerisation. The default is None.
        linked_PGN : Peptidoglycan, optional
            The default is None.

        Returns
        -------
        None.

        """
        super().__init__(smiles)
        self.glycan = copy.deepcopy(glycan)  # Glycan object
        self.peptide = copy.deepcopy(peptide)  # Peptide object
        if modifications is None:
            self.modifications = {
                mod: False for mod in Peptidoglycan.modifications}
        elif type(modifications) == dict:  # inherit mods
            self.modifications = copy.deepcopy(modifications)
        elif type(modifications) == list:  # set all to False except listed
            self.modifications = {}
            for mod in Peptidoglycan.modifications:
                if mod in modifications:
                    self.modifications[mod] = True
                else:
                    self.modifications[mod] = False
        else:
            raise ValueError

        self.update_mz()
        self.polymerisation_type = polymerisation_type
        self.linked_PGN = copy.deepcopy(linked_PGN)  # PGN object
        if self.linked_PGN is None:
            self.index = 1
        else:
            self.index = self.linked_PGN.index+1

    # %%%% Magic Methods

    def __str__(self):  # works for monomers only
        if self.polymerisation_type is None:
            return f"{str(self.glycan)}{BRANCH_CHAR}{str(self.peptide)}"
        else:
            return f'''{str(self.glycan)}{BRANCH_CHAR}{str(self.peptide)};({self.polymerisation_type}){BRANCH_CHAR}{str(self.linked_PGN)}'''

    # %%%% Properties

    @property
    def ontology(self):
        components = [self.glycan, self.peptide, self.linked_PGN]
        ont_components = [c.ontology for c in components if c is not None]
        ont = BRANCH_CHAR.join(ont_components)
        return ont

    @property
    def degree_acetylation(self):
        # get self dAc
        if self.glycan is None:
            dAc = 0
        else:
            dAc = self.glycan.degree_acetylation
        # check if polymer
        if self.index == 1:
            return dAc
        else:
            return dAc + self.linked_PGN.degree_acetylation

    @property
    def degree_amidation(self):
        # get self dAc
        if self.peptide is None:
            dAmi = 0
        else:
            dAmi = self.peptide.degree_amidation
        # check if polymer
        if self.index == 1:
            return dAmi
        else:
            return dAmi + self.linked_PGN.degree_amidation

    @property
    def gly_code(self):
        if self.glycan is None:
            return [None]
        else:
            return self.glycan.gly_code

    @property
    def peptide_code(self):
        if self.peptide is None:
            return [None]
        else:
            return self.peptide.peptide_code

    # %%%% Methods

    def convert_to_dict(self):
        """
        Returns dictionary representation of PGN.

        Returns
        -------
        d : dictionary
            Dictionary representation of PGN.

        """
        if self.linked_PGN is None:
            linked_PGN_name = "None"
        else:
            linked_PGN_name = self.linked_PGN.name
        d = {
            "Name": self.name,
            "Synonym": self.synonym,
            "Mol": self.mol,
            "Formula": self.formula,
            "Monoisotopic Mass": self.mMass,
            "Modifications": self.mods,
            "Degree Amidation": self.degree_amidation,
            "Degree Acetylation": self.degree_acetylation,
            "Ontology": self.ontology,
            "PGN Units": self.index,
            "Glycan Units": self.glycan_len(),
            "Peptide Length": self.peptide_len(),
            "Peptide Units": self.peptide_units(),
            "Polymerised Unit": linked_PGN_name,
            "Polymerisation Type": self.polymerisation_type,
            "Glycan": self.glycan_name(),
            "Peptide": self.peptide_name(),
            "SMILES": self.smiles,
            "INCHIKEY": self.InchiKey,
            "clogP": self.clogP,
            "RT": self.retention_time
        }
        d.update({adduct: self.mz[adduct] for adduct in OUTPUT_ADDUCTS})
        return d

    def refresh(self):
        """
        Refreshes molecule after reaction in-place.

        Returns
        -------
        None.

        """
        super().refresh()
        for x in self.peptide, self.glycan, self.linked_PGN:
            if x:
                x.refresh()
        self.update_mz()

    def state_peptide_str(self, index=None):
        """
        Prints the structure of the peptide.
        index = None : prints self
        index = self.index: prints self
        index = int: prints for linked PGN with that index

        Parameters
        ----------
        index : integer, optional
            See above. The default is None.

        Returns
        -------
        None.

        """
        if index not in (None, self.index):
            self.linked_PGN.state_peptide_str(index=index)
        elif self.peptide is None:
            pass
        else:
            self.peptide.state_peptide_str()

    def peptide_units(self, index=None):
        """
        Returns the no. of non-null residues in the peptide(s).
        index = None : returns total
        index = self.index: returns self
        index = int: returns for PGN with that index

        Parameters
        ----------
        index : integer, optional
            See above. The default is None.

        Returns
        -------
        integer
            Number of non-null residues in peptide.

        """
        if index not in (None, self.index) and index < self.index:
            return self.linked_PGN.peptide_len(index=index)
        elif index is None:
            return sum(self.peptide_len(index=i)
                       for i in range(1, self.index+1))
        elif self.peptide is None:
            return 0
        else:
            return self.peptide.peptide_units

    def peptide_len(self, index=None):
        """
        Returns the length of the peptide(s) - includes null residues but excludes
        bridge peptides.
        index = None : returns total
        index = self.index: returns self
        index = int: returns for PGN with that index

        Parameters
        ----------
        index : integer, optional
            See above. The default is None.

        Returns
        -------
        integer
            Length of peptide(s)

        """
        if index not in (None, self.index) and index < self.index:
            return self.linked_PGN.peptide_len(index=index)
        elif index is None:
            return sum(self.peptide_len(index=i)
                       for i in range(1, self.index+1))
        elif self.peptide is None:
            return 0
        else:
            return len(self.peptide)

    def glycan_len(self, index=None):
        """
        Returns the length of the glycan(s).
        index = None : returns total (self + all linked PGN)
        index = self.index: returns self
        index = int: returns for PGN with that index

        Parameters
        ----------
        index : integer, optional
            See above. The default is None.

        Returns
        -------
        integer
            Length of glycan(s)

        """
        if index not in (None, self.index) and index < self.index:
            return self.linked_PGN.glycan_len(index=index)
        elif index is None:
            return sum(self.glycan_len(index=i)
                       for i in range(1, self.index+1))
        elif self.glycan is None:
            return 0
        else:
            return len(self.glycan)

    def peptide_name(self, index=None):
        """
        Returns the name of stem peptide.
        index = None : returns self
        index = self.index: returns self
        index = int: returns for PGN with that index

        Parameters
        ----------
        index : integer, optional
            See above. The default is None.

        Returns
        -------
        string
            Name of stem peptide.

        """
        if index not in (None, self.index) and index < self.index:
            return self.linked_PGN.peptide_name(index=index)
        elif self.peptide is None:
            return None
        else:
            return self.peptide.name

    def glycan_name(self, index=None):
        """
        Returns the name of glycan chain.
        index = None : returns self
        index = self.index: returns self
        index = int: returns for PGN with that index

        Parameters
        ----------
        index : integer, optional
            See above. The default is None.

        Returns
        -------
        string
            Name of glycan chain.

        """
        if index not in (None, self.index) and index < self.index:
            return self.linked_PGN.glycan_name(index=index)
        elif self.glycan is None:
            return None
        else:
            return self.glycan.name

    def has_bridge(self, pos, index=None):
        """
        Checks if peptide has bridge peptide (side chain) at "pos".
        index = None : returns self
        index = self.index: returns self
        index = int: returns for PGN with that index

        Parameters
        ----------
        pos : integer
            Position of peptide to check.
        index : integer, optional
            See above. The default is None.

        Returns
        -------
        boolean

        """
        if index not in (None, self.index) and index < self.index:
            return self.linked_PGN.has_bridge(pos, index=index)
        elif self.peptide is None:
            return False
        else:
            #peptide retains side chain terminology (instead of bridge peptide)
            return self.peptide.has_side_chain(pos)

    def check_AA_identity(self, pos, identity, index=None):
        """
        Checks identity of amino acid in peptide at "pos".
        index = None : returns self
        index = self.index: returns self
        index = int: returns for PGN with that index

        Parameters
        ----------
        pos : integer
            Position of peptide to check.
        identity : list
            List of identities to check against. Can be names of amino acids
            or groups or None. Valid entries for groups include: "diamino",
            "dicarboxy" and "all". "all" does not include null amino acids. None
            indicates shorter peptides are accepted; i.e. length < pos.
        index : integer, optional
            See above. The default is None.

        Returns
        -------
        boolean

        """
        if index not in (None, self.index) and index < self.index:
            return self.linked_PGN.check_AA_identity(pos, identity, index=None)
        if self.peptide is None:
            return None in identity
        else:
            return self.peptide.check_AA_identity(pos, identity)

    def check_bridge_identity(self, pos, identity, index=None):
        """
        Checks if peptide has specific bridge peptide (side chain).
        index = None : returns self
        index = self.index: returns self
        index = int: returns for PGN with that index

        Parameters
        ----------
        pos : integer
            Position of peptide to check.
        identity : list
            List of valid bridge peptides (side chains) to check against.
        index : integer, optional
            See above. The default is None.

        Returns
        -------
        boolean

        """
        if index not in (None, self.index) and index < self.index:
            return self.linked_PGN.check_bridge_identity(pos, identity, index=None)
        if self.peptide is None:
            return None in identity
        else:
            #peptide retains side chain terminology (instead of bridge peptide)
            return self.peptide.check_side_chain_identity(pos, identity)

    def peptide_comp_diffcalc(self, pep_range, index=None):
        """
        Calculates differences in peptide composition.
        index = None : returns self
        index = self.index: returns self
        index = int: returns for PGN with that index

        Parameters
        ----------
        pep_range : list
            List of strings which describes valid amino acids.
        index : integer, optional
            See above. The default is None.

        Returns
        -------
        integer

        """
        if index not in (None, self.index) and index < self.index:
            return self.linked_PGN.peptide_comp_diffcalc(pep_range, index=index)
        elif self.peptide is None:
            return 0
        else:
            return self.peptide.peptide_comp_diffcalc(pep_range)

    def bridge_comp_diffcalc(self, bridge_range, index=None):
        """
        Calculates differences in bridge peptide composition.
        index = None : returns self
        index = self.index: returns self
        index = int: returns for PGN with that index

        Parameters
        ----------
        bridge_range : list
            List of lists which describes valid bridge peptides.
        index : integer, optional
            See above. The default is None.

        Returns
        -------
        integer

        """
        if index not in (None, self.index) and index < self.index:
            return self.linked_PGN.bridge_comp_diffcalc(bridge_range, index=index)
        elif self.peptide is None:
            return 0
        else:
            #peptide retains side chain terminology (instead of bridge peptide)
            return self.peptide.side_chain_comp_diffcalc(bridge_range)

    def glycan_comp_diffcalc(self, gly_range, index=None):
        """
        Calculates differences in glycan composition
        index = None : returns self
        index = self.index: returns self
        index = int: returns for PGN with that index

        Parameters
        ----------
        glycan_range : List
            List of strings which describes valid glycans.
        index : integer, optional
            See above. The default is None.

        Returns
        -------
        integer

        """
        if index not in (None, self.index) and index < self.index:
            return self.linked_PGN.glycan_comp_diffcalc(gly_range,
                                                        index=index)
        elif self.glycan is None:
            return 0
        else:
            return self.glycan.glycan_comp_diffcalc(gly_range)

    # %%% Modifications

    def get_modifications(self):
        """
        Returns list of valid (active) modifications in PGN.

        Returns
        -------
        list
            Valid (active) modifications in PGN.

        """
        return [x for x in self.modifications if self.modifications[x]]

    @property
    def mods(self):
        m = self.get_modifications()
        if len(m) > 0:
            return ", ".join(m)
        else:
            return ""

    def create_amidase_product(self):
        """
        Creates a PGN amidase product in-place (stem peptide absent lactate).

        Returns
        -------
        boolean
            Returns True if modification successful.

        """
        rxn = PEPTIDE_RXN_TYPES["Lac Removal"]
        if rxn.RunReactantInPlace(self.mol):
            self.refresh()
            self.modifications["Amidase"] = True
            return True
        else:
            return False

    def substitute_terminal_Ala_Lac(self):
        """
        Creates a modified PGN in-place with a terminal lactate (instead of Ala).

        Returns
        -------
        boolean
            Returns True if modification successful.

        """
        rxn = TER_ALA_LAC_SUB_RXN
        self.modifications["Alanine/Lactate Substitution"] = rxn.RunReactantInPlace(
            self.mol)
        if self.modifications["Alanine/Lactate Substitution"]:
            self.peptide.peptide_code[-1] = "Lac"
            self.refresh()
        return self.modifications["Alanine/Lactate Substitution"]

    def form_terminal_muramic_lactam(self):
        """
        Creates a modified PGN in-place with a muramic lactam (absent peptide).

        Returns
        -------
        boolean
            Returns True if modification successful.

        """
        # should this work with MurNGcl?
        rxn = GLYCAN_RXN_TYPES["Muramic Lactam"]
        self.modifications["Muramic Lactam"] = rxn.RunReactantInPlace(self.mol)
        if self.modifications["Muramic Lactam"]:
            self.glycan.gly_code[-1] += "[lactam]"
            self.peptide = None
            self.polymerisation_type = None
            self.linked_PGN = None
            self.refresh()
        return self.modifications["Muramic Lactam"]

    def check_endopeptidase(self):
        """
        Checks if PGN is an endopeptidase product (P1 or P2) and adds
        modification to self.

        Returns
        -------
        bool
            Returns True if endopeptidase product.

        """
        if self.peptide.len_null == 2:
            self.modifications["EPase P2"] = True
            return True
        elif self.peptide.len_null == 1:
            self.modifications["EPase P1"] = True
            return True
        else:
            return False


# %% Importing

AMINO_ACID_DF = pd.read_excel(
    DATABASE, sheet_name="AminoAcids", index_col=1)
GLYCAN_DF = pd.read_excel(
    DATABASE, sheet_name="Glycans", index_col=1)

###AAs###
NULL_AA = AminoAcid()
AMINO_ACID_DB = {}
AMINO_ACID_DIAMINO = AMINO_ACID_DF[AMINO_ACID_DF["SIDE_N_SMARTS"].notna(
)]["Code_Ext"].to_list()
AMINO_ACID_DICARBOXY = AMINO_ACID_DF[AMINO_ACID_DF["SIDE_C_SMARTS"].notna(
)]["Code_Ext"].to_list()
ALL_AMINO_ACID_LST = AMINO_ACID_DF[AMINO_ACID_DF["Type"]
                                   != "Dummy"]["Code_Ext"].tolist()
NICK_DB = {}
temp_dict = AMINO_ACID_DF.to_dict(orient="index")

for aa_code in temp_dict:
    smiles = temp_dict[aa_code]["SMILES"]
    aa_code_ext = temp_dict[aa_code]["Code_Ext"]  # uses 3 letter coding
    eupeptide_smarts = {"COOH":     temp_dict[aa_code]["EUPEPTIDE_C_SMARTS"],
                        "NH2":      temp_dict[aa_code]["EUPEPTIDE_N_SMARTS"]}
    side_chain_smarts = {"COOH":    temp_dict[aa_code]["SIDE_C_SMARTS"],
                         "NH2":      temp_dict[aa_code]["SIDE_N_SMARTS"]}
    AMINO_ACID_DB[aa_code_ext] = AminoAcid(smiles=smiles,
                                           peptide_code=aa_code_ext,
                                           eupeptide_smarts=eupeptide_smarts,
                                           side_chain_smarts=side_chain_smarts)
    NICK_DB[aa_code_ext] = aa_code


###Glycans###
temp_dict = GLYCAN_DF.to_dict(orient="index")
ALL_GLC_LST = GLYCAN_DF[GLYCAN_DF["Type"] == "Glc"]["Code_Ext"].tolist()
ALL_MUR_LST = GLYCAN_DF[GLYCAN_DF["Type"] == "Mur"]["Code_Ext"].tolist()
ALL_REDMUR_LST = GLYCAN_DF[GLYCAN_DF.index.str.contains(
    "[r]")]["Code_Ext"].tolist()
ALL_ANHYDROMUR_LST = GLYCAN_DF[GLYCAN_DF.index.str.contains(
    "an")]["Code_Ext"].tolist()

for gly_code in temp_dict:
    smiles = temp_dict[gly_code]["SMILES"]
    gly_type = temp_dict[gly_code]["Type"]
    gly_code_ext = temp_dict[gly_code]["Code_Ext"]  # uses 3 letter coding
    GLYCAN_DB[gly_code_ext] = Glycan(smiles=smiles,
                                     gly_code=gly_code_ext,
                                     init_gly_type=gly_type,
                                     fin_gly_type=gly_type)
    NICK_DB[gly_code_ext] = f"({gly_code})"

del temp_dict

TER_ALA_LAC_SUB_RXN = AllChem.ReactionFromSmarts(
    PEPTIDE_RXN_TYPES["Template N-O Sub"].replace(
        "*NH2*",
        AMINO_ACID_DB["Ala"].eupeptide_smarts["NH2"]))


def shorten_molecule_name(synonym):
    """
    Shortens molecule's synonym and returns molecule's name.
    Substitutes BRANCH_CHAR to shBRANCH_CHAR.

    Parameters
    ----------
    synonym : string
        Systematic name, long form

    Returns
    -------
    name : string
        Systematic name, short form


    """
    def replace(part):
        for symbol in ["[", "]", BRANCH_CHAR, ";","<", ">", "(", ")", ","]:
            if symbol in part:
                return symbol.join(
                    replace(subpart) for subpart in part.split(symbol))
        if part in NICK_DB:
            return NICK_DB[part]
        else:
            return part
    name_parts = [replace(part) for part in synonym.split(BOND_CHAR)]
    name = "".join(name_parts)
    name = name.replace(BRANCH_CHAR,shBRANCH_CHAR)
    return name
