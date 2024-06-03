'''
Custom Exceptions
env = chem
'''
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw  # drawing functions
from pathlib import Path

CWD = Path.cwd()

# %% Molecules


def get_substructs_from_rxn(rxn):
    return rxn.GetReactantTemplate(0), rxn.GetReactantTemplate(1)


def highlight_reacting_groups(cpd1, cpd2, type_rxn, rxn, error_name):
    substruct1, substruct2 = get_substructs_from_rxn(rxn)
    match1 = sum(cpd1.mol.GetSubstructMatches(substruct1), ())
    match2 = sum(cpd2.mol.GetSubstructMatches(substruct2), ())
    img1 = Draw.MolToImage(cpd1.mol, size=(800, 800), highlightAtoms=match1)
    img1.save(CWD/"img"/f"{error_name}_cpd1.png")
    img2 = Draw.MolToImage(cpd2.mol, size=(800, 800), highlightAtoms=match2)
    img2.save(CWD/"img"/f"{error_name}_cpd2.png")


class NoReactionError(Exception):

    def __init__(self, cpd1, cpd2, type_rxn, rxn=None):
        """
        Exception raised when mo products are formed.

        Parameters
        ----------
        cpd1 : Molecule
        cpd2 : Molecule
        type_rxn : string
            Name of reaction.
        rxn : Chem.Reaction, optional
            Chem.Reaction object. The default is None.

        Returns
        -------
        None.

        """
        message = f"No {type_rxn} reaction."
        cpd1_message = f"\nself: \t{cpd1}, {cpd1.smiles}"
        cpd2_message = f"\nother: \t{cpd2}, {cpd2.smiles}"
        self.message = message+cpd1_message+cpd2_message
        if type_rxn != "invalid" and rxn != None:
            highlight_reacting_groups(cpd1, cpd2, type_rxn, rxn, "NoReaction")
        super().__init__(self.message)


class TooManyPdtsError(Exception):

    def __init__(self, cpd1, cpd2, type_rxn, rxn=None, pdt_smiles=[]):
        """
        Exception raised when too many products are formed.

        Parameters
        ----------
        cpd1 : Molecule
        cpd2 : Molecule
        type_rxn : string
            Name of reaction.
        rxn : Chem.Reaction, optional
            Chem.Reaction object. The default is None.
        pdt_smiles : list, optional
            List of product SMILES. The default is [].

        Returns
        -------
        None.

        """
        message = f"Too many pdts for {type_rxn} reaction."
        cpd1_message = f"\nself: \t{cpd1}, {cpd1.smiles}"
        cpd2_message = f"\nother: \t{cpd2}, {cpd2.smiles}"
        pdt_message = "".join(f"\n{i}:\t{smiles}" for i,
                              smiles in enumerate(pdt_smiles))
        self.message = message+cpd1_message+cpd2_message+pdt_message
        if type_rxn != "invalid" and rxn != None:
            highlight_reacting_groups(cpd1, cpd2, type_rxn, rxn, "TooManyPdts")
        super().__init__(self.message)


# %% Generator

class InputError(Exception):

    def __init__(self,component,msg):
        self.message = f"Invalid Generator input for {component}.\n{msg}"
        super().__init__(self.message)

# %% Illustrator

# Not necessary

# %% Fragmenter

# to add

# %% MSPMaker

# to add

# %% UI

# to add
