'''
Checks bonding for AA, Peptide, Glycan
env = chem
'''

# %% Load Molecules

import yaml

from Common import flatten
from Molecules import Molecule
from Molecules import AMINO_ACID_DB, GLYCAN_DB
from Reactions import PEPTIDE_FRAG_TYPES, GLYCAN_FRAG_TYPES
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt
from rdkit import Chem
from rdkit.Chem import Draw  # drawing functions
from rdkit.Chem.Descriptors import MolLogP
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.molSize = 600,400

Lac = Molecule("OC(C)C(O)=O")
A = AMINO_ACID_DB["Ala"]  # Alanine
E = AMINO_ACID_DB["Glu"]  # Glutamic Acid
Q = AMINO_ACID_DB["Gln"]  # Glutamine
R = AMINO_ACID_DB["Arg"]  # Arginine
K = AMINO_ACID_DB["Lys"]  # Lysine
F = AMINO_ACID_DB["Phe"]  # Phenylalanine
G = AMINO_ACID_DB["Gly"]  # Glycine
q = AMINO_ACID_DB["isoGln"]  # Isoglutamine
e = AMINO_ACID_DB["isoGlu"]  # Isoglutamic Acid
m = AMINO_ACID_DB["mDAP"]  # mDAP
a = AMINO_ACID_DB["mDAP(NH2)"]  # Amidated mDAP

C = AMINO_ACID_DB["Cys"]  # Cysteine
D = AMINO_ACID_DB["Asp"]  # Aspartic Acid
H = AMINO_ACID_DB["His"]  # Histidine
W = AMINO_ACID_DB["Trp"]  # Tryptophan
o = AMINO_ACID_DB["Orn"]  # Ornithine
d = AMINO_ACID_DB["isoAsp"]  # Isoaspartic Acid[?]
s = AMINO_ACID_DB["Hse"]  # Homoserine
t = AMINO_ACID_DB["Lan"]  # Lanthionine

GLU = GLYCAN_DB["GlcNAc"]
MUR = GLYCAN_DB["MurNAc"]
ANH = GLYCAN_DB["anMurNAc"]
MUR_N = GLYCAN_DB["MurN"]
MUR_GCL = GLYCAN_DB["MurNGlyc"]
MUR_DIAC = GLYCAN_DB["MurNAcOAc"]
GLU_N = GLYCAN_DB["GlcN"]

GLU_MUR = GLU+MUR
GLU_ANH = GLU+ANH

# %% Load Test SMILES
SMILES_FILE = "data/check_cpds.yaml"
with open(SMILES_FILE,"r") as file:
    SMILES_CHECK = yaml.safe_load(file)

FRAG_MZ_FILE = "data/check_frag.yaml"
with open(FRAG_MZ_FILE,"r") as file:
    MZ_CHECK = yaml.safe_load(file)

# %% Test Functions - Bond Formation

def test_peptide_bond_formation():
    '''
    Test peptide bond formation.
    Sequences: 
    [1] Ala.β-Asp.Orn.Ala.Ala (3–NH2–Ala.Ala)
    [2] Ala.Gly.mDAP(NH2).Trp (3–NH2–Ala.Gly.Lys.Ala)
    [3] Ala.isoGln.Lys.His (3–NH2–Gly)
    [4] Gly.isoGlu.Hse.Ala (2–COOH–Orn)
    [5] Ala.isoGlu.mDAP.Hse (2–COOH–Ala,3–NH2–Cys)

    Returns
    -------
    peptide_str : List
        List of peptides.

    '''

    # Test1 Orn as DAA, Isoaspartic Acid
    test1_lateral = A+A
    test1_tri_w_lateral = test1_lateral - (Lac+A+d+o)
    test1_full = (test1_tri_w_lateral+A)+A

    # Test2 Amidated mDAP as DAA, AGKA lateral chain
    test2_lateral = ((A+K)+G)+A  # has to be built H2N ---> COOH
    test2_tri_w_lateral = test2_lateral - (Lac+A+G+a)
    test2_full = test2_tri_w_lateral + W

    # Test3 Lys as DAA, with His as AA4, Gly lateral chain
    test3_tri_w_lateral = G - (Lac+A+q+K)
    test3_full = test3_tri_w_lateral + H

    # Test 4 Branching at isoGlu (AA2) with Orn lateral
    test4_di_w_lateral = (Lac+G+e) & o
    test4_full = (test4_di_w_lateral + s)+A

    # Test 5 Branching Asp-COOH, mDAP-NH2 with laterals A,C resp.
    test5_full = (C - (((A+e) & A)+m))+s

    peptide_str = [test1_full, test2_full,
                   test3_full, test4_full, test5_full]
    peptide_str_smiles = [x.smiles for x in peptide_str]

    print("\nChecking with various peptides [1] to [5]...")
    for i, smiles in enumerate(peptide_str_smiles):
        test_bool = smiles == SMILES_CHECK["peptide"][i]
        print(f"\nTest peptide [{i+1}] \n{peptide_str[i].name} \nSMILES: \t{test_bool}")
        if not test_bool:
            print(f"PREV: {SMILES_CHECK['peptide'][i]} \nNOW: {smiles}")

    return peptide_str


def test_glycosidic_bond_formation():
    '''
    Test glycosidic bond formation.
    Sequences:    
    [1] GlcNAc.MurNAc
    [2] GlcNAc.anMurNAc
    [3] MurNAcOAc.MurNAc.GlcNAc.anMurNAc
    [4] GlcN.MurNGlyc
    [5] MurN.MurNAcOAc
    
    Returns
    -------
    glycan_str : List
        List of glycans.

    '''

    test1 = GLU+MUR  # standard
    test2 = GLU+ANH  # anhydro
    test3 = MUR_DIAC+MUR+GLU+ANH  # tetra saccharide
    test4 = GLU_N+MUR_GCL
    test5 = MUR_N+MUR_DIAC

    glycan_str = [test1, test2, test3, test4, test5]
    glycan_str_smiles = [x.smiles for x in glycan_str]

    print("\nChecking with various glycans [1] to [5]...")
    for i, smiles in enumerate(glycan_str_smiles):
        test_bool = smiles == SMILES_CHECK["glycan"][i]
        print(f"\nTest glycan [{i+1}] \n{glycan_str[i].name} \nSMILES: \t{test_bool}")
        if not test_bool:
            print(f"PREV: {SMILES_CHECK['glycan'][i]} \nNOW: {smiles}")
    return glycan_str


def test_peptidoglycan_formation():
    '''
    Test peptidoglycan formation.
    Sequences:
    [1] GlcNAc.anMurNAc-Ala.isoGlu.mDAP.Ala
    [2] anMurNAc-Ala.isoGlu.mDAP.Ala
    [3] GlcNAc.MurNAc-Ala.isoGlu.mDAP.Ala
    [4] GlcNAc.MurNAc-Ala.isoGln.Lys.Ala.Ala[3-NH2-Gly.Gly.Gly.Gly.Gly]
    [5] MurNAc-Ala.isoGln.Lys.Ala.Ala[3-NH2-Gly.Gly.Gly.Gly.Gly]
    [6] GlcNAc.MurNAc-Ala.isoGln.Orn.Ala[3-NH2-Ala.Orn.isoGln.Ala]
    [7] GlcNAc.MurNAc-Ala.isoGlu.Lan.Ala[3-NH2-Lan.isoGlu.Ala]

    Returns
    -------
    PGN_str : List
        List of peptidoglycans.

    '''
    # Tan et al structure [a-c]
    a_tetra = Lac+A+e+m+A

    Tan_str_a = GLU_ANH & a_tetra
    Tan_str_b = ANH & a_tetra
    Tan_str_c = GLU_MUR & a_tetra

    # Tan et al structure [d-e]
    d_tri = Lac+A+q+K
    d_side = G+G+G+G+G
    d_full = (d_side-d_tri)+A+A

    Tan_str_d = GLU_MUR & d_full
    Tan_str_e = MUR & d_full
    
    # Ornithine as DAA (Bifidobacterium)
    peptide_6a = A+q+o+A
    peptide_6 = peptide_6a-(Lac+A+q+o)+A #3-4 cross-link
    PGN_6 = GLU_MUR & peptide_6
    
    # Lanthione as DAA (Fusobacterium)
    peptide_7a = A+e+t
    peptide_7 = (peptide_7a - (Lac+A+e+t))+A #3-3 cross-link
    PGN_7 = GLU_MUR & peptide_7 
    

    PGN_str = [Tan_str_a, Tan_str_b, Tan_str_c, Tan_str_d, Tan_str_e,
               PGN_6, PGN_7]
    PGN_str_smiles = [x.smiles for x in PGN_str]

    print("\nChecking with various PGN [1] to [7]...")
    for i, smiles in enumerate(PGN_str_smiles):
        test_bool = smiles == SMILES_CHECK['peptidoglycan'][i]
        print(f"\nTest PGN [{i+1}] \n{PGN_str[i].name} \nSMILES: \t{test_bool}")
        if not test_bool:
            print(
                f"PREV: {SMILES_CHECK['peptidoglycan'][i]} \nNOW: {smiles}")

    return PGN_str

# %% Test Functions - Bond Fragmentation

def test_peptide_fragmentation():
    print("Not coded yet!")
    pass

def test_glycan_fragmentation():
    '''
    Test glycan fragmentation.
    Sequences:
    [1] GlcNAc.MurNAc-Ala
    [2] GlcNAc.MurNAc[r]-Ala
    [3] GlcNAc.anMurNAc-Ala
    [4] GlcN.MurNGlyc
    [5] MurN.MurNAcOAc
    [6] MurNAc
    [7] MurN

    Returns
    -------
    results : dict
        Dict of fragment m/z.

    '''
    
    test1 = GLU+MUR&(Lac+A)
    test2 = GLU+MUR&(Lac+A)
    test2.reduce_terminal_muramic_acid()
    test3 = GLU+ANH&(Lac+A)
    test4 = GLU_N+MUR_GCL
    test5 = MUR_N+MUR_DIAC
    test6 = MUR
    test7 = MUR_N
    
    results = {}
    molecules = [test1,test2,test3,test4,test5,test6,test7]
    for molecule in molecules:
        results[molecule.smiles] = {}
        for frag_type in GLYCAN_FRAG_TYPES:
            frag_rxn = GLYCAN_FRAG_TYPES[frag_type]
            pdts = list(flatten(frag_rxn.RunReactants((molecule.mol,))))
            for mol in pdts:
                Chem.SanitizeMol(mol)
            mw = [round(CalcExactMolWt(mol),3) for mol in pdts]                
            results[molecule.smiles][frag_type]=mw
            
    print("\nChecking fragmentation with various glycans [1] to [7]...")
    for i, m in enumerate(molecules):
        test_bool = results[m.smiles] == MZ_CHECK['glycan'][m.smiles]
        print(f"\nTest Glycan [{i+1}] \n{m.name}; {m.mMass} \nMatch: \t{test_bool}")
        if not test_bool:
            print("\nPREV:")
            for rxn, fragments in MZ_CHECK['glycan'][m.smiles].items():
                print(f"{rxn}:\t{fragments}")
            print("\nNOW:")   
            for rxn, fragments in results[m.smiles].items():
                print(f"{rxn}:\t{fragments}")      
    return results

# peptide_test = test_peptide_bond_formation()
# glycan_test = test_glycosidic_bond_formation()
# pgn_test = test_peptidoglycan_formation()

# peptide_frags = test_peptide_fragmentation()
glycan_frags = test_glycan_fragmentation()
