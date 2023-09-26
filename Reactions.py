'''
Reactions and SMARTS
env = chem
'''
from rdkit.Chem.Draw import rdDepictor
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.molSize = 600, 400
rdDepictor.SetPreferCoordGen(True)


def make_rxn(smarts):
    return AllChem.ReactionFromSmarts(smarts)

# %% Bonding Reactions

# %%% Peptides

# %%%% Lactate Removal

PEP_RXN_LAC_REMOVAL_SMARTS = "[CH3X4:99]-[CHX4:98](-[OH:97])-[CX3:96](=[O:95])-[NX3:94]>>[NX3:94]"
PEP_RXN_LAC_REMOVAL = make_rxn(PEP_RXN_LAC_REMOVAL_SMARTS)

# %%%% Peptide Bond

PEP_RXN_PEPTIDE_BOND_SMARTS = "[C:1](=[O:2])-[OD1].[NX3;H2,H1;!$(NC=O):50]>>[C:1](=[O:2])[N:50]"
PEP_RXN_PEPTIDE_BOND = make_rxn(PEP_RXN_PEPTIDE_BOND_SMARTS)

# %%%% Peptide Templates

PEP_RXNTEMP_EUPEPTIDE_SMARTS = "*COOH*-[OH].[NX3;H2,H1;!$(NC=O):50]-*NH2*>>*COOH*-[NX3;!$(NC=O):50]-*NH2*"
PEP_RXNTEMP_EUPEPTIDE_NH2_SMARTS = "[C:1](=[O:2])-[OD1].[NX3;H2,H1;!$(NC=O):50]-*NH2*>>[C:1](=[O:2])-[NX3;!$(NC=O):50]-*NH2*"
PEP_RXNTEMP_EUPEPTIDE_COOH_SMARTS = "*COOH*-[OH].[NX3;H2,H1;!$(NC=O):50]>>*COOH*-[NX3;!$(NC=O):50]"

# terminal only
PEP_RXNTEMP_N_O_SUB_SMARTS = "[ND2:50]-*NH2*-[OH:51]>>[OX2:50]-*NH2*-[OH:51]"

# %%%% Rxn Types

PEPTIDE_RXN_TYPES = {"Lac Removal": PEP_RXN_LAC_REMOVAL,
                     "Peptide": PEP_RXN_PEPTIDE_BOND,
                     "Template Eupeptide": PEP_RXNTEMP_EUPEPTIDE_SMARTS,
                     "Template Eupeptide NH2": PEP_RXNTEMP_EUPEPTIDE_NH2_SMARTS,
                     "Template Eupeptide COOH": PEP_RXNTEMP_EUPEPTIDE_COOH_SMARTS,
                     "Template N-O Sub": PEP_RXNTEMP_N_O_SUB_SMARTS}

# %%% Glycans

# [OH:25]
GLY_DONOR = "[CH1X4:21]1(-[C:23]-[OX2:24])-[CX4:20]-[CX4:19]-[CH1X4:18](-[NX3:26])-[CX4:17](-[O:22]-1)-[OH:25]"
# [OH:16]
GLY_ACCEPTOR = "[OH:16]-[CH1X4:4](-[CX4:5]-[C:7]-[O:8])-[CH1X4:3](-[O:10])-[CX4:2](-[NX3:9])-[CX4;$([C]-[OX2]),CR2:1]"
GLY_ACCEPTED = "[OX2H0:25]-[CH1X4:4](-[CX4:5]-[C:7]-[O:8])-[CH1X4:3](-[O:10])-[CX4:2](-[NX3:9])-[CX4;$([C]-[OX2]),CR2:1]"
GLY_COMB = GLY_ACCEPTOR.replace("[OH:16]", GLY_DONOR)
GLY_COMB = GLY_COMB.replace("[OH:25]", "[O:25]")
DIGLY_COMB = '[CH1X4:21]1(-[C:23]-[OX2:24])-[CX4:20](-[OX2:26]-[CH1X4:27]2-[CH1X4:28](-[NX3:35])-[CH1X4:29]-[CH1X4:30]-[CH1X4:31](-[C:33]-[OX2:34])-[OX2:32]-2)-[CX4:19]-[CH1X4:18](-[NX3:26])-[CX4:17](-[O:22]-1)-[O:25]-[CH1X4:4](-[CX4:5]-[C:7]-[O:8])-[CH1X4:3](-[O:10])-[CX4:2](-[NX3:9])-[CX4;$([C]-[OX2]),CR2:1]'

# %%%% Glycosidic Bond

GLY_RXN_GLYCOSIDIC_BOND_SMARTS = f"{GLY_DONOR}.{GLY_ACCEPTOR}>>{GLY_COMB}"
GLY_RXN_GLYCOSIDIC_BOND = make_rxn(GLY_RXN_GLYCOSIDIC_BOND_SMARTS)

# %%%% Muramic Reduction

MURAMIC_RED_0 = "[OX2:16]-[CH1X4:4]2-[CH1X4:3]-[CX4:2](-[NX3:9])-[CH1X4:1](-[OH:27])-[O:6]-[CX4:5](-[C:7]-[O:8])-2"
MURAMIC_RED_1 = "[OH:6]-[CX4:5](-[C:7]-[O:8])-[CH1X4:4](-[OX2:16])-[CH1X4:3]-[CX4:2](-[NX3:9])-[CH2X4:1]-[OH:27]"
GLY_RXN_MURAMIC_RED_SMARTS = f"{MURAMIC_RED_0}>>{MURAMIC_RED_1}"
GLY_RXN_MURAMIC_RED = make_rxn(GLY_RXN_MURAMIC_RED_SMARTS)

# %%%% Muramic Lactamization

MURAMIC_LACTAM_0 = "[CH3:1]-[C:2](=[O:3])-[NH:4]-[CX4:5](-[CH1X4;$([C]-[OH]),CR2:6])-[CX4:7]-[OD2:8]-[CH:9]-[CH3:10]"
MURAMIC_LACTAM_1 = "[CH2:1]1-[C:2](=[O:3])-[NH:4]-[CX4:5](-[CH1X4;$([C]-[OH]),CR2:6])-[CX4:7]-[OD2:8]-1"
GLY_RXN_MURAMIC_LACTAM_SMARTS = f"{MURAMIC_LACTAM_0}>>{MURAMIC_LACTAM_1}"
GLY_RXN_MURAMIC_LACTAM = make_rxn(GLY_RXN_MURAMIC_LACTAM_SMARTS)

# %%%% Lactate Merger (Forming Glycopeptide)

LAC_MERGER_GLYCAN = "[CX4;$([C]-[OH]),CR2:1]-[CH1X4;$(CN):2]-[CH1X4:3]-[O:4]-[CH1X4:5](-[CH3:6])-[CX3:7](=[O:8])-[OH:9]"
LAC_MERGER_PEPTIDE = "[CH3X4:99]-[CHX4:98](-[OX2:97])-[CX3:96](=[O:95])-[NX3:94]"
LAC_MERGER_GLYCOPEPTIDE = "[CX4;$([C]-[OX2]),CR2:1]-[CH1X4;$(CN):2]-[CH1X4:3]-[O:4]-[CH1X4:5](-[CH3:6])-[CX3:7](=[O:8])-[NX3:94]"

GLY_RXN_LAC_MERGER_SMARTS = f"{LAC_MERGER_GLYCAN}.{LAC_MERGER_PEPTIDE}>>{LAC_MERGER_GLYCOPEPTIDE}"
GLY_RXN_LAC_MERGER = make_rxn(GLY_RXN_LAC_MERGER_SMARTS)

# %%%% Rxn Types

GLYCAN_RXN_TYPES = {"Glycosidic": GLY_RXN_GLYCOSIDIC_BOND,
                    "Lac Merger": GLY_RXN_LAC_MERGER,
                    "Muramic Reduction": GLY_RXN_MURAMIC_RED,
                    "Muramic Lactam": GLY_RXN_MURAMIC_LACTAM}


# %% Molecular Patterns (SMARTS)

# %%% Peptides

# both side chains True
PEP_PATT_PEPTIDE_SMARTS = "[C;!H3][CX3](=[OX1])-[NX3]-[#6]"
PEP_PATT_PEPTIDE = Chem.MolFromSmarts(PEP_PATT_PEPTIDE_SMARTS)

# both side chains True
PEP_PATT_EUPEPTIDE_SMARTS = "[NX3]-"+PEP_PATT_PEPTIDE_SMARTS+"[CX3](=[OX1])"
PEP_PATT_EUPEPTIDE = Chem.MolFromSmarts(PEP_PATT_EUPEPTIDE_SMARTS)

# COOH_side_chain False
PEP_PATT_ALPHA_C_EUPEPTIDE_SMARTS = "[#6][CX3](=[OX1])-" + \
    PEP_PATT_EUPEPTIDE_SMARTS
PEP_PATT_ALPHA_C_EUPEPTIDE = Chem.MolFromSmarts(
    PEP_PATT_ALPHA_C_EUPEPTIDE_SMARTS)

# NH2_side_chain False
PEP_PATT_ALPHA_N_EUPEPTIDE_SMARTS = PEP_PATT_EUPEPTIDE_SMARTS+"-[NX3]-[#6]"
PEP_PATT_ALPHA_N_EUPEPTIDE = Chem.MolFromSmarts(
    PEP_PATT_ALPHA_N_EUPEPTIDE_SMARTS)

# di eupeptide
PEP_PATT_DIPEPTIDE_SMARTS = PEP_PATT_EUPEPTIDE_SMARTS + \
    "[NX3]-[#6]-[CX3](=[OX1])"
PEP_PATT_DIPEPTIDE = Chem.MolFromSmarts(PEP_PATT_DIPEPTIDE_SMARTS)

# tri eupeptide
PEP_PATT_TRIPEPTIDE_SMARTS = PEP_PATT_DIPEPTIDE_SMARTS + \
    "[NX3]-[#6]-[CX3](=[OX1])"
PEP_PATT_TRIPEPTIDE = Chem.MolFromSmarts(PEP_PATT_TRIPEPTIDE_SMARTS)

# tetra eupeptide
PEP_PATT_TETRAPEPTIDE_SMARTS = PEP_PATT_TRIPEPTIDE_SMARTS + \
    "[NX3]-[#6]-[CX3](=[OX1])"
PEP_PATT_TETRAPEPTIDE = Chem.MolFromSmarts(PEP_PATT_TETRAPEPTIDE_SMARTS)

# di peptide lax
PEP_PATT_LAX_DIPEPTIDE_SMARTS = PEP_PATT_DIPEPTIDE_SMARTS[6:-14]
PEP_PATT_LAX_DIPEPTIDE = Chem.MolFromSmarts(PEP_PATT_LAX_DIPEPTIDE_SMARTS)

# tri peptide lax
PEP_PATT_LAX_TRIPEPTIDE_SMARTS = PEP_PATT_TRIPEPTIDE_SMARTS[6:-14]
PEP_PATT_LAX_TRIPEPTIDE = Chem.MolFromSmarts(PEP_PATT_LAX_TRIPEPTIDE_SMARTS)

# tetra peptide lax
PEP_PATT_LAX_TETRAPEPTIDE_SMARTS = PEP_PATT_TETRAPEPTIDE_SMARTS[6:-14]
PEP_PATT_LAX_TETRAPEPTIDE = Chem.MolFromSmarts(
    PEP_PATT_LAX_TETRAPEPTIDE_SMARTS)

# isoGln/isoGlu
PEP_PATT_ISOGLX_DIPEPTIDE_SMARTS = "[NX3]-[#6](-[CX3](=[OX1]))-[#6]-[#6]-[CX3](=[OX1])-[NX3]-[#6]-[CX3](=[OX1])-[NX3]"
PEP_PATT_ISOGLX_DIPEPTIDE = Chem.MolFromSmarts(
    PEP_PATT_ISOGLX_DIPEPTIDE_SMARTS)

PEPTIDE_PATT_TYPES = {"Peptide": PEP_PATT_PEPTIDE,
                      "Alpha C Eupeptide": PEP_PATT_ALPHA_C_EUPEPTIDE,
                      "Alpha N Eupeptide": PEP_PATT_ALPHA_N_EUPEPTIDE,
                      "Eupeptide": PEP_PATT_EUPEPTIDE,
                      "Dipeptide": PEP_PATT_DIPEPTIDE,
                      "Tripeptide": PEP_PATT_TRIPEPTIDE,
                      "Tetrapeptide": PEP_PATT_TETRAPEPTIDE,
                      "Lax Dipeptide": PEP_PATT_LAX_DIPEPTIDE,
                      "Lax Tripeptide": PEP_PATT_LAX_TRIPEPTIDE,
                      "Lax Tetrapeptide": PEP_PATT_LAX_TETRAPEPTIDE,
                      "isoGlx Dipeptide": PEP_PATT_ISOGLX_DIPEPTIDE}

# %%% Misc.

# Gobbi, A. & Poppinger, D. “Genetic optimization of combinatorial libraries.” Biotechnology and Bioengineering 61, 47-54 (1998)

MISC_PATT_BASIC = Chem.MolFromSmarts(
    "[#7;+,$([N;H2&+0][$([C,a]);!$([C,a](=O))]),$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]")
MISC_PATT_ACIDIC = Chem.MolFromSmarts("[$([C,S](=[O,S,P])-[O;H1,-1])]")
MISC_PATT_HYDROXYL = Chem.MolFromSmarts("[CX4]-[OH1]")
MISC_PATT_AMINO = Chem.MolFromSmarts("[CX4]-[NH2]")
MISC_PATT_GLY_BOND = Chem.MolFromSmarts(GLY_COMB)
MISC_PATT_DIGLY_BOND =  Chem.MolFromSmarts(DIGLY_COMB)
MISC_PATT_GLY_DONOR = Chem.MolFromSmarts(GLY_DONOR)
MISC_PATT_GLY_ACCEPTOR = Chem.MolFromSmarts(GLY_ACCEPTOR)
MISC_PATT_LACTOYL_BOND = Chem.MolFromSmarts("[#6]-[OX2]-[CHX4](-[CH3])-[CX3](=[OX1])-[NX3]")
MISC_PATT_LACTATE = Chem.MolFromSmarts("[OH]-[CHX4](-[CH3])-[CX3](=[OX1])-[NX3]")

MISC_PATT_TYPES = {"Peptide": PEP_PATT_PEPTIDE,
                   "Gly Bond": MISC_PATT_GLY_BOND,
                   "DiGly Bond": MISC_PATT_DIGLY_BOND,
                   "Gly Donor": MISC_PATT_GLY_DONOR,
                   "Gly Acceptor": MISC_PATT_GLY_ACCEPTOR,
                   "Basic": MISC_PATT_BASIC,
                   "Acidic": MISC_PATT_ACIDIC,
                   "Hydroxyl": MISC_PATT_HYDROXYL,
                   "Amino": MISC_PATT_AMINO,
                   "Lactoyl Bond": MISC_PATT_LACTOYL_BOND,
                   "Lactate": MISC_PATT_LACTATE,}

# %% Fragmentations

# %%% Peptides

# Should not cleave Ac or create NH3+
# KIV have different B,Y reactions corresponding to different env.

# %%%% Peptide by

PEP_FRAG_Y_SMARTS = "[C!H3;!$(C=O):4]-[C:1](=[O:2])[NX3!H2;!$(NC=C):3]>>[NX4+:3].[O+]"
PEP_FRAG_Y = make_rxn(PEP_FRAG_Y_SMARTS)

PEP_FRAG_B_SMARTS = "[C!H3;!$(C=O):4]-[C:1](=[O:2])[NX3!H2;!$(NC=C):3]>>[C:4]-[C:1](#[O+:2]).[O+]"
PEP_FRAG_B = make_rxn(PEP_FRAG_B_SMARTS)

# %%%% AA specific reactions

# isoGln-1 and isoGln-2
ISOGLN_0 = "[CH2X4:1](-[CHX4:2](-[CX3:3](=[O:4])-[NH2X3:5])-[NH2:6])-[CH2X4:7]-[CX3:8]=[O:10]"
ISOGLN_1 = "[CH2X4:1](-[CHX4:2](-[CX3:3](#[O+:4]))-[NH2:6])-[CH2X4:7]-[CX3:8]=[O:10]"
ISOGLN_2 = "[CX2:1](#[CX2:2])-[CH2X4:7]-[CX3:8]=[O+:10]"
PEP_FRAG_ISOGLN_SMARTS = f"{ISOGLN_0}>>{ISOGLN_1}.{ISOGLN_2}"
PEP_FRAG_ISOGLN = make_rxn(PEP_FRAG_ISOGLN_SMARTS)

# isoGlu-1 and isoGlu-2
ISOGLU_0 = "[CH2X4:1](-[CHX4:2](-[CX3:3](=[O:4])-[OHX2:5])-[NH2:6])-[CH2X4:7]-[CX3:8]=[O:10]"
ISOGLU_1 = "[CH2X4:1](-[CHX4:2](-[CX3:3](#[O+:4]))-[NH2:6])-[CH2X4:7]-[CX3:8]=[O:10]"
ISOGLU_2 = "[CX2:1](#[CX2:2])-[CH2X4:7]-[CX3:8]=[O+:10]"
PEP_FRAG_ISOGLU_SMARTS = f"{ISOGLU_0}>>{ISOGLU_1}.{ISOGLU_2}"
PEP_FRAG_ISOGLU = make_rxn(PEP_FRAG_ISOGLU_SMARTS)

# Lys-1
LYS_0 = "[CH1X4:1](-[NH2:9])(-[CH2:2]-[CH2:3]-[CH2:4]-[CH2:5]-[NX3:6])-[C:7](=[O:8])"
LYS_1 = "[CX2:1](#[C:2]-[CH2:3]-[CH2:4]-[CH2:5]-[NX4+:6])"
PEP_FRAG_LYS_SMARTS = f"{LYS_0}>>{LYS_1}.[O+]"
PEP_FRAG_LYS = make_rxn(PEP_FRAG_LYS_SMARTS)

# %%%% Dehydration of COOH/CONH2

PEP_FRAG_ACYL_N_SMARTS = "[#6:1]-[CX3:2](=[OX1:3])-[NH2:4]>>[#6:1]-[CX3:2](#[OX1+:3]).[N+]"
PEP_FRAG_ACYL_N = make_rxn(PEP_FRAG_ACYL_N_SMARTS)

PEP_FRAG_ACYL_O_SMARTS = "[#6:1]-[CX3:2](=[OX1:3])-[OH1:4]>>[#6:1]-[CX3:2](#[OX1+:3]).[O+]"
PEP_FRAG_ACYL_O = make_rxn(PEP_FRAG_ACYL_O_SMARTS)

PEP_FRAG_r5_ACYL_N_SMARTS = "[CX3:1](=[O:2])(-[NH2:7])-[CX4:3]-[NH1X3:4]-[CX3:5]=[O:6]>>[CX3:1]1(=[O:2])-[CX4:3]-[NH1X3+:4]=[CX3:5]-[O:6]-1.[NH3+:7]"
PEP_FRAG_r5_ACYL_N = make_rxn(PEP_FRAG_r5_ACYL_N_SMARTS)

PEP_FRAG_r5_ACYL_O_SMARTS = "[CX3:1](=[O:2])(-[OH:7])-[CX4:3]-[NH1X3:4]-[CX3:5]=[O:6]>>[CX3:1]1(=[O:2])-[CX4:3]-[NH1X3+:4]=[CX3:5]-[O:6]-1.[OH2+:7]"
PEP_FRAG_r5_ACYL_O = make_rxn(PEP_FRAG_r5_ACYL_O_SMARTS)

# %%%% Frag Types

PEPTIDE_FRAG_TYPES = {"Peptide y": PEP_FRAG_Y,
                      "Peptide b": PEP_FRAG_B,
                      "-H2O": PEP_FRAG_ACYL_O,
                      "-NH3": PEP_FRAG_ACYL_N,
                      "-H2O(r5)": PEP_FRAG_r5_ACYL_O,
                      "-NH3(r5)": PEP_FRAG_r5_ACYL_N,
                      "Peptide q": PEP_FRAG_ISOGLN,
                      "Peptide e": PEP_FRAG_ISOGLU,
                      "Peptide K": PEP_FRAG_LYS}

# %%% Glycans

NONRED_GLY_COMB = GLY_COMB.replace(
    "[CX4;$([C]-[OX2]),CR2:1]",
    "[CX4;$([C]-[OX2H0]),CR2:1]")
RED_GLY_COMB = GLY_COMB.replace(
    "[CX4;$([C]-[OX2]),CR2:1]",
    "[CH2X4;$([C]-[OX2]):1]")

GLY_DONOR_B_FRAG = GLY_DONOR.replace(
    "[CH1X4:18](-[NX3:26])-[CX4:17](-[O:22]-1)-[OH:25]",
    "[CH1X4:18](-[NX3:26])-[CX4:17](=[O+:22]-1)")
GLY_ACCEPTOR_Y_FRAG = GLY_ACCEPTED.replace("[OX2H0:25]",
                                           "[OH2+:25]")
NONRED_GLY_ACCEPTOR_Y_FRAG = GLY_ACCEPTOR_Y_FRAG.replace(
    "[CX4;$([C]-[OX2]),CR2:1]",
    "[CX4;$([C]-[OX2H0]),CR2:1]")
RED_GLY_ACCEPTOR_Y_FRAG = GLY_ACCEPTOR_Y_FRAG.replace(
    "[CX4;$([C]-[OX2]),CR2:1]",
    "[CH2X4;$([C]-[OX2]):1]")

GLY_DONOR_C_FRAG = GLY_DONOR.replace("[OH:25]", "[OH2+:25]")

GLY_ACCEPTOR_Z_FRAG = GLY_ACCEPTED.replace(
    "[OX2H0:25]-[CH1X4:4](-[CX4:5]-[C:7]-[O:8])-[CH1X4:3](-[O:10])",
    "[CH1X3:4](-[CX4:5]-[C:7]-[O:8])=[CH0X3:3](-[O+:10])")
NONRED_GLY_ACCEPTOR_Z_FRAG = GLY_ACCEPTOR_Z_FRAG.replace(
    "[CX4;$([C]-[OX2]),CR2:1]",
    "[CX4;$([C]-[OX2H0]),CR2:1]")
RED_GLY_ACCEPTOR_Z_FRAG = GLY_ACCEPTOR_Z_FRAG.replace(
    "[CX4;$([C]-[OX2]),CR2:1]",
    "[CH2X4;$([C]-[OX2]):1]")

# %%%% Glycan BY

# Donors, Acceptors
NONRED_GLY_FRAG_BY_SMARTS = f"{NONRED_GLY_COMB}>>{GLY_DONOR_B_FRAG}.{NONRED_GLY_ACCEPTOR_Y_FRAG}"
NONRED_GLY_FRAG_BY = make_rxn(NONRED_GLY_FRAG_BY_SMARTS)
RED_GLY_FRAG_BY_SMARTS = f"{RED_GLY_COMB}>>{GLY_DONOR_B_FRAG}.{RED_GLY_ACCEPTOR_Y_FRAG}"
RED_GLY_FRAG_BY = make_rxn(RED_GLY_FRAG_BY_SMARTS)

# %%%% Glycan ZC

# Acceptor Dehydrated, H3O as dummy
NONRED_GLY_FRAG_CZ_SMARTS = f"{NONRED_GLY_COMB}>>{GLY_DONOR_C_FRAG}.{NONRED_GLY_ACCEPTOR_Z_FRAG}"
NONRED_GLY_FRAG_CZ = make_rxn(NONRED_GLY_FRAG_CZ_SMARTS)
RED_GLY_FRAG_CZ_SMARTS = f"{RED_GLY_COMB}>>{GLY_DONOR_C_FRAG}.{RED_GLY_ACCEPTOR_Z_FRAG}"
RED_GLY_FRAG_CZ = make_rxn(RED_GLY_FRAG_CZ_SMARTS)

# %%%% Lac Moiety Cleavage

# Cleavage about PEP: PEP-Glycan bond, PEP-Peptide bond
# More prominent in non-reduced Mur
FRAG_LAC_DEHYD = "[CX3:5](=[CH2:6])-[CX3:7](=[OH+:8])-[NX3:94]"
FRAG_LAC_PEPTIDE = "[OH2+:4]-[CH1X4:5](-[CH3:6])-[CX3:7](=[O:8])-[NX3:94]"
FRAG_LAC_GLYCAN = "[CX4;$([C]-[OX2H+]),CR2:1]-[CX3;$(CN):2](-[NX3H2+:9])=[CX3:3]" #test

RED_GLYCAN_LAC = "[CH2X4;$([C]-[OX2]):1]-[CH1X4;$(CN):2](-[NX3H1:9])-[CH1X4:3]-[O:4]-[CH1X4:5](-[CH3:6])-[CX3:7](=[O:8])-[NX3:94]"
RED_GLY_FRAG_LAC_SMARTS = f"{RED_GLYCAN_LAC}>>{FRAG_LAC_PEPTIDE}.{FRAG_LAC_GLYCAN}" # only red
RED_GLY_FRAG_LAC = make_rxn(RED_GLY_FRAG_LAC_SMARTS)

NONRED_GLYCAN_LAC = "[CX4;$([C]-[OX2H0]),CR2:1]-[CH1X4;$(CN):2](-[NX3H1:9])-[CH1X4:3]-[O:4]-[CH1X4:5](-[CH3:6])-[CX3:7](=[O:8])-[NX3:94]"
NONRED_GLY_FRAG_LAC_SMARTS = f"{NONRED_GLYCAN_LAC}>>{FRAG_LAC_PEPTIDE}.{FRAG_LAC_GLYCAN}"
NONRED_GLY_FRAG_LAC = make_rxn(NONRED_GLY_FRAG_LAC_SMARTS)

# %%%% Muramic Fragmentation

# This makes 138...
MUR_BREAK_0 = "[NX3:1]-[CX4:2]-[CX4:3](-[O:4]-[CX4:5]-[CH3:6])-[CX4:7](-[O:8])-[CX4:9](-[OX2:10])-[CX4:11]"
MUR_BREAK_138 = "[CH]1:[C](-[OH2+]):[CH]:[N]:[CH]:[C]:1-[O]-[CH]=[CH2]"
GLY_FRAG_MUR_BREAK_SMARTS = f"{MUR_BREAK_0}>>{MUR_BREAK_138}.[O]"
GLY_FRAG_MUR_BREAK = make_rxn(GLY_FRAG_MUR_BREAK_SMARTS)

# %%%% Glucosidic Fragmentation 1 and 2

# This makes 168 and 126...
# GLC_BREAK_0 = "CC(CCC(CO)O)[N:26]"
GLC_BREAK_0 = "[OX2:1]-[CX4:2]-[CX4:3](-[OX2:4])-[CX4:5](-[O:6])-[CX4:7](-[OX2:8])-[CX4:9](-[NX3:10])-[CX4:11]-[OX2:12]"

# GLC_BREAK_186 = "[C](=[C])-[C](-[OH2+])-[C](-[OH])-[C](-[N:26])=[C](-[O])"
GLC_BREAK_186 = "[OX2:1]-[CX4:2]-[CX4:3](-[OX2:4])-[CX3:5](-[OH2+])=[CX3:7]-[CX3:9](-[NX3:10])=[CX3:11]"

# GLC_BREAK_168 = "[C](=[C])-[C](-[OH2+])=[C]-[C](-[N:26])=[C](-[O])"
GLC_BREAK_168 = "[CX3:2]=[CX3:3](-[OX2:4])-[CX3:5](-[OH2+])=[CX3:7]-[CX3:9](-[NX3:10])=[CX3:11]"

GLC_BREAK_126 = "[CX3:2]=[CX3:3](-[OX2:4])-[CX3:5](-[OH2+])=[CX3:7]-[CX3:9](-[NH2])=[CX3:11]"
GLY_FRAG_GLC_BREAK1_SMARTS = f"{GLC_BREAK_0}>>{GLC_BREAK_186}.{GLC_BREAK_168}"
GLY_FRAG_GLC_BREAK2_SMARTS = f"{GLC_BREAK_168.replace('OH2+','OH')}>>{GLC_BREAK_126}.O"
GLY_FRAG_GLC_BREAK1 = make_rxn(GLY_FRAG_GLC_BREAK1_SMARTS)
GLY_FRAG_GLC_BREAK2 = make_rxn(GLY_FRAG_GLC_BREAK2_SMARTS)

# %%%% Metal Adducts

# Lactate Metal

LAC_MOIETY = '[CX4;$([C]-[OH]),CR2:1]-[CH1X4;$(CN):2]-[CH1X4:3]-[O:4]-[CH1X4:5](-[CH3:6])-[CX3:7](=[O:8])-[OX2,NX3:9]'

def create_lactate_metal_rxn(metal):
    adduct = f"({FRAG_LAC_PEPTIDE}.[{metal}+])".replace("[OH2+:4]","[OH:4]")
    GLY_FRAG_M_LAC_SMARTS = f'''{LAC_MERGER_GLYCOPEPTIDE}>>{adduct}.O'''
    return make_rxn(GLY_FRAG_M_LAC_SMARTS)


def create_metal_adduct_rxn(metal):
    adduct = f"({LAC_MOIETY}.[{metal}+])"
    GLY_FRAG_M_ADDUCT_SMARTS = f'''{LAC_MOIETY}>>{adduct}.O'''
    return make_rxn(GLY_FRAG_M_ADDUCT_SMARTS)


GLY_FRAG_LAC_Na = create_lactate_metal_rxn("Na")
GLY_FRAG_LAC_K = create_lactate_metal_rxn("K")
GLY_FRAG_ADDUCT_Na = create_metal_adduct_rxn("Na")
GLY_FRAG_ADDUCT_K = create_metal_adduct_rxn("K")

# Glycan BY - omit

# %%%% Frag Types

GLYCAN_FRAG_TYPES = {"Gly. B/Y": NONRED_GLY_FRAG_BY, # non-reducing
                     "Gly. C/Z": NONRED_GLY_FRAG_CZ, # non-reducing
                     "Lac": NONRED_GLY_FRAG_LAC, # non-reducing
                     "Gly. B/Y[r]": RED_GLY_FRAG_BY, # reducing
                     "Gly. C/Z[r]": RED_GLY_FRAG_CZ, # reducing
                     "Lac[r]": RED_GLY_FRAG_LAC, # reducing
                     "Mur 1": GLY_FRAG_MUR_BREAK,
                     "Glc 1": GLY_FRAG_GLC_BREAK1,
                     "Glc 2": GLY_FRAG_GLC_BREAK2,
                     "Lac [Na]": GLY_FRAG_LAC_Na, # not in use
                     "Lac [K]": GLY_FRAG_LAC_K, # not in use
                     "Adduct [Na]": GLY_FRAG_ADDUCT_Na, # not in use
                     "Adduct [K]": GLY_FRAG_ADDUCT_K # not in use
                     }

# %% Uncharge Reactions

# Uncharges by forming a 5mbr ring across amino acid
MISC_RXN_UNCHARGE_CO_A_SMARTS = "[CX2:1](#[O+:2])-[CX4:3]-[NH1X3:4]-[CX3:5]=[O:6]>>[CX3:1]1(=[O+0:2])-[CX4:3]-[NX3:4]=[CX3:5]-[O:6]-1"
MISC_RXN_UNCHARGE_CO_A = make_rxn(MISC_RXN_UNCHARGE_CO_A_SMARTS)

# Uncharges by forming a 6mbr ring across isoGlu/isoGln
MISC_RXN_UNCHARGE_CO_B_SMARTS = "[CX2:1](#[O+:2])-[CX4:3]-[CX4:4]-[CX4:5]-[CX3:6](=[O:7])-[NX3;H2,H1:8]>>[CX3:1]1(=[O+0:2])-[CX4:3]-[CX4:4]-[CX4:5]-[CX3:6](-[O:7]-1)=[NX3:8]"
MISC_RXN_UNCHARGE_CO_B = make_rxn(MISC_RXN_UNCHARGE_CO_B_SMARTS)

# Uncharges by removing metal ion
MISC_RXN_UNCHARGE_M_SMARTS = "([*:1].[Na,K:2])>>[*:1]"
MISC_RXN_UNCHARGE_M = make_rxn(MISC_RXN_UNCHARGE_M_SMARTS)

# %% Outputs

# Rxns
PEPTIDE_RXN_TYPES
GLYCAN_RXN_TYPES

# Molecular Patterns
PEPTIDE_PATT_TYPES
MISC_PATT_TYPES

# Fragmentation Rxns
PEPTIDE_FRAG_TYPES
GLYCAN_FRAG_TYPES
