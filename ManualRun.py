'''
Testing
env = chem
'''

from MSPMaker import MSPMaker
from Molecules import Molecule, AminoAcid, Peptide, Glycan, Peptidoglycan
from Generator import Generator

# %% Globals

CURRENT_NAME = "Test"  # No spacing

# %% Generator

adducts = ['[M+H]+', '[M+Na]+', '[M+2Na]2+', '[M+K]+',
           '[M+2H]2+', '[M+3H]3+']  # adducts
gen = Generator(CURRENT_NAME,
                output_adducts=adducts)

# %%% Modifications

# all available modifications shown
gen.modifications["Alanine/Lactate Substitution"] = False
gen.modifications["Muramic Lactam"] = False
gen.modifications["Muramic Reduction"] = False
gen.modifications["Lactoyl Peptide"] = False
gen.modifications["Braun LPP"] = False
gen.modifications["EPase P1"] = False
gen.modifications["EPase P2"] = False

# %%% Lengths

gen.set_length("peptide", 0, 5)
gen.set_length("glycan", 0, 2)

# %%% Glycan

# see PGN.xlsx for valid glycan units
gen.set_glycan_units("Glc",
                     ["GlcNAc",
                      "GlcNAcOAc",
                      "GlcN"
                      ])
gen.set_glycan_units("Mur",
                     ["MurNAc",
                      "MurN",
                      "MurNAcOAc",
                      "anMurN",
                      "anMurNAc"])
# %%% Peptide

# see PGN.xlsx for valid peptide units
AA_list_1 = ["Ala", "Gly", "Lys", "Phe",
             "Cys", "Leu", "Pro", "Met",
             "Trp", "Tyr", "His", "Arg",
             "Ser", "Asp", "Asn"]
AA_list_2 = ["Ala", "Gly", "Lys"]

gen.set_peptide_residues(1, ["Ala"])
gen.set_peptide_residues(2, ["γ-Glu", "γ-isoGln"])
# gen.set_peptide_residues(3, ["mDAP", "mDAP(NH2)"])
gen.set_peptide_residues(3, ["Lys", "Orn", "Lan", "Lan(NH2)",
                             "mDAP", "mDAP(NH2)"])
# gen.set_peptide_residues(3, ["Lys", "Orn"])
gen.set_peptide_residues(4, AA_list_2)
gen.set_peptide_residues(5, AA_list_2)

# %%% Crosslinks

gly_allowed = [
    # "GlcN",
    "GlcNAc",
    # "GlcNAcOAc",
    # "MurN",
    "MurNAc",
    # "MurNAcOAc",
    "anMurNAc"]

gen.set_num_polymerisations(1)
gen.set_diffcalc_units(1, [1, 0])
gen.set_diffcalc_param(gly_len=[2],
                       gly_range=[None,
                                  "GlcNAc",
                                  "MurNAc"
                                  ],
                       pep_len=[4, 5],
                       pep_range=["Ala",
                                  "γ-Glu", "γ-isoGln",
                                  "mDAP", "mDAP(NH2)",
                                  "Lys", "Orn",
                                  None],
                       bridge_range=None,
                       polymerisation_range=None)

# lateral chains from innermost AA to outermost AA
chains = [
    ["Ser", "Ala", "Thr", "Ala"],
    ["Gly", "Gly", "Gly", "Gly", "Gly"],
    ["β-isoAsn",],
    ["β-Asp",],
    ["Ala",],
    ["Ala", "Ala"],
    ["Ala", "Ser"],
    ["Ala", "Ala", "Ala"],
    ["Ala", "Ala", "Ser"]
]

gen.set_bridge_peptides(3, "NH2",
                        chains,
                        valid_AAs=["Lys", "Orn"])

# %%%% 'G>G'

# gen.set_polymerisation_types(
#     P1={2: {"allowed": ["any"]},
#         3: {"allowed": ["mDAP", "Lys", None]},
#         4: {"allowed": ["Ala", None]},
#         5: {"allowed": ["Ala", None]},
#         "pep_len": [0, 3, 4, 5],
#         "gly_len": [2],
#         #"gly_allowed": gly_allowed,
#         "gly_rejected": ["anMurNAc", "anMurN", "anMurNGlyc"]},
#     P2={2: {"allowed": ["any"]},
#         3: {"allowed": ["mDAP", "Lys", None]},
#         4: {"allowed": ["Ala", None]},
#         5: {"allowed": ["Ala", None]},
#         "pep_len": [0, 3, 4, 5],
#         "gly_len": [2],
#         #"gly_allowed": gly_allowed
#         },
#     P1_bond=[0, "glycan", "glycan", "Acc"],
#     P2_bond=[0, "glycan", "glycan", "Dnr"],
#     kmer_range=[1, 2])

# %%%% '3s-4'

# Lys
gen.set_polymerisation_types(
    P1={2: {"allowed": ["any"]},
        3: {"allowed": ["Lys"],
            "bridge": False},
        4: {"allowed": ["Ala", "Gly", None]},
        5: {"allowed": ["Ala", "Gly", None]},
        "pep_len": [3, 4, 5],
        "gly_len": [0, 2],
        "gly_allowed": gly_allowed},
    P2={2: {"allowed": ["any"]},
        3: {"allowed": ["Lys"],
            "bridge": False},
        4: {"allowed": ["Ala", "Gly"]},
        "pep_len": [4],
        "gly_len": [0, 2],
        "gly_allowed": gly_allowed},
    P1_bond=[3, "main", "side", "NH2"],
    P2_bond=[4, "main", "eupeptide", "COOH"])

# mDAP
gen.set_polymerisation_types(
    P1={2: {"allowed": ["any"]},
        3: {"allowed": ["mDAP"],
            "bridge": False},
        4: {"allowed": ["Ala", "Gly", None]},
        5: {"allowed": ["Ala", "Gly", None]},
        "pep_len": [3, 4, 5],
        "gly_len": [0, 2],
        "gly_allowed": gly_allowed},
    P2={2: {"allowed": ["any"]},
        3: {"allowed": ["mDAP"],
            "bridge": False},
        4: {"allowed": ["Ala", "Gly"]},
        "pep_len": [4],
        "gly_len": [0, 2],
        "gly_allowed": gly_allowed},
    P1_bond=[3, "main", "side", "NH2"],
    P2_bond=[4, "main", "eupeptide", "COOH"])

# mDAP(NH2)
gen.set_polymerisation_types(
    P1={2: {"allowed": ["any"]},
        3: {"allowed": ["mDAP(NH2)"],
            "bridge": False},
        4: {"allowed": ["Ala", "Gly", None]},
        5: {"allowed": ["Ala", "Gly", None]},
        "pep_len": [3, 4, 5],
        "gly_len": [0, 2],
        "gly_allowed": gly_allowed},
    P2={2: {"allowed": ["any"]},
        3: {"allowed": ["mDAP(NH2)"],
            "bridge": False},
        4: {"allowed": ["Ala", "Gly"]},
        "pep_len": [4],
        "gly_len": [0, 2],
        "gly_allowed": gly_allowed},
    P1_bond=[3, "main", "side", "NH2"],
    P2_bond=[4, "main", "eupeptide", "COOH"])

# %%%% '3s-3'
# gen.set_polymerisation_types(
#     P1={2: {"allowed": ["any"]},
#         3: {"allowed": ["Lan"],
#             "bridge": False},
#         4: {"allowed": ["Ala", "Gly", "Lys", None]},
#         5: {"allowed": ["Ala", None]},
#         "pep_len": [3, 4, 5],
#         "gly_len": [0, 2],
#         "gly_allowed": gly_allowed},
#     P2={2: {"allowed": ["any"]},
#         3: {"allowed": ["Lan"],
#             "bridge": False},
#         "pep_len": [3],
#         "gly_len": [0, 2],
#         "gly_allowed": gly_allowed},
#     P1_bond=[3, "main", "side", "NH2"],
#     P2_bond=[3, "main", "eupeptide", "COOH"])

# %%%% '3br-4'
# Ala-Ala
gen.set_polymerisation_types(
    P1={2: {"allowed": ["any"]},
        3: {"allowed": ["Lys"],
            "bridges_allowed": [
                ["Ala", "Ala"]]},
        4: {"allowed": ["Ala", None]},
        5: {"allowed": ["Ala", None]},
        "pep_len": [3, 4, 5],
        "gly_len": [2],
        "gly_allowed": gly_allowed},
    P2={2: {"allowed": ["any"]},
        3: {"allowed": ["Lys"],
            "bridges_allowed": [
                None,
                ["Ala", "Ala"]]},
        4: {"allowed": ["Ala", None]},
        "pep_len": [4],
        "gly_len": [0, 2],
        "gly_allowed": gly_allowed},
    P1_bond=[3, "bridge", "eupeptide", "NH2"],
    P2_bond=[4, "main", "eupeptide", "COOH"],
    kmer_range=[0, 99])

# Ala-Ala-Ala
# gen.set_polymerisation_types(
#     P1={2: {"allowed": ["any"]},
#         3: {"allowed": ["Lys"],
#             "bridges_allowed":[
#                 ["Ala","Ala","Ala"]]},
#         4: {"allowed": ["Ala", None]},
#         5: {"allowed": ["Ala", None]},
#         "pep_len": [3, 4, 5],
#         "gly_len": [2],
#         "gly_allowed": gly_allowed},
#     P2={2: {"allowed": ["any"]},
#         3: {"allowed": ["Lys"],
#             "bridges_allowed":[
#                 None,
#                 ["Ala","Ala","Ala"]]},
#         4: {"allowed": ["Ala", None]},
#         "pep_len": [4],
#         "gly_len": [0,2],
#         "gly_allowed": gly_allowed},
#     P1_bond=[3, "bridge", "eupeptide", "NH2"],
#     P2_bond=[4, "main", "eupeptide", "COOH"],
#     kmer_range=[0, 99])

# Ala
gen.set_polymerisation_types(
    P1={2: {"allowed": ["any"]},
        3: {"allowed": ["Lys"],
            "bridges_allowed": [
                ["Ala"]]},
        4: {"allowed": ["Ala", None]},
        5: {"allowed": ["Ala", None]},
        "pep_len": [3, 4, 5],
        "gly_len": [2],
        "gly_allowed": gly_allowed},
    P2={2: {"allowed": ["any"]},
        3: {"allowed": ["Lys"],
            "bridges_allowed": [
                None,
                ["Ala"]]},
        4: {"allowed": ["Ala", None]},
        "pep_len": [4],
        "gly_len": [0, 2],
        "gly_allowed": gly_allowed},
    P1_bond=[3, "bridge", "eupeptide", "NH2"],
    P2_bond=[4, "main", "eupeptide", "COOH"],
    kmer_range=[0, 99])

# Asx
gen.set_polymerisation_types(
    P1={2: {"allowed": ["any"]},
        3: {"allowed": ["Lys"],
            "bridges_allowed": [
                ["β-isoAsn"],
                ["β-Asp"]]},
        4: {"allowed": ["Ala", None]},
        5: {"allowed": ["Ala", None]},
        "pep_len": [3, 4, 5],
        "gly_len": [2],
        "gly_allowed": gly_allowed},
    P2={2: {"allowed": ["any"]},
        3: {"allowed": ["Lys"],
            "bridges_allowed": [
                None,
                "β-isoAsn",
                "β-Asp"]},
        4: {"allowed": ["Ala", None]},
        "pep_len": [4],
        "gly_len": [0, 2],
        "gly_allowed": gly_allowed},
    P1_bond=[3, "bridge", "eupeptide", "NH2"],
    P2_bond=[4, "main", "eupeptide", "COOH"],
    kmer_range=[0, 99])

# %%%% '3br-3'
# gen.set_polymerisation_types(
#     P1={2: {"allowed": ["any"]},
#         3: {"allowed": ["Lys"],
#             "bridge": True},
#         4: {"allowed": ["Ala", "Lys", "Gly", None]},
#         5: {"allowed": ["Ala", "Lys", "Gly", None]},
#         "pep_len": [3, 4, 5],
#         "gly_len": [0, 2]},
#     P2={2: {"allowed": ["any"]},
#         3: {"allowed": ["Lys"]},
#         "pep_len": [3],
#         "gly_len": [0, 2]},
#     P1_bond=[3, "bridge", "eupeptide", "NH2"],
#     P2_bond=[3, "main", "eupeptide", "COOH"])

# %%% Execution
gen.generate_PGN()
gen.export_settings_as_yaml()
gen.export_settings_as_image()
pkl = gen.export_PGN_as_pickle(kmers=[1, 2])
gen.export_dataframes(["xlsx"])

# %% MSPMaker

maker = MSPMaker(output_adducts=['[M+H]+', '[M+2H]2+', '[M+3H]3+'])
maker.import_PGN_as_pickle(pkl, export="all")
