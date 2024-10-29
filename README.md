# PGN_MS2: _In silico_ MS/MS prediction for peptidoglycan profiling

## Introduction
PGN_MS2 is a computational tool that generates a customizable peptidoglycan (PGN) database from user-defined parameters.
Furthermore, it can simulate MS/MS spectra for each PGN and compile these predicted MS/MS spectra to a spectral library in the NIST format (.msp).
The spectral library (.msp) is compatible with open-access and vendor software, e.g. MS-DIAL, for automated matching and scoring of experimental MS/MS peaks, facilitating automated PGN identification.
Read the open access paper [here.](https://doi.org/10.1039/D3SC05819K) 

![Figure 1](https://github.com/jerickwan/PGN_MS2_private/assets/95602149/b0c08ea9-5efd-430c-8118-95cad216cbab)

![Workflow](https://github.com/jerickwan/PGN_MS2_private/assets/95602149/12a82d7c-4bff-4a2d-a836-efb28f84a7d7)

## Installation
PGN_MS2 is written in Python 3.9 and uses [RDKit](https://www.rdkit.org/) to manipulate molecules. A graphical user interface (built with [easygui](https://github.com/robertlugg/easygui)) is available.
The following Python packages are required:
```
Common
easygui
IPython
joblib
matplotlib
numpy
openpyxml
pandas
pyyaml
rdkit
```
The Python environment can be created with [Conda](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) using the following command:
```
conda env create -f /PATH/OF/PGN_MS2/environment.yaml
```
## Running PGN_MS2
### GUI
The GUI of PGN_MS2 can be run from command line with the following command:
```
python /PATH/OF/PGN_MS2/UserInterface.py
```
A more detailed user guide for the GUI can be found [here.](https://github.com/user-attachments/files/15529917/SuppInfo1_User.Guide.to.PGN_MS2.v2.pdf)
### Running from IDE
Alternatively, PGN_MS2 can be run with an IDE.* Sample code is provided with ManualRun.py

*Spyder 3.9 is not compatible with RDKit and must be ran from a separate environment. See Spyder's [FAQ](https://docs.spyder-ide.org/current/faq.html#using-existing-environment) for more information.
### Output
Output is stored in /output.Each file is named with a prefix comprising the starting datetime and a user-given name (e.g. 20240605_Ecoli).
The various outputs are divided among the three subfolders as such:
| Subfolder     | Filename                                   | Description                                                                                   |
| :------------ |:-------------------------------------------| :---------------------------------------------------------------------------------------------|
| compounds     | [prefix].xlsx                              | MS1 database in spreadsheet format. Monomers, dimers and trimers are shown on separate sheets.|
|               | [prefix].pickle                            | MS1 Database in pickle format (in batches of 5,000 compounds).                                |
|               | [prefix].yaml                              | User-defined settings saved in yaml format.                                                   |
|               | [prefix]_graphical_summary.svg             | MS1 Database in pickle format (in batches of 5,000 compounds).                                |
| msp           | [prefix].msp                               | MS2 database. Different adduct forms are given as separate entries. Fragment ions with same _m/z_ (but different structures) are combined and only the top 50 most intense ions are saved.|
|               | [prefix]_graphical_summary.svg             | MS1 Database in pickle format (in batches of 5,000 compounds).                                |
|               | [prefix]_[number].pickle                   | MS2 Database in pickle format. Saved in batches of 5,000 compounds, which is indicated by [number].|
| peaklists     | [prefix]_spectradata.xlsx                  | MS2 database in spreadsheet format. Each batch has its own sheet. Each compound is presented as its own table containing the top 200 most intense ions.|
|               | [prefix]_iondata.xlsx                      | All ions and their respective structures are tabulated in this fil                            |
### Supported PGN Chemotypes
PGN_MS2 imports chemical information from an internal library located at:
```
data/PGN.xlsx
```
PGN_MS2 was designed to accomodate most [PGN chemotypes.](https://doi.org/10.1111/j.1574-6976.2007.00094.x) It is able to generate PGN with:
* modified glycans: acetylation (increase/decrease), glycolylation (anMurNGlyc) and dehydration (anMurNAc).
* stem peptide sequences up to eight amino acids long. Supported amino acids include the canonical amino acids as well as non-canonical amino acids commonly found in PGN (mDAP, Orn, γ-isoGln).
* bridge peptides (i.e. branch peptides, side chains) that are attached to either diamino/dicarboxy amino acids in the stem peptide.
* a wide variety of modifications such as lactamization, endopeptidase digestion.
* two different polymerisation modes: either through glycosidic bonds or peptide bonds.
## Misc / Other Links
This [tool](https://www.yqiaolab.com/pgn_ms2-tool) was built by members of [Qiao Lab.](https://www.yqiaolab.com)
MS/MS spectra for all identified PGN from the a/m paper is also available as a download on [MoNA.](https://mona.fiehnlab.ucdavis.edu/downloads)
PGN_MS2 was used in combination with MS-DIAL, an open source MS analysis software, available [here.](https://systemsomicslab.github.io/compms/msdial/main.html)

The following can be found in the [Supplementary Information](https://www.rsc.org/suppdata/d3/sc/d3sc05819k/d3sc05819k1.pdf) of our paper:
* Nomenclature (Table S1)
* Overview of GUI (Table S2)

Read the open access paper [here.](https://doi.org/10.1039/D3SC05819K) 
