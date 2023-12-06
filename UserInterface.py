'''
User Interface, Testing
env = chem
'''
# %% Imports
import easygui
import yaml
import sys
import os
from pathlib import Path

# %%% Resolve Paths
frozen = getattr(sys, 'frozen', False)
if frozen:
    wd = Path(sys.executable).parent
    os.chdir(wd)
    print('os.getcwd is', os.getcwd() )

# %%% Import Customs

from MSPMaker import MSPMaker, Fragmenter
from Molecules import Molecule, AminoAcid, Peptide, Glycan, Peptidoglycan
from Generator import Generator
from Exceptions import InputError
from Common import TIME_STRING, MEMORY

# %% GUI


def get_preselect_integers(full_range, preselect_range):
    """
    Helper function for easygui.multchoicebox(). Returns indexes of all strings
    from preselect_range in full_range (if found).

    Parameters
    ----------
    full_range : list
        List of options (as strings)
    preselect_range : list
        List of pre-selected options (as strings)

    Returns
    -------
    preselect : list
        List of indexes. If there is no overlap, return 0.

    """
    if preselect_range == 0:
        return 0
    preselect = []
    for i, obj in enumerate(full_range):
        if obj in preselect_range:
            preselect.append(i)
    if len(preselect) == 0:
        return 0
    else:
        return preselect


class GUI():
    PARAM_FOLDER = Path("settings/ui")
    PARAM_FILEPATH = Path("settings/ui/parameters.yaml")
    HISTORY_MAX = 10

    def __init__(self):
        self.init_title = "PGN_MS2"
        self.title = self.init_title
        self._history = [self.init_title]
        self.generator = None
        self.maker = None
        self.settings = {}
        self.load_GUI_settings()
        self.name = None
        self.retries = 0
        self.start_step = 0

    def rename(self):
        """
        Renames the Generator and the GUI window.

        Returns
        -------
        None.

        """
        name = easygui.enterbox("Enter a name.",
                                self.title)
        self.name = name
        self.history = f"Name: {self.name}"
        if self.generator is not None:
            self.generator.refresh_name(self.name)
        else:
            self.create_generator()

    @property
    def history(self):
        return "\n".join(self._history)+"\n\n"

    @history.setter
    def history(self, text):
        if text not in (None, ""):
            self._history.append(text)
        if len(self._history) > GUI.HISTORY_MAX:
            self._history.pop(0)

    def clear_history(self):
        """
        Clears self.history.

        Returns
        -------
        None.

        """
        self._history = [self.title]

    def start(self):
        """
        Starts the GUI.

        Returns
        -------
        None.

        """
        if self.generator is None:
            self.create_generator()
        num_errors = 0
        prog_run = True
        while num_errors < 3 and prog_run:
            try:
                prog_run = self.run()
            except Exception as exc:
                print(exc)
                easygui.exceptionbox()
                num_errors += 1
        sys.exit()

    def run(self):
        """
        Main loop of GUI.

        Returns
        -------
        boolean
            Returns False to kill loop.

        """
        choice = 0
        while choice not in ["Exit", None]:
            msg = '''Select an option.'''
            choices = [
                "Build", "Load",
                "Edit", "Clear",
                "Run", "More Options",
                "Exit"]
            choice = easygui.choicebox(
                self.history+msg, self.title, choices=choices)
            if choice == "Build":
                self.gen_build_settings()
            elif choice == "Load":
                self.gen_load_settings()
            elif choice == "Edit":
                self.gen_edit_settings()
            elif choice == "Clear":
                self.create_generator(force=True)
            elif choice == "Run":
                return self.gen_execute()
            elif choice == "More Options":
                self.run_extra_options()
        return False

    def run_extra_options(self):
        """
        Runs additional, less used options.

        Raises
        ------
        Exception
            Raises exception (for testing only)..

        Returns
        -------
        boolean
            Returns False to kill loop.

        """
        choice = 0
        while choice not in ["Exit", None]:
            msg = '''Select an option.'''
            choices = [
                "Fragment Single Compound",
                "GUI Settings", "Raise Exception",
                "Reduce Cache", "Delete Cache",
                "Exit"]
            choice = easygui.choicebox(
                self.history+msg, self.title, choices=choices)
            if choice == "Fragment Single Compound":
                self.frag_single_compound()
            elif choice == "GUI Settings":
                self.menu_GUI_settings()
            elif choice == "Raise Exception":
                raise Exception("Testing!")
            elif choice == "Reduce Cache":
                print("Reducing cache... please wait!")
                MEMORY.reduce_size()
            elif choice == "Delete Cache":
                print("Deleting cache... please wait!")
                MEMORY.clear(warn=False)
            else:
                pass
        return False

    # %%% GUI Settings

    def menu_GUI_settings(self):
        """
        Submenu for editting GUI settings.

        Returns
        -------
        None.

        """
        msg = "Select an option."
        choice = 0
        choices = [
            "Load",
            "Load Default",
            "Save",
            "Return"]
        while choice not in ["Return", None]:
            choice = easygui.choicebox(
                self.history+msg, self.title, choices=choices)
            if choice == "Load":
                self.load_GUI_settings(use_default=False)
            elif choice == "Load Default":
                self.load_GUI_settings(use_default=True)
            elif choice == "Save":
                self.save_GUI_settings()
            else:
                pass

    def load_GUI_settings(self, use_default=True):
        """
        Loads settings for the GUI.

        Parameters
        ----------
        use_default : boolean, optional
            If True, uses the hardcoded settings; else loads settings from
            a file The default is True.

        Returns
        -------
        None.

        """
        if use_default:
            self.load_default_GUI_settings()
        else:
            msg = "Open a settings file (.yaml)"
            file = easygui.fileopenbox(
                msg,
                default=GUI.PARAM_FOLDER/"*.yaml",
                filetypes=["*.yaml"])
            filepath = Path(file)
            with open(filepath, "r") as f:
                self.settings = yaml.safe_load(f)
            self.history = f"Loaded settings from: {filepath}"

    def load_default_GUI_settings(self):
        """
        Loads the default hardcoded settings.

        Returns
        -------
        None.

        """
        self.settings = {}
        self.settings["default_Glc_types"] = "GlcNAc,GlcN,GlcNAcOAc"
        self.settings["default_Mur_types"] = "MurNAc,MurN,MurNAcOAc,anMurNAc"
        self.settings["default_tripeptide"] = {1: "Ala",
                                               2: "γ Glu,γ isoGln",
                                               3: "mDAP,mDAP(NH2),Lys,Orn"}
        self.settings["full_AA_range"] = \
            "Ala,Gly,Lys,Phe,Cys,Leu,Pro,Met,Trp,Tyr,His,Arg,Ser,Asp,Asn"
        self.settings["limited_AA_range"] = '''Ala,Gly,Lys'''
        self.history = "Loaded default settings"

    def save_GUI_settings(self, overwrite_default=False):
        """
        Saves the GUI settings.

        Parameters
        ----------
        overwrite_default : boolean, optional
            If True, overwrites, the default parameter file at GUI.PARAM_FILEPATH.
            The default is False.

        Returns
        -------
        None.

        """
        filepath = None
        if overwrite_default:
            filepath = GUI.PARAM_FILEPATH
        else:
            if self.name is None:
                self.rename()
            filepath = Path(f"settings/ui/{self.name}_parameters.yaml")

        with open(filepath, "w") as f:
            yaml.safe_dump(self.settings, f)
        self.history = f"Saved settings as {filepath}"

    # %%% Generator

    def create_generator(self, force=False):
        """
        Creates a Generator with same name as GUI.

        Parameters
        ----------
        force : boolean, optional
            If False, creates a Generator only if it doesn't exist or if its name
            doesn't match. The default is False.

        Returns
        -------
        bool
            Returns True if a new Generator is created.

        """
        if self.generator is None or self.generator.name != self.name or force:
            self.generator = Generator(self.name)
            self.history = "Generator created."
            return True
        else:
            return False

    def gen_load_settings(self):
        """
        Loads settings for the Generator from a .yaml file.

        Returns
        -------
        boolean
            Returns True if settings loaded successfully.

        """
        try:
            default = "*.yaml"
            settings_file = Path(easygui.fileopenbox(default=default))
            if settings_file.suffix == ".yaml":
                self.generator.import_settings_as_yaml(settings_file)
            elif settings_file.suffix == ".pickle":
                self.generator.import_settings_as_pickle(settings_file)
            else:
                self.history = "Not a valid settings file."
                return False
        except TypeError:
            return False

        self.history = f"Generator settings loaded from \n{settings_file}"
        return True

    def gen_build_settings(self, exit_step=12,
                           max_retries=3):
        """
        Builds settings for the Generator.
        Settings are generated in 12 steps:
           1: Peptide, Glycan Lengths
           2: Glc-Type Glycans
           3: Mur-Type Glycans
           4-8: Peptide Residues for positions 1-5
           9: Bridge Peptide
           10: Modifications
           11: Polymerisation
           12: Difference Calculator

        Parameters
        ----------
        exit_step : int, optional
            Indicates the integer of the last step.
            Use 12 to exit after Difference Calculator. The default is 12.
        max_retries : int, optional
            Maximum no. of times the user is allowed to backtrack.

        Returns
        -------
        None.

         """
        if self.start_step == 0:
            self.start_step += 1
            self.retries = 0
        defaults = [
            "Ala",
            "γ-Glu,γ-isoGln",
            "mDAP,mDAP(NH2),Lys,Orn",
            self.settings["full_AA_range"],
            self.settings["full_AA_range"]]
        # Loop
        while self.retries <= max_retries and exit_step >= self.start_step:
            # Length
            if self.start_step == 1:
                self.rename()
                self.gen_set_length(mod_values=[0, 2, 0, 5, 1])
            # Glycans - Glc-type
            if self.start_step == 2:
                self.gen_set_glycan("Glc",
                                    default=self.settings["default_Glc_types"])
            # Glycans - Mur-type
            if self.start_step == 3:
                self.gen_set_glycan("Mur",
                                    default=self.settings["default_Mur_types"])
                if self.start_step == 4:
                    # Check for reduction
                    reduced_muro = False
                    for glycan in self.generator.glycan_units["Mur"]:
                        if "[r]" in glycan:
                            reduced_muro = True
                            print("Reduced Mur-type sugars preferred\
                                  in subsequent settings.")
            # Peptide
            if 4 <= self.start_step < 9:
                idx = self.start_step-3
                self.gen_set_peptide(idx, default=defaults[idx-1])
            # Bridge Peptides
            if self.start_step == 9:
                self.gen_set_bridge_peptides()
            # Modifications
            if self.start_step == 10:
                self.gen_set_modifications()
            # Crosslinks
            if self.start_step == 11:
                self.gen_choose_default_polymerisations(
                    prefer_reduced=reduced_muro)
            # DiffCalc
            if self.start_step == 12:
                self.history = self.generator.set_diffcalc_units(1, [1])
                if reduced_muro:
                    gly_range = [None, "GlcNAc", "MurNAc[r]", "anMurNAc"]
                else:
                    gly_range = [None, "GlcNAc", "MurNAc", "anMurNAc"]
                if self.generator.modifications["EPase P1"] or \
                        self.generator.modifications["EPase P2"]:
                    gly_len = [0, 2]
                else:
                    gly_len = [2]
                self.gen_set_diffcalc_params(gly_len=gly_len,
                                             gly_range=gly_range,
                                             pep_len=[4, 5])
        # Terminate
        self.start_step = 0
        self.retries = 0

    def gen_edit_settings(self):
        """
        Opens a menu that allows the user to edit the settings for Generator
        after the initial sequence.

        Returns
        -------
        None.

        """
        msg = "Select an option."
        choice = 0
        choices = ["Modifications", "Length",
                   "Glycan Units", "Peptide Residues",
                   "Bridge Peptides", "Polymerisations",
                   "Diffcalc Units", "Diffcalc Params",
                   "Return"]
        while choice not in ("Return", None):
            choice = easygui.choicebox(
                self.history+msg, self.title, choices=choices)
            if choice == "Modifications":
                self.gen_set_modifications()
            elif choice == "Length":
                self.gen_set_length()
            elif choice == "Glycan Units":
                self.gen_set_cpds("glycan")
            elif choice == "Peptide Residues":
                self.gen_set_cpds("peptide")
            # Bridge Peptides
            elif choice == "Bridge Peptides":
                self.gen_set_bridge_peptides()
            # Crosslinks
            elif choice == "Polymerisations":
                self.gen_choose_default_polymerisations()
            elif choice == "Diffcalc Units":
                self.gen_set_diffcalc_units()
            elif choice == "Diffcalc Params":
                self.gen_set_diffcalc_params()
            else:
                pass

    def register_step(self, success=True):
        """
        Helper function to help keep track of gen_build_settings.

        Parameters
        ----------
        success : boolean, optional
            If True, current step ran successfully. The default is True.

        Returns
        -------
        None.

        """
        if success:
            self.start_step += 1
            print(f"Advancing to step [{self.start_step}]")
        else:
            self.retries += 1
            if self.start_step > 1:
                self.start_step -= 1
            trace_msg = f"--- Returning to previous step [{self.start_step}] ---"
            self.history = trace_msg

    # %%%% Modifications

    def gen_set_modifications(self):
        """
        Selects modifications for PGN Generator.

        Returns
        -------
        None.

        """
        msg = "Select modifications."
        choices = list(self.generator.modifications.keys())
        selected = [i for i, mod in enumerate(
            choices) if self.generator.modifications[mod]]
        choices.append("No modifications")
        picked = easygui.multchoicebox(
            self.history+msg, self.title, choices=choices, preselect=selected)
        if picked is None:
            self.register_step(False)
        elif "No modifications" in picked:
            picked = []
            self.history = "No modifications added."
            self.register_step(True)
        else:
            self.history = self.generator.add_modifications(picked)
            self.register_step(True)

    # %%%% Length

    def gen_set_length(self, mod_values=None):
        """
        Selects glycan length, peptide length and number of polymerisations.

        Parameters
        ----------
        mod_values : list
            User-modified values for glycan length, peptide length and number
            of crosslinks.

        Returns
        -------
        None.

        """
        if mod_values is None:
            msg = "Set parameters."
            fields = ["Glycan min", "Glycan max",
                      "Peptide min", "Peptide max",
                      "Num. Polymerisations"]
            field_values = [self.generator.glycan_lmin,
                            self.generator.glycan_lmax,
                            self.generator.peptide_lmin,
                            self.generator.peptide_lmax,
                            self.generator.num_polymerisations]
            mod_values = easygui.multenterbox(self.history+msg, self.title,
                                              fields=fields, values=field_values)
        if mod_values is None:
            self.register_step(False)  # Cancel
        else:
            self.history = self.generator.set_length(
                "glycan", *map(int, mod_values[0:2]))
            self.history = self.generator.set_length(
                "peptide", *map(int, mod_values[2:4]))
            self.history = self.generator.set_num_polymerisations(
                int(mod_values[4]))
            self.register_step(True)

    # %%%% Glycan / Peptide

    def gen_set_cpds(self, cpd):
        """
        Initial stage of selection for glycan/peptides.

        Parameters
        ----------
        cpd : string
            Type of compound. "peptide" or "glycan" accepted.

        Raises
        ------
        InputError
            Raised if cpd is not valid.

        Returns
        -------
        None.

        """
        if cpd == "glycan":
            msg = "Select a glycan type"
            gtype = easygui.choicebox(
                self.history+msg,
                self.title,
                choices=["Glc", "Mur"])
            self.gen_set_glycan(gtype)

        elif cpd == "peptide":
            msg = "Select position in stem peptide"
            idx = int(easygui.choicebox(
                self.history+msg,
                self.title,
                choices=range(1, 6)))
            self.gen_set_peptide(idx)

        else:
            raise InputError("Glycan/Peptide",
                             f"{cpd} not recognised")

    def gen_set_glycan(self, gtype, default=None):
        """
        Next stage of selection for glycans.

        Parameters
        ----------
        gtype : string
            Type of glycan. "Glc" or "Mur" accepted.
        default : list, optional
            Default list of compounds. The default is None.

        Raises
        ------
        InputError
            Raised if gtype is not valid.

        Returns
        -------
        None.

        """
        # type mode
        if gtype not in ["Glc", "Mur"]:
            raise InputError("Glycan",
                             f"{gtype} not recognised.")
        msg = f'''
        Type {gtype}-type glycan units, separated by a ','.
        e.g. {gtype}NAc, {gtype}N, {gtype}NAcOAc
        Reduced 'Mur'-type glycans are denoted with [r]; e.g. MurNAc[r]
        '''
        selected = easygui.enterbox(
            self.history+msg, self.title, default=default)
        if selected is None:
            self.register_step(False)
        else:
            selected = selected.split(",")
            selected = set(map(str.strip, selected))
            self.history = self.generator.set_glycan_units(gtype, selected)
            self.register_step(True)

    def gen_set_peptide(self, idx, default=None):
        """
        Next stage of selection for peptides.

        Parameters
        ----------
        idx : int
            Position of amino acid(s) to be added.
        default : list, optional
            Default list of compounds. The default is None.

        Returns
        -------
        None.

        """
        # type mode
        if idx < 1 or idx > 20:
            raise InputError("Peptide",
                             "{idx} outside of range 1-20.")
        msg = f'''
        Type AAs for position {idx}, separated by a ','.
        e.g. Ala, Ser
        '''
        selected = easygui.enterbox(
            self.history+msg, self.title, default=default)
        if selected is None:
            self.register_step(False)
        else:
            selected = selected.split(",")
            selected = set(map(str.strip, selected))
            condition = None
            if idx == 5:
                conditional_bool = easygui.boolbox(
                    msg="Use conditional addition?")
                if conditional_bool:
                    msg = f'''
                    Add AAs at {idx} only if preceding AA is from ...
                    '''
                    condition = easygui.enterbox(
                        self.history+msg, self.title).strip().split(",")
                    condition = set(map(str.strip, condition))
                if condition == [""]:
                    condition = None

            self.history = self.generator.set_peptide_residues(
                idx, selected, precAA_lst=condition)
            self.register_step(True)

    def gen_set_bridge_peptides(self):
        """
        Selects bridge peptides.

        Returns
        -------
        None.

        """
        msg = "Select position of bridge peptide or 0 to skip."
        idx = easygui.choicebox(
            self.history+msg,
            self.title,
            choices=range(0, 6))
        if idx is None:
            self.register_step(False)
            return
        idx = int(idx)
        if idx == 0:
            self.register_step(True)
            return

        msg = "Select connection type"
        grp = easygui.choicebox(
            self.history+msg,
            self.title,
            choices=["COOH", "NH2"])
        # type mode
        msg = f'''
        Type AAs which have a chain at position {idx}-{grp}, separated by a ','.
        'any'\t--> any AA can be used.
        'diamino'\t--> any diamino AA can be used (e.g. Lys, mDAP).
        'dicarboxy'\t--> any dicarboxy AA can be used (e.g. Glu).
        '''
        if grp == "COOH":
            default = "dicarboxy"
        elif grp == "NH2":
            default = "diamino"
        else:
            self.register_step(False)
            return

        valid = easygui.enterbox(
            self.history+msg, self.title, default=default)
        if valid is None:
            self.register_step(False)
        else:
            valid = valid.split(",")
            valid = set(map(str.strip, valid))

        msg = f'''
        Type chains that connect to {idx}-{valid}-{grp}, separated by a ','.
        Amino acids are typed from nearest to stem peptide to furthest
        and separated by a '>'
        e.g. β-Asp, Ala>Ala, Ser>Ala>Thr>Ala
        '''
        chains = easygui.enterbox(
            self.history+msg, self.title)
        if chains is None:
            self.register_step(False)
            return

        chains = chains.split(",")
        chain_lst = [chain.split(">") for chain in chains]

        self.history = self.generator.set_bridge_peptides(
            idx, grp, chain_lst, valid_AAs=valid)
        self.register_step(True)

    def gen_retrieve_AAs(self, idx="all"):
        """
        Retrieve AAs from Generator for positions {idx} and returns a combined,
        sorted list.

        Parameters
        ----------
        idx : list, optional
            List of positions to retrieve. The default is "all", which retrieves
            for positions 1 - 6.

        Returns
        -------
        AAs : list
            Sorted list of amino acids.

        """
        AAs = []
        if idx == "all":
            idx = [1, 6]
        for i in range(*idx):
            if i in self.generator.peptide_residues:
                AAs.extend(self.generator.peptide_residues[i]["AAs"])
            else:
                AAs = Generator.valid_AAs
        AAs = sorted(list(set(AAs)))
        return AAs

    # %%%% Polymerisations

    def gen_choose_default_polymerisations(self, prefer_reduced=False):
        """
        Select polymerisations with default options.

        Parameters
        ----------
        prefer_reduced : boolean, optional
            If True, polymerisations will prefer the reduced Mur-type sugars
            over non-reduced sugars. The default is False.

        Returns
        -------
        None.

        """
        if self.generator.num_polymerisations == 0:
            self.history = "No polymerisations picked."
            self.register_step(True)
        msg = "Pick polymerisations."
        crosslink_options = ["G-G",
                             "3s-3", "3s-4",
                             "3br-3", "3br-4",
                             "Clear existing polymerisations",
                             "No polymerisations"]
        polymerisations = easygui.multchoicebox(
            self.history+msg, self.title, choices=crosslink_options)

        preselect = 0
        if polymerisations is None:
            self.register_step(False)
            return
        else:
            if "No polymerisations" in polymerisations:
                self.history = "No polymerisations picked."
                self.register_step(True)
            if "Clear existing polymerisations" in polymerisations:
                self.history = self.generator.clear_polymerisations()
            if "G-G" in polymerisations:
                preselect = \
                    self.gen_add_default_glycan_polymerisation(
                        "G-G", prefer_reduced, preselect)
            if "3s-3" in polymerisations:
                preselect = \
                    self.gen_add_default_peptide_polymerisation(
                        "3s-3", 3, False, prefer_reduced, preselect)
            if "3s-4" in polymerisations:
                preselect = \
                    self.gen_add_default_peptide_polymerisation(
                        "3s-4", 4, False, prefer_reduced, preselect)
            if "3br-3" in polymerisations:
                preselect = \
                    self.gen_add_default_peptide_polymerisation(
                        "3br-3", 3, True, prefer_reduced, preselect)
            if "3br-4" in polymerisations:
                preselect = \
                    self.gen_add_default_peptide_polymerisation(
                        "3br-4", 4, True, prefer_reduced, preselect)

            self.register_step(True)

    def gen_add_default_glycan_polymerisation(self, name,
                                              prefer_reduced,
                                              preselect=0):
        """
        Select default glycosidic polymerisation (G-G) with preset options.

        Parameters
        ----------
        name : string
            Name of polymerisation.
        prefer_reduced : boolean, optional
            If True, polymerisations will prefer the reduced Mur-type sugars
            over non-reduced sugars. The default is False.
        preselect : list, optional
            List of integers that correspond to AAs. These will appear as the
            default option. The default is 0 (No amino acids).

        Returns
        -------
        preselect : list
            List of integers that correspond to AAs. These will appear as the
            default option in future selections.

        """
        msg = f"Set acceptable AAs in position 4 and 5 in the stem peptide\
            of {name} polymerisation."
        AAs = self.gen_retrieve_AAs(idx=[4, 5])
        if self.generator.modifications["Alanine/Lactate Substitution"]:
            AAs.append("Lac")  # add Lactate if modification used
        AA_lst = easygui.multchoicebox(
            self.history+msg,
            self.title,
            choices=AAs,
            preselect=preselect)
        if AA_lst is None:
            self.history = f"{name} polymerisation not added - no AAs provided."
            return preselect
        AA_lst.append(None)
        P1 = {"pep_len": [3, 4, 5],
              4: {"allowed": AA_lst},
              5: {"allowed": AA_lst},
              "gly_len": [2],
              "gly_rejected": ["anMurNAc", "anMurN", "anMurNGlyc",  # anhydro
                               "MurN[r]", "MurNAc[r]", "MurNAcP[r]",  # reduced
                               "MurNGlyc[r]", "MurNAcOAc[r]"]}

        if prefer_reduced:
            gly_accepted = ["GlcNAc", "MurNAc[r]", "anMurNAc"]
        else:
            gly_accepted = ["GlcNAc", "MurNAc", "anMurNAc"]

        P2 = {"pep_len": [3, 4, 5],
              4: {"allowed": AA_lst},
              5: {"allowed": AA_lst},
              "gly_len": [2],
              "gly_accepted": gly_accepted}
        P1_bond = [0, "glycan", "glycan", "Acc"]
        P2_bond = [0, "glycan", "glycan", "Dnr"]

        self.history = self.generator.set_polymerisation_types(P1, P2,
                                                               P1_bond, P2_bond)
        preselect = get_preselect_integers(AAs, AA_lst)
        return preselect

    def gen_add_default_peptide_polymerisation(self, name, length_p2,
                                               bridge, prefer_reduced,
                                               preselect=0):
        """
        Adds default peptide crosslinks, assuming the acyl acceptor is the
        third amino acid on P1, with or without a bridge.


        Parameters
        ----------
        name : string
            Name of polymerisation.
        length_p2 : integer
            Position of the amino acid that acts as an acyl donor on P2.
        bridge : boolean
            If True, the peptide bridges are required on the acyl acceptor.
        prefer_reduced : boolean, optional
            If True, polymerisations will prefer the reduced Mur-type sugars
            over non-reduced sugars. The default is False.
        preselect : list, optional
            List of integers that correspond to AAs. These will appear as the
            default option. The default is 0 (No amino acids).

        Returns
        -------
        preselect : list, optional
            List of integers that correspond to AAs. These will appear as the
            default option in future selections.

        """
        msg = f"Set acceptable AAs in position 4 and 5 in the stem peptide\
            of {name} polymerisation."
        AAs = self.gen_retrieve_AAs(idx=[4, 5])
        if self.generator.modifications["Alanine/Lactate Substitution"]:
            AAs.append("Lac")  # add Lactate if modification used
        AA_lst = easygui.multchoicebox(
            self.history+msg,
            self.title,
            choices=AAs,
            preselect=preselect)
        if AA_lst is None:
            self.history = f"{name} polymerisation not added - no AAs provided."
            return preselect
        AA_lst.append(None)
        P1 = {3: {"bridge": bridge},
              4: {"allowed": AA_lst},
              5: {"allowed": AA_lst},
              "pep_len": [3, 4, 5],
              "gly_len": [0, 2]}
        P2 = {"pep_len": [length_p2],
              "gly_len": [0, 2]}
        for i in range(4, 6):
            if i <= length_p2:
                P2[i] = {"allowed": AA_lst}

        if bridge:
            P1_bond = [3, "bridge", "eupeptide", "NH2"]
        else:
            P1_bond = [3, "main", "side", "NH2"]
        P2_bond = [length_p2, "main", "eupeptide", "COOH"]

        self.history = self.generator.set_polymerisation_types(P1, P2,
                                                               P1_bond, P2_bond)
        preselect = get_preselect_integers(AAs, AA_lst)
        return preselect

    # %%%% Diffcalc

    def gen_set_diffcalc_units(self):
        """
        Selects the difference calculator overall parameters:
            Min. polymerisation no.: Reduce differences only when no.
            polymerisations exceeds this threshold.
            Max. differences: Remove monomers which exceed this much difference.

        Returns
        -------
        None.

        """
        msg = '''Set diffcalc parameters.
        Min. polymerisation no.: Reduce only when no. polymerisations exceeds this.
        Max. differences: Remove monomers which exceed this much difference.'''
        fields = ["Min. polymerisation no.", "Max. differences"]
        field_values = [1,
                        1]
        mod_values = easygui.multenterbox(self.history+msg, self.title,
                                          fields=fields, values=field_values)
        min_cl = int(mod_values[0])
        max_diffs = mod_values[1].split(",")
        max_diffs = map(float,max_diffs)
        self.history = self.generator.set_diffcalc_units(min_cl, max_diffs)

    def gen_set_diffcalc_params(self,
                                gly_len=None, gly_range=None,
                                pep_len=None, pep_range=None):
        """
        Selects the difference calculator parameters if parameter currently
        is None. Selects acceptable glycan and peptide lengths, acceptable glycans,
        amino acids and bridge peptides.
        !!! To amend: Acceptable polymerisation types can't be chosen

        Parameters
        ----------
        gly_len : list, optional
            List of acceptable glycan lengths. The default is None.
        gly_range : list, optional
            List of acceptable glycans. The default is None.
        pep_len : list, optional
            List of acceptable peptide lengths. The default is None.
        pep_range : list, optional
            List of acceptable amino acids in stem peptide. The default is None.

        Returns
        -------
        None.

        """
        try:
            if gly_len is None:
                gly_len = self.gen_set_diffcalc_gly_len()
            if gly_range is None:
                gly_range = self.gen_set_diffcalc_gly_range()
            if pep_len is None:
                pep_len = self.gen_set_diffcalc_pep_len()
            if pep_range is None:
                pep_range = self.gen_set_diffcalc_pep_range()
            bridge_range = self.gen_set_diffcalc_bridge_range()
            polymerisation_range = self.gen_set_diffcalc_poly_range()
            self.history = self.generator.set_diffcalc_param(gly_len, gly_range,
                                                             pep_len, pep_range,
                                                             bridge_range,
                                                             polymerisation_range)
        except Exception:
            self.register_step(False)
        finally:
            self.register_step(True)

    def gen_set_diffcalc_gly_len(self):
        """
        Selects acceptable glycan lengths.

        Returns
        -------
        gly_len : list, optional
            List of acceptable glycan lengths.

        """
        msg = "Set canonical number of glycans per monomer."
        if self.generator.glycan_lmax > 0:
            max_gly_len = self.generator.glycan_lmax
            max_gly_len = 2
        preselect = self.generator.diffcalc_params.get("gly_len", 0)
        choices = range(0, max_gly_len+1)
        preselect = get_preselect_integers(choices, preselect)
        gly_len = easygui.multchoicebox(
            self.history+msg,
            self.title,
            choices=choices,
            preselect=preselect)
        gly_len = [int(x) for x in gly_len]
        return gly_len

    def gen_set_diffcalc_gly_range(self):
        """
        Selects list of acceptable glycans.

        Returns
        -------
        gly_range : list, optional
            List of acceptable glycans.

        """
        msg = "Set canonical glycans in monomer."
        glycans = []
        for t in ("Glc", "Mur"):
            if t in self.generator.glycan_units:
                glycans.extend(self.generator.glycan_units[t])
            else:
                glycans.extend(Generator.valid_glycans[t])
        glycans = list(set(glycans))
        preselect = self.generator.diffcalc_params.get("gly_range", 0)
        preselect = get_preselect_integers(glycans, preselect)
        gly_range = easygui.multchoicebox(
            self.history+msg,
            self.title,
            choices=glycans,
            preselect=preselect)
        gly_range.append(None)
        return gly_range

    def gen_set_diffcalc_pep_len(self):
        """
        Selects acceptable peptide lengths.

        Returns
        -------
        pep_len : list, optional
            List of acceptable peptide lengths.

        """
        msg = "Set canonical number of AAs in stem peptide of monomer."
        if self.generator.peptide_lmax > 0:
            max_pep_len = self.generator.peptide_lmax
        else:
            max_pep_len = 5
        preselect = self.generator.diffcalc_params.get("pep_len", 0)
        choices = range(0, max_pep_len+1)
        preselect = get_preselect_integers(choices, preselect)
        pep_len = easygui.multchoicebox(
            self.history+msg,
            self.title,
            choices=choices,
            preselect=preselect)
        pep_len = [int(x) for x in pep_len]
        return pep_len

    def gen_set_diffcalc_pep_range(self):
        """
        Selects acceptable amino acids in stem peptide.

        Returns
        -------
        pep_range : list, optional
            List of acceptable amino acids in stem peptide.

        """
        msg = "Set canonical AAs in stem peptide of monomer."
        AAs = self.gen_retrieve_AAs()
        preselect = self.generator.diffcalc_params.get("pep_range", 0)
        preselect = get_preselect_integers(AAs, preselect)
        pep_range = easygui.multchoicebox(
            self.history+msg,
            self.title,
            choices=AAs,
            preselect=preselect)
        if pep_range is not None:
            pep_range.append(None)
        else:
            pep_range = []
        return pep_range

    def gen_set_diffcalc_bridge_range(self):
        """
        Selects acceptable bridge peptides.

        Returns
        -------
        bridge_range : list
            List of acceptable bridge peptides.

        """
        if len(self.generator.bridge_peptides) == 0:
            return None
        bridges = []
        for idx in self.generator.bridge_peptides:
            for grp in self.generator.bridge_peptides[idx]:
                bridges.extend(
                    self.generator.bridge_peptides[idx][grp]["bridges"])
        if len(bridges) < 2:
            return None
        msg = "Set canonical bridge(s) in monomer."
        mapping = {str(i): x for i, x in enumerate(bridges)}
        msg += f"\n{mapping}"
        bridge_range = easygui.multchoicebox(
            self.history+msg,
            self.title,
            choices=mapping.keys())
        bridge_range = [mapping[k] for k in bridge_range]
        return bridge_range

    def gen_set_diffcalc_poly_range(self):
        """
        !!! To amend: Acceptable polymerisation types can't be chosen

        Returns
        -------
        None.

        """
        return None

    # %%%% Run

    def gen_execute(self):
        """
        Runs Generator and MSPMaker in stages. Four options are available:
            "Muropeptide Library":  Creates muropeptide library.
            "Spectral Library":     Creates spectral library. Fails if muropeptide
                                    library is not run.
            "Run with no outputs":  Runs but no output is created. This option
                                    is for troubleshooting only.
            "Close after run":      Closes the GUI window when complete.
                                    !!!Does this work?

        Returns
        -------
        boolean
            If True, closes the GUI window.

        """
        msg = "What would you like to generate?"
        options = ["Muropeptide Library",
                   "Spectral Library",
                   "Run with no outputs",
                   "Close after run"]
        run_options = easygui.multchoicebox(self.history+msg,
                                            self.title,
                                            choices=options,
                                            preselect=[0, 1])
        pkl = None
        if run_options is None:
            return True
        else:
            output = "Run with no outputs" not in run_options
            if not output:
                print("No outputs will be saved!")
            print("\n####### RUNNING.... #######\n")
            if "Muropeptide Library" in run_options:
                pkl = self.gen_create_muropeptide_library(save_output=output)
            if "Spectral Library" in run_options:
                self.gen_create_spectral_library(pkl, save_output=output)
            print("\n####### COMPLETED! #######\n")
            return "Close after run" not in run_options

    def gen_create_muropeptide_library(self, save_output=True):
        """
        Creates the muropeptide library by running Generator.

        Parameters
        ----------
        save_output : boolean, optional
            If True, saves all outputs. The default is True.

        Returns
        -------
        pkl : pickle
            Pickle file with chemical structures of PGN monomers and dimers.
            PGN trimers and above are excluded as their spectra takes a long
            time to generate and are generally not useful for identification.

        """
        gen = self.generator
        self.history = gen.generate_PGN()
        pkl = None
        if save_output:
            gen.export_settings_as_yaml()
            gen.export_settings_as_image()
            pkl = gen.export_PGN_as_pickle(kmers=[1, 2])
            self.history = f"Pickle file saved at: {pkl}"
            gen.export_dataframes(["xlsx"])
            gen.export_dataframes(["txt"], kmers=[3, 4])
        return pkl

    def gen_create_spectral_library(self, pkl,
                                    adducts=None,
                                    export="all",
                                    save_output=True):
        """
        Creates a spectral library. PGN trimers and above are excluded by default
        as their spectra takes a long time to generate and are generally not useful
        for identification.

        Parameters
        ----------
        pkl : pickle
            Pickle file containing chemical structures of PGN
        adducts : list, optional
            List of adducts to include. If None, the default adducts are used:
            '[M+H]+', '[M+2H]2+', '[M+3H]3+'
        export : list, optional
            List of strings that describe output files.
            "pickle" : Pickled compound data (fragments)
            "msp" : mspfile
            "peaklist" : Detailed compound data in xlsx.
            "ionlist" : Consolidated list of ions for all compounds in xlsx.
            "all" : All outputs.
            The default is "all".
        save_output : boolean, optional
            if True, saves all outputs. The default is True.

        Returns
        -------
        None.

        """
        if adducts is None:
            adducts = ['[M+H]+', '[M+2H]2+', '[M+3H]3+']
        if pkl is None:
            msg = "Pick a pkl file."
            pkl = easygui.fileopenbox(msg,
                                      self.title,
                                      default="/output/*.pickle",
                                      filetypes=".pickle"
                                      )

        maker = MSPMaker(output_adducts=adducts)
        if save_output:
            maker.import_PGN_as_pickle(pkl, export=export)
        else:
            maker.import_PGN_as_pickle(pkl, export=[])

    # %%% Fragmenter

    def frag_single_compound(self):
        """
        Fragments a single compound.

        Returns
        -------
        None.

        """
        mol_name, molecule = self.frag_get_molecule_details()
        if molecule is None:
            return
        msg = "Pick adduct types."
        adducts = ['[M+H]+', '[M+2H]2+', '[M+3H]3+']
        adduct = easygui.multchoicebox(self.history+msg,
                                       self.title,
                                       choices=adducts)
        if adduct is None:
            return
        else:
            if mol_name is None:
                mol_name = molecule.InchiKey
            frag = Fragmenter("single_cpd")
            for a in adduct:
                self.history = f"Fragmenting {mol_name}_{a}"
                frag.load_PGN(molecule, mol_name,
                              adduct_type=a, auto_process=True)
            os.startfile(frag.dir)

    def frag_get_molecule_details(self):
        """
        Gets details about the molecule from the user.

        Returns
        -------
        name : string
            Name of molecule. If not provided by user, uses the InChlKey.
        molecule : string
            SMILES string describing molecule.

        """
        choice = 0
        smiles = None
        name = None
        molecule = None
        options = ["Enter SMILES", "Clear SMILES",
                   "Enter Name", "Clear Name",
                   "View Molecule", "Return"]
        while choice not in ["Return", None]:
            msg = "Choose an option."
            choice = easygui.choicebox(self.history+msg,
                                       self.title,
                                       choices=options)
            if choice == "Enter SMILES":
                msg = "Enter SMILES."
                smiles = easygui.textbox(self.history+msg,
                                         self.title)
                molecule = Molecule(smiles)
                self.history = f"Added {molecule.InchiKey} with formula: {molecule.formula}"
            elif choice == "Clear SMILES":
                smiles = None
                molecule = None
                self.history = "Cleared SMILES"
            elif choice == "View Molecule" and molecule is not None:
                img = molecule.draw()
                key = molecule.InchiKey
                img_path = f"img/{key}.png"
                img.save(img_path)
                msg = f'''This is your molecule:
                    Formula: {molecule.formula}
                    InChiKey: {molecule.InchiKey}'''
                choices = ["I see!"]
                easygui.buttonbox(msg, image=img_path, choices=choices)
            elif choice == "Enter Name":
                msg = "Add a name (optional)."
                name = easygui.textbox(self.history+msg,
                                       self.title)
                if name is not None:
                    self.history = f"Added name: {name}"
            elif choice == "Clear Name":
                name = None
        return name, molecule


if __name__ == "__main__":
    myGUI = GUI()
    myGUI.start()
