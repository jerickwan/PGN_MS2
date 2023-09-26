'''
Outputs a diagram of Generator
env = chem
'''
import matplotlib.pyplot as pyplot
import copy
from math import cos, sin, pi
from matplotlib.colors import to_rgba


ILL_PARAMS = {
    "COLORS": {
        "GLU": "blue", "MUR": "green",
        "AA": "black", "C": "red", "N": "blue"},
    "SIZES": {
        "A5": (5.8, 8.3), "A4": (8.3, 11.7)},
    "RADII": {
        "GLYCAN": 0.09, "PEPTIDE": 0.025, "GAP": 0.0175}
}


class Illustrator:
    '''Generates a graphical summary of the generated PGN.'''

    def __init__(self, size, fontsize, bond_char):
        '''Sets the size ("A4","A5"), fontsize and bonding character.'''
        self.height = ILL_PARAMS["SIZES"][size][0]
        self.width = ILL_PARAMS["SIZES"][size][1]
        self.fontsize = fontsize
        self.max_vertical_size = 5
        self.title_multiplier = 1.25
        self.subtitle_multiplier = 0.7
        self.fontsize_title = self.fontsize*self.title_multiplier
        self.fontsize_subtitle = self.fontsize*self.subtitle_multiplier
        self.fig, self.ax = pyplot.subplots(figsize=[self.width, self.height],
                                            dpi=300)
        self.ax.tick_params(bottom=False, top=False, left=False, right=False)
        self.ax.tick_params(labelbottom=False, labeltop=False,
                            labelleft=False, labelright=False)
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)
        # Size of Plot
        self.glycan_radius = ILL_PARAMS["RADII"]["GLYCAN"]
        self.peptide_radius = ILL_PARAMS["RADII"]["PEPTIDE"]  # ideal for 5
        self.peptide_wrap_width = 70  # ideal for 2 polymerisations
        self.gap = ILL_PARAMS["RADII"]["GAP"]
        self.AA_xs = None
        self.AA_ys = None
        self.pep_x = None
        self.pep_y = None
        self.bond_char = bond_char

    def adjust_sizes_for_peptide_length(self, pep_length, ideal_length=5):
        """
        Shrinks peptide radius, gap and fontsize if maximum peptide length
        exceeds 5.

        Parameters
        ----------
        pep_length : integer
            Maximum peptide length used in library generation.
        ideal_length : integer, optional
            Ideal peptide length. The default is 5.

        Returns
        -------
        None.

        """
        if pep_length <= ideal_length:
            pass
        else:
            adjustment = 1 - 0.05*(pep_length-ideal_length)
            self.peptide_radius *= adjustment
            self.gap *= adjustment
            self.fontsize *= adjustment

    def adjust_sizes_for_polymerisation_types(self, polymer_types, ideal_no=2):
        """
        Shrinks the text wrap width if the number of polymerisations exceeds 2.

        Parameters
        ----------
        polymer_types : list
            List of polymerisations used in library generation.
        ideal_no : TYPE, optional
            Ideal number of polymerisations. The default is 2.

        Returns
        -------
        None.

        """
        polymer_types_no = len(polymer_types)
        if polymer_types_no <= ideal_no:
            self.peptide_wrap_width *= 0.5
        else:
            self.peptide_wrap_width *= 0.2

    # %% Drawing

    def draw_polygon(self, ax, num_sides, x, y, r, ec, fc, alpha=1, init_angle=0):
        """
        Draws n-gon (num_sides) on ax with centre (x,y), radius (r),
        edge color (ec) and face color (fc).

        Parameters
        ----------
        ax : matplotlib.Ax
            Axes of shape
        num_sides : integer
            Number of sides.
        x : float
           x position as a fraction
        y : float
            y position as a fraction
        r : float
            Radii of shape
        ec : string
            matplotlib color for shape edge
        fc : string
            matplotlib color for shape face
        alpha : float, optional
            Transparency of shape. The default is 1.
        init_angle : float, optional
            Initial angle of shape (in degrees). The default is 0.

        Returns
        -------
        None.

        """
        xs = [x+r*cos(init_angle+n*2*pi/num_sides)
              for n in range(0, num_sides+1)]
        ys = [y+r*sin(init_angle+n*2*pi/num_sides)
              for n in range(0, num_sides+1)]
        fc = self.add_alpha(fc, alpha)
        ax.fill(xs, ys, ec=ec, fc=fc, lw=2, zorder=2)

    def draw_rectangle(self, ax, x, y, hw, hh, ec):
        """
        Draws an unfilled rectangle on ax with centre (x,y), half-width (hw),
        half-height(hh) and edge color (ec).

        Parameters
        ----------
        ax : matplotlib.Ax
            Axes of shape
        num_sides : integer
            Number of sides.
        x : float
           x position as a fraction
        y : float
            y position as a fraction
        hw : float
            Half-width of shape
        hh : float
            Half-height of shape
        ec : string
            matplotlib color for shape edge

        Returns
        -------
        None.

        """
        ''''''
        xs = [x+hw, x+hw, x-hw, x-hw]
        ys = [y+hh, y-hh, y-hh, y+hh]
        ax.fill(xs, ys, fill=False, edgecolor=ec, lw=2, zorder=2)

    # %% Formatting

    def add_alpha(self, c, alpha):
        """
        Adds alpha attribute to color (c).

        Parameters
        ----------
        c : string
            matplotlib color
        alpha : float
            Transparency of shape. The default is 1.

        Returns
        -------
        rgba : tuple
            Tuple describing color in RGBA format.

        """
        rgba = to_rgba(c, alpha)
        return rgba

    def highlight_text(self, text_list, highlights_list):
        """
        Highlights text in text_list that are present in highlights_list by
        adding suitable formatting codes.

        Parameters
        ----------
        text_list : list
            List of strings.
        highlights_list : list
            List of strings to be highlighted.

        Returns
        -------
        text_list : list
            Modified list with suitable formatting codes.

        """
        if None in highlights_list:
            highlights_list.remove(None)
        replacements_list = [r"$\bf{" + hl + "}$" for hl in highlights_list]
        for hl, rep in zip(highlights_list, replacements_list):
            if hl in text_list:
                idx = text_list.index(hl)
                text_list.pop(idx)
                text_list.insert(0, rep)  # shift to front
        return text_list

    def hor_text(self, text_list, wrap_width=70):
        """
        Formats text horizontally.
        Wraps text if number of characters exceed threshold (wrap_width).

        Parameters
        ----------
        text_list : list
            List of strings.
        wrap_width : integer, optional
            Maximum character length allowed per line. The default is 70.

        Returns
        -------
        txt : string
            Formatted and wrapped text.

        """
        txt = ", ".join(str(x) for x in text_list)
        pure_txt = txt
        formatting_strings = ["$\bf{", "}$"]
        for format_string in formatting_strings:
            if format_string in pure_txt:
                pure_txt = pure_txt.replace(format_string, "")
        if len(pure_txt) > wrap_width:
            divider = len(text_list)//2
            wrap_widthped = "\n".join((self.hor_text(text_list[:divider]),
                                       self.hor_text(text_list[divider:])))
            return wrap_widthped
        return txt

    def vert_text(self, text_list):
        """
        Formats text vertically.

        Parameters
        ----------
        text_list : list
            List of strings.

        Returns
        -------
        txt : string
            Formatted and wrapped text.

        """
        return "\n".join(str(x) for x in text_list)

    def get_fontsize(self, orient, length):
        """
        Adjusts fontsize if vertical text exceeds maximum vertical size.

        Parameters
        ----------
        orient : string
            Orientation of text - "vert" or "horizontal".
        length : integer
            Length of text.

        Returns
        -------
        float
            Fontsize.

        """
        if orient == "vert" and length > self.max_vertical_size:
            return self.fontsize*self.max_vertical_size/(length+0.5)
        else:
            return self.fontsize

    def format_glycan_range_text(self, gly_lst, wrap_width=20):
        """
        Formats and summarised range of glycans.
        i.e. (GlcNAc,GlcNAcOAc, GlcN) --> Glc: NAc, NAcOAc, N
        ***Not in use.

        Parameters
        ----------
        gly_lst : list
            List of glucsoamine/muramic acid sugars - each represented by a string.
        wrap_width : integer, optional
            Maximum character length allowed per line. The default is 20.

        Returns
        -------
        string
            Summarised text showing range of glycans.

        """
        gly_lst = copy.deepcopy(gly_lst)
        txts = []
        if None in gly_lst:
            txts.append("None")
            gly_lst.remove(None)
        for sugar in "Glc", "Mur":
            selected = [gly for gly in gly_lst if sugar in gly]
            summarised_txt = sorted([gly.replace(sugar, "")
                                    for gly in selected])
            summarised_txt = self.hor_text(summarised_txt, wrap_width)
            txts.append(f"{sugar}: {summarised_txt}")
        return "\n".join(txts)

    # %% Plotting

    def plot_sugars(self,
                    glc_lst=["GlcNAc"],
                    mur_lst=["MurNAc"],
                    lmin=0, lmax=2,
                    highlight_sugars=None):
        """
        Plots the glycan chain.

        Parameters
        ----------
        glc_lst : list, optional
            List of glucsoamine sugars - each represented by a string.
            The default is ["GlcNAc"].
        mur_lst : TYPE, optional
            List of muramic acid sugars - each represented by a string.
            The default is ["MurNAc"].
        lmin : integer, optional
            Minimum glycan length. The default is 0.
        lmax : integer, optional
            Maximum glycan length. The default is 2.
        highlight_sugars : list, optional
            List of sugars to highlight. The default is None.

        Returns
        -------
        None.

        """
        # Draw Glc
        glc_x = 0.5-2*self.glycan_radius-self.glycan_radius*cos(4*pi/3)
        glc_y = 0.85
        if highlight_sugars is not None:
            glc_lst = self.highlight_text(glc_lst, highlight_sugars)
        glc_text = self.vert_text(glc_lst)
        self.draw_polygon(self.ax, 6,
                          glc_x, glc_y,
                          self.glycan_radius,
                          "black",
                          ILL_PARAMS["COLORS"]["GLU"],
                          alpha=0.5)
        self.ax.text(glc_x, glc_y,
                     glc_text,
                     ha="center", va="center",
                     fontsize=self.get_fontsize("vert", len(glc_lst)))
        # Draw Mur
        mur_x, mur_y = glc_x+2*self.glycan_radius, glc_y
        if highlight_sugars is not None:
            mur_lst = self.highlight_text(mur_lst, highlight_sugars)
        mur_text = self.vert_text(mur_lst)
        self.draw_polygon(self.ax, 6,
                          mur_x, mur_y,
                          self.glycan_radius,
                          "black",
                          ILL_PARAMS["COLORS"]["MUR"],
                          alpha=0.5)
        self.ax.text(mur_x, mur_y,
                     mur_text,
                     ha="center", va="center",
                     fontsize=self.get_fontsize("vert", len(mur_lst)))
        # Labels
        len_text = f"Glycan length: G{lmin} - G{lmax}"
        self.ax.text(glc_x-self.glycan_radius-self.gap, mur_y,
                     len_text,
                     ha="right", va="center",
                     fontsize=self.get_fontsize("vert", 1))
        # Title
        self.ax.text(glc_x, 1-self.gap,
                     "Glycan Chain",
                     ha="center", va="top",
                     fontsize=self.fontsize_title)

        self.pep_x = mur_x+self.glycan_radius*cos(4*pi/3)
        self.pep_y = mur_y+self.glycan_radius*sin(4*pi/3)

    def plot_peptide_chain(self, x0, y0, lmax, AA_dict,
                           label_ha="right",
                           descending=True, subtitle=False,
                           highlight_AAs=None,
                           wrap_width=100):
        """
        Plots an ascending/descending peptide chain of length (lmax) from
        initial position (x0,y0). Amino acids are labelled from AA_dict.

        Parameters
        ----------
        x0 : float
           initial x position as a fraction
        y0 : float
           initial y position as a fraction
        lmax : integer
            Maximum peptide length.
        AA_dict : Dictionary
            Dictionary which shows list of amino acids for each position.
        label_ha : string, optional
            Horizontal alignment for peptide labels.
            The default is "right".
        descending : boolean, optional
            If True, peptide is drawn in descending orientation.
            The default is True.
        subtitle : boolean, optional
            Determines fontsize of labels. If True, uses the smaller fontsize.
            The default is False.
        highlight_AAs : list, optional
            List of amino acids to be highlighted. The default is None.
        wrap_width : integer, optional
            Maximum character length allowed per line. The default is 100.

        Returns
        -------
        xs : list
            list of floats representing x position of each amino acid.
        ys : list
            list of floats representing y position of each amino acid.
        N_terminus : Lines.Line2D
            matplotlib line object representing N-terminus of peptide
        C_terminus : Lines.Line2D
            matplotlib line object representing C-terminus of peptide

        """
        xs, ys = [0]*lmax, [0]*lmax
        if subtitle:
            peptide_fs = self.fontsize_subtitle
        else:
            peptide_fs = self.fontsize
        for i in range(lmax):
            xs[i] = x0
            y_diff = (2*self.peptide_radius+self.gap)*i
            if descending:
                ys[i] = y0-y_diff
            else:
                ys[i] = y0+y_diff
            # Draw Square
            self.draw_polygon(self.ax, 4,
                              xs[i], ys[i],
                              self.peptide_radius,
                              "black", "white",
                              init_angle=pi/4)
            # Label Position
            self.ax.text(xs[i], ys[i],
                         i+1,
                         ha="center",
                         va="center",
                         fontsize=peptide_fs)
            # Label AAs
            AA_lst = AA_dict.get(i+1, [""])
            if highlight_AAs is not None:
                AA_lst = self.highlight_text(AA_lst, highlight_AAs)
            AA_txt = self.hor_text(AA_lst, wrap_width)
            if label_ha == "right":
                label_x = xs[i]-self.peptide_radius
            else:
                label_x = xs[i]+self.peptide_radius
            self.ax.text(label_x,
                         ys[i],
                         AA_txt,
                         ha=label_ha,
                         va="center",
                         wrap=True,
                         fontsize=peptide_fs)
        # Peptide Chain
        self.ax.vlines(x0, y0,
                       ys[-1],
                       color=ILL_PARAMS["COLORS"]["AA"],
                       zorder=1)
        # Terminals
        if descending:
            y1 = y0+(self.peptide_radius+self.gap)
        else:
            y1 = y0-(self.peptide_radius+self.gap)
        N_terminus = self.ax.vlines(x0, y1,
                                    y0,
                                    color=ILL_PARAMS["COLORS"]["N"],
                                    zorder=1)
        if descending:
            y1 = ys[-1]-(self.peptide_radius+self.gap)
        else:
            y1 = ys[-1]+(self.peptide_radius+self.gap)
        C_terminus = self.ax.vlines(x0, ys[-1],
                                    y1,
                                    color=ILL_PARAMS["COLORS"]["C"],
                                    zorder=1)
        return xs, ys, N_terminus, C_terminus

    def plot_stem_peptide(self, AA_dict, lmin=0, lmax=5, highlight_AAs=None):
        """
        Plots the stem peptide of PGN.

        Parameters
        ----------
        AA_dict : Dictionary
            Dictionary which shows list of amino acids for each position.
        lmin : integer, optional
            Minimum peptide length. The default is 0.
        lmax : integer, optional
            Maximum peptide length. The default is 5.
        highlight_AAs : list, optional
            List of amino acids to be highlighted. The default is None.

        Returns
        -------
        None.

        """
        pep_x = self.pep_x
        pep_y = self.pep_y-self.peptide_radius-self.gap
        self.AA_xs, self.AA_ys, _, _ = self.plot_peptide_chain(
            pep_x, pep_y, lmax, AA_dict,
            highlight_AAs=highlight_AAs,
            wrap_width=self.peptide_wrap_width)
        # Stem Peptide Length
        len_text = f"Stem length: P{lmin} - P{lmax}"
        self.ax.text(pep_x, self.AA_ys[-1]-self.peptide_radius-self.gap,
                     len_text,
                     ha="center", va="top",
                     fontsize=self.get_fontsize("vert", 1))
        # Stem Peptide Title
        self.ax.text(self.gap, 0.6,
                     "Stem Peptide",
                     va="center",
                     ha="left",
                     rotation=90,
                     fontsize=self.fontsize_title)

    def plot_bridge_peptides(self, bridge_peptides, highlight_bridges=None):
        """
        Plots the bridge peptides.

        Parameters
        ----------
        bridge_peptides : dictionary
            Dictionary object with indexes of bridge peptides as keys and
            lists of bridge peptides as entries.
        highlight_bridges : list, optional
            List of bridge peptides to be highlighted. The default is None.

        Returns
        -------
        None.

        """
        # Title
        if len(bridge_peptides) > 0:
            self.ax.text(1-self.gap, 0.6,
                         "Bridge Peptides",
                         va="center",
                         ha="right",
                         rotation=270,
                         fontsize=self.fontsize_title)

        for idx in bridge_peptides:
            for grp in bridge_peptides[idx]:
                x1 = self.AA_xs[idx-1]+self.gap+self.peptide_radius
                y1 = self.AA_ys[idx-1]
                if grp == "NH2":
                    y1 -= self.peptide_radius
                    color = "blue"
                    bridges = [self.bond_char.join(x)
                               for x in bridge_peptides[idx][grp]['bridges']]
                else:
                    y1 += self.peptide_radius
                    color = "red"
                    bridges = [self.bond_char.join(x)
                               for x in reversed(bridge_peptides[idx][grp]['bridges'])]
                if highlight_bridges is not None:
                    highlights = [self.bond_char.join(x)
                                  for x in highlight_bridges]
                    bridges = self.highlight_text(bridges, highlights)

                bridges_txt = self.hor_text(bridges)
                self.ax.text(x1, y1,
                             bridges_txt,
                             ha="left",
                             va="center",
                             wrap=True,
                             fontsize=self.fontsize_subtitle)
                x0, y0 = self.AA_xs[idx-1], self.AA_ys[idx-1]
                self.ax.plot([x0, x1], [y0, y1],
                             color=color,
                             ls="--",
                             zorder=1)

    def plot_polymerisations(self, polymer_types, num_polymers):
        """
        Plots the polymerisation types (glycosidic or peptide crosslinks).

        Parameters
        ----------
        polymer_types : list
            List of lists; with each list containing the polymerisation name,
            type and conditions.
            See Generator.polymerisation_types for more details.
        num_polymers : integer
            Total number of polymerisations.

        Returns
        -------
        None.

        """
        total = len(polymer_types)
        if total == 0 or num_polymers == 0:
            # Polymerisation Num. text
            self.ax.text(0.5, self.gap,
                         f"No. Polymerizations: None",
                         va="bottom",
                         ha="center",
                         fontsize=self.fontsize_title)
            return
        # Determines the location of each polymerisation illustration
        origins = [i/(2*total) for i in range(1, total*2, 2)]
        center_y = self.AA_ys[-1]-5*self.peptide_radius
        self.gap /= 2
        self.peptide_radius *= self.subtitle_multiplier
        # Polymerisation Num. text
        self.ax.text(0.5, self.gap,
                     f"No. Polymerizations: 0-{num_polymers}",
                     va="bottom",
                     ha="center",
                     fontsize=self.fontsize_title)
        # Plots each polymerisation type
        for i, info in enumerate(polymer_types):
            # n_cond: conditions for linked chains
            # n_bond: conditions for crpss-link bond
            name, P1_cond, P2_cond, P1_bond, P2_bond, kmer_range = info
            center_x = origins[i]
            self.ax.text(center_x, center_y+self.gap,
                         name,
                         va="center",
                         ha="center",
                         fontsize=self.fontsize_subtitle)
            # Plotting parameters for N and C terminus
            both_x = []  # x position
            both_y = []  # y position
            colors = []  # red for C; blue for N
            styles = []  # solid for eupeptide; dashed for side chain
            special = False  # peptide bridge or glycosidic
            # iterate over P1, P2
            for cond_dict, bond in ((P1_cond, P1_bond), (P2_cond, P2_bond)):
                AA_dict = {i: cond_dict[i]["allowed"] for i in cond_dict
                           if "allowed" in cond_dict[i]}
                # shift P1 to the left
                if bond == P1_bond:
                    center_x -= 2.5*self.peptide_radius
                    label_ha = "right"
                # shift P2 to the right
                else:
                    center_x += 5*self.peptide_radius
                    label_ha = "left"
                bond_pos, chain, rxn_type, grp = bond
                # "pep_len": peptide_length is necessary
                lmin = min(cond_dict["pep_len"])
                lmax = max(cond_dict["pep_len"])
                xs, ys, N_terminus, C_terminus = self.plot_peptide_chain(
                    center_x,
                    center_y,
                    lmax,
                    AA_dict,
                    label_ha=label_ha,
                    subtitle=True,
                    wrap_width=self.peptide_wrap_width)
                # Setting parameters for terminus
                both_x.append(xs[bond_pos-1])
                both_y.append(ys[bond_pos-1])
                if rxn_type == "eupeptide":
                    styles.append("-")
                else:
                    styles.append("--")
                if grp == "NH2":
                    colors.append("blue")
                    if rxn_type == "eupeptide":
                        N_terminus.set(visible=False)
                elif grp == "COOH":
                    colors.append("red")
                    if rxn_type == "eupeptide":
                        C_terminus.set(visible=False)
                else:
                    colors.append("black")
                if chain == "bridge":
                    special = "bridged"
                elif chain == "glycan":
                    special = "glyco"
                # Peptide Length Label
                current_y = ys[-1]
                current_y -= (self.peptide_radius+self.gap)
                len_text = f"P{lmin} - P{lmax}"
                self.ax.text(xs[0], current_y,
                             len_text,
                             ha="center", va="top",
                             fontsize=self.fontsize_subtitle)
                # Glycan Length Label
                if "gly_len" in cond_dict:
                    current_y -= (self.peptide_radius+self.gap)
                    gmin = min(cond_dict["gly_len"])
                    gmax = max(cond_dict["gly_len"])
                    glen_text = f"G{gmin} - G{gmax}"
                    self.ax.text(xs[0], current_y,
                                 glen_text,
                                 ha="center", va="top",
                                 fontsize=self.fontsize_subtitle)
                # Glycan Allowed Label
                # this is ugly
                # if "gly_allowed" in cond_dict:
                #     current_y -= (self.peptide_radius+self.gap)
                #     gly_lst = cond_dict["gly_allowed"]
                #     gly_text = self.format_glycan_range_text(gly_lst)
                #     self.ax.text(xs[0], current_y,
                #                  gly_text,
                #                  ha=label_ha, va="top",
                #                  wrap=True,
                #                  fontsize=self.fontsize_subtitle)

            # Determine meeting point of terminals
            bond_x = sum(both_x)/2
            bond_y = sum(both_y)/2
            # Bridge label
            if special:
                self.ax.text(bond_x, bond_y+self.peptide_radius,
                             special,
                             va="center",
                             ha="center",
                             fontsize=self.fontsize_subtitle)
            # Plot N and C terminus
            if special == "glyco":
                continue  # no need to plot
            for j in range(2):
                self.ax.plot([both_x[j], bond_x],
                             [both_y[j], bond_y],
                             ls=styles[j],
                             color=colors[j],
                             zorder=1)

    def plot_modifications(self, modifications):
        """
        Plots a list of modifications.

        Parameters
        ----------
        modifications : list
            List of strings representing each modification.

        Returns
        -------
        None.

        """
        lst_mods = [mod for mod in modifications if modifications[mod]]
        if len(lst_mods) > 0:
            mods_txt = self.vert_text(lst_mods)
            # Title
            self.ax.text(1-self.gap,
                         1-self.gap,
                         "Modifications",
                         va="top",
                         ha="right",
                         fontsize=self.fontsize_title)
            # Mods
            self.ax.text(1-self.gap,
                         1-self.glycan_radius,
                         mods_txt,
                         va="top",
                         ha="right",
                         fontsize=self.get_fontsize("vert", len(lst_mods)))

    def plot_all(self, gen):
        """
        Plots all components of PGN created by Generator (gen)

        Parameters
        ----------
        gen : Generator
            Generator object.

        Returns
        -------
        None.

        """
        self.adjust_sizes_for_peptide_length(gen.peptide_lmax)
        self.plot_sugars(gen.glycan_units["Glc"],
                         gen.glycan_units["Mur"],
                         gen.glycan_lmin,
                         gen.glycan_lmax,
                         highlight_sugars=gen.diffcalc_params.get("gly_range", None))
        self.plot_stem_peptide(gen.get_peptide_residues(),
                               gen.peptide_lmin,
                               gen.peptide_lmax,
                               highlight_AAs=gen.diffcalc_params.get("pep_range", None))
        self.plot_bridge_peptides(gen.bridge_peptides,
                                  highlight_bridges=gen.diffcalc_params.get("side_chain_range", None))
        self.plot_modifications(gen.modifications)
        self.adjust_sizes_for_polymerisation_types(gen.polymerisation_types)
        self.plot_polymerisations(
            gen.polymerisation_types, gen.num_polymerisations)

    # %% Export

    def save_figure(self, filename):
        """
        Saves generated plot.

        Parameters
        ----------
        filename : Path
            Save path.

        Returns
        -------
        None.

        """
        self.fig.savefig(filename, dpi=300, bbox_inches="tight")
