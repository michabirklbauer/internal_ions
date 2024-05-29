# Standard Library Imports
import json
import math
import itertools
from itertools import combinations, chain
from collections import Counter

# Third-party Library Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import rel_entr, softmax
from numpy.linalg import norm
from sklearn.mixture import GaussianMixture
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
from sklearn.metrics.pairwise import cosine_similarity
import plotly.graph_objects as go
import networkx as nx
from pyteomics import mass, parser as pyteomics_parser
import brainpy as bp
import fraggraph.constant as constant
import psm_utils
from scipy.optimize import minimize
from time import sleep
from tqdm import tqdm
import plotly.express as px

class FragGraph(nx.DiGraph):
    default_node_attributes = {
        # fragment attributes
        "frag_codes": [],
        "peptidoforms": [],
        "node_type": None,
        "start_ioncap": "t",
        "end_ioncap": "t",
        "start_pos": None,
        "end_pos": None,
        "frag_dir": None,
        "mz": None,
        "charge": -1,
        "isotope": -1,
        "isotope_prob": None,
        "intensity": 0,
        "depth": None,
        "neutral_loss": "",
        "kl_divergence": None,
        "cosine_similarity": None,
        # peak attributes
        "its_stdev": None,
        "its_mean": None,
        "its_median": None,
        "n_peaks": None,
        "weight": 1,
        # # visualisation attributes
        "color": "white",
        "size": 0,
        "y": 0,
        "x": 0,
        "physics": True,
        "title": "",
    }

    def __init__(self, **kwargs):
        super().__init__()

        # # read parameters from json file and add as attributes
        # frag_params = json.load(open("./fraggraph/fragmentation_parameters.json"))
        # if fragmentation_parameters in frag_params:
        #     frag_param = frag_params[fragmentation_parameters]
        #     # add attributes
        #     self.__dict__.update(frag_param)
        # else:
        #     raise NameError(
        #         f"Fragmentation parameters '{fragmentation_parameters}' not found in json file."
        #     )

        # parameters
        self.terminals_only = False
        self.min_weight_threshold = 0.01
        self.min_intermediate_2_evidence = 1
        self.min_intermediate_1_evidence = 1
        self.min_it_threshold_factor = (
            1  # Factor multiplied to the minimal itensity threshold
        )
        self.min_cosine_similarity = (
            0.7  # minimum cosine similarity to keep a intermediate 2 node
        )
        self.I1_rounding = (
            3  # n decimals in intermediate 1 nodes (determines node merging)
        )
        self.internal_over_terminal = (
            False  # if True, internal fragments can overlap with terminal fragments
        )
        self.fake_peak_factor = 100  # "fake peaks" are set to min_it / fake_peak_factor
        self.min_frag_length = 3  # minimum fragment length to consider
        self.filter_charge = False
        self.max_overlap_group_size = 250  # for overlaping group above this size, internal ions are discarded and only terminal ions are kept

        # Constants
        self.spectrum_node_radius = 10055
        self.circle_radius_leaf = 10000
        self.circle_radius_intermediate_0_I = 4000
        self.circle_radius_intermediate_0_C = 3000
        self.circle_radius_intermediate_0_N = 2000
        self.circle_radius_intermediate_1 = 6000
        self.peak_scaling_factor = 1000

        # viz constants
        self.colors = {
            "root": "red",
            "intermediate_0": "black",
            "intermediate_1_C": "#ed9aba",
            "intermediate_1_N": "#ed9aba",
            "intermediate_1_I": "#d0b5f5",
            "intermediate_2": "#ab8ef5",
            "leaf_first_isotope": "#ffbf00",
            "leaf_other_isotopes": "#f2d8a0",
            "leaf_not_matched": "#77abe6",
            "peak_fake": "#77abe6",
            "peak": "darkgrey",
        }

        # update kwargs
        self.__dict__.update(kwargs)

        # display the paramters in a printed table
        # print("Fragmentation parameters:")
        # for key, value in self.__dict__.items():
        #     print(f"{key}: {value}")

        # store I_1 nodes mz:
        self.I_1_nodes_mz = [-1000]

        # matching success
        self.matching_success = pd.DataFrame(
            columns=[
                "fragment_code",
                "intensity",
                "cosine_similarity",
                "kl_divergence",
                "charge",
                "length",
                "cterm_intens",
                "nterm_intens",
                "cterm_intens_diff",
                "nterm_intens_diff",
            ]
        )

    def generate_graph(self, peptidoforms, mz_list, intensities_list):
        # convert peptidoforms to psm_utils.Peptidoform objects
        peptidoforms = [
            psm_utils.Peptidoform(peptidoform) for peptidoform in peptidoforms
        ]
        self.peptidoforms = peptidoforms

        # Determine the maximum charge for the precursor peptide
        if self.max_prec_charge == "auto":
            # if peptidoform.precursor charge is not an integer
            if not all(
                isinstance(peptidoform.precursor_charge, int)
                for peptidoform in peptidoforms
            ):
                raise ValueError(
                    "Could not determine maximum precursor charge. Please specify it manually."
                )
            max_charge = max(
                peptidoform.precursor_charge for peptidoform in peptidoforms
            )

        else:
            max_charge = self.max_prec_charge

        # Calculate the maximum charge considering charge loss
        self.max_charge = max_charge - 1 if self.charge_loss else max_charge

        # Store mz_list and intensities_list
        self.mzs = mz_list
        self.its = intensities_list

        # Remove outlier peaks in mz_list
        # self.remove_outlier_peaks_mz()

        # Calculate min and max values for mz and intensity
        self.min_mz = min(self.mzs)
        self.max_mz = max(self.mzs)
        self.min_it = min(self.its)
        self.max_it = max(self.its)

        # Estimate the maximum number of isotopic peaks to consider
        # get the isotopic probabilities of the longest fragment (precursor)
        iso_probs_prec = self.sequence_to_isotopic_probabilities(
            peptidoforms[0].sequence, max_isotope=100
        )

        # from the most likely to the least likely, multiply the intensity of the most intense peak by the isotopic probability until it is below the minimum intensity
        iso_probs_prec = sorted(iso_probs_prec, reverse=True)

        # auto determined max isotopic peaks
        if self.max_isotope == "auto":
            for i in range(1, len(iso_probs_prec)):
                if self.max_it * iso_probs_prec[i] < self.min_it:
                    self.max_isotope = i
                    print("max isotope: ", self.max_isotope)
                    break
        else:
            pass  # keep the max_isotope value from the json file

        # print(iso_probs_prec)

        # Add spectrum peaks to the graph
        self.add_spectrum_peaks_to_graph(self.spectrum_node_radius)

        for peptidoform in peptidoforms:
            # Create the root node corresponding to the precursor peptide
            root_node_id = peptidoform.proforma
            self.add_node(
                root_node_id, node_type="root", color=self.colors["root"], size=15
            )
            # Create intermediate nodes corresponding to possible fragments of the sequence
            self.add_intermediate_0_nodes(root_node_id, peptidoform)

        #self.print_number_each_node_type()
        self.compute_cosine_similarity()
        self.filter_nodes_cosine_similarity(threshold=self.min_cosine_similarity)
        #self.print_number_each_node_type()
        # self.calculate_number_full_overlap_leaf()
        #self.print_number_each_node_type()
        if self.filter_charge:
            self.filter_not_consistent_charge()

        self.filter_fragment_best_isotopic_fit()

        self.propagate_intensity()
        self.remove_non_leaf_branches()
        self.set_node_title_as_attributes()
        self.set_position_nodes()

        print("Graph generated")

    def add_node(self, node_id, parent=None, **kwargs):
        # Copy the default node attributes
        node_attributes = self.default_node_attributes.copy()

        # Update the node attributes with the parent attributes and create an edge
        if parent != None:
            node_attributes.update(self.nodes[parent])
            self.add_edge(parent, node_id)

        # Override the parent attributes with the kwargs
        node_attributes.update(kwargs)

        super().add_node(node_id, **node_attributes)

    def set_node_attributes(self, node_id, **kwargs):
        super().nodes[node_id].update(kwargs)
        self.nodes[node_id].update(kwargs)

    def add_edge(self, parent, child, **kwargs):
        super().add_edge(parent, child, **kwargs)

    def get_frag_code(
        self,
        start_pos,
        end_pos,
        start_ioncap="t",
        end_ioncap="t",
        neutral_loss=[],
        charge=-1,
        nth_isotope=-1,
    ):
        """
        Returns the fragment code of a single fragment.
        """

        # if neutral lost empty
        if len(neutral_loss) == 0:
            NL_code = ""
        else:
            NL_code = "[" + str(neutral_loss) + "]"

        # if charge is not determined:
        if charge == -1:
            charge_code = ""
        else:
            charge_code = "(" + str(charge) + ")"

        # if isotope is not determined:
        if nth_isotope == -1:
            isotope_code = ""
        else:
            isotope_code = "{" + str(nth_isotope) + "}"

        frag_code = (
            str(start_ioncap)
            + ":"
            + str(end_ioncap)
            + "@"
            + str(start_pos)
            + ":"
            + str(end_pos)
            + charge_code
            + NL_code
            + isotope_code
        )

        return frag_code

    def update_frag_code_charge(self, frag_code, charge=None):
        # from a frag_code corresponding to charge 0, return the frag_code corresponding to charged fragments
        if charge is None:
            raise ValueError("Charge must be specified")
        else:
            return frag_code + "(" + str(charge) + ")"

    def update_frag_code_isotope(self, frag_code, isotope=None):
        # from a frag_code corresponding to isotope 0, return the frag_code corresponding to isotopic fragments
        if isotope is None:
            raise ValueError("Isotope must be specified")
        else:
            return frag_code + "{" + str(isotope) + "}"

    def get_fragment_direction(self, start_pos, end_pos, peptide_sequence):
        """
        Returns the fragment direction of a single fragment.
        """
        if start_pos != 1 and end_pos != len(peptide_sequence):
            frag_dir = "I"
        elif start_pos == 1:
            frag_dir = "N"
        elif end_pos == len(peptide_sequence):
            frag_dir = "C"
        else:
            raise ValueError("Fragment direction not found")

        return frag_dir

    def get_fragment_isotope_probabilities(self, node, peptidoform, max_isotope):
        """
        Returns the isotope probabilities of a single fragment.
        """
        # get the node attributes
        node_dict = self.nodes[node]

        # get the isotope probabilities
        sequence = []
        mods = []
        for aa, mod in peptidoform.parsed_sequence[
            node_dict["start_pos"] : node_dict["end_pos"]
        ]:
            sequence.append(aa)
            if not mod is None:
                mods.extend([m.mass for m in mod])
        sequence = "".join(sequence)

        return self.sequence_to_isotopic_probabilities(sequence, max_isotope)

    def sequence_to_isotopic_probabilities(self, sequence, max_isotope):
        composition = mass.Composition(sequence=sequence)
        theoretical_isotopic_cluster = bp.isotopic_variants(
            composition, npeaks=max_isotope, charge=0
        )
        monoisotopic_mass = mass.calculate_mass(composition=composition, charge=0)
        iso_prob = [round(peak.intensity, 4) for peak in theoretical_isotopic_cluster]

        return iso_prob

    def get_fragment_theoretical_mz(
        self,
        peptidoform,
        start,
        end,
        start_ioncap,
        end_ioncap,
        formula,
        charge,
        isotope,
    ):
        # peptide and modification mass
        sequence = []
        mods = []
        for aa, mod in peptidoform.parsed_sequence[start - 1 : end]:
            sequence.append(aa)
            if not mod is None:
                mods.extend([m.mass for m in mod])

        formula = "".join(formula)

        # mass AA sequence
        ps = pyteomics_parser.parse("".join(sequence), show_unmodified_termini=True)
        P = mass.calculate_mass(parsed_sequence=ps)
        # mass modifications
        M = sum(mods)
        # mass start ion cap
        SI = constant.ion_cap_delta_mass[start_ioncap]
        # mass end ion cap
        EI = constant.ion_cap_delta_mass[end_ioncap]
        # hydrogen mass
        H = 1.00784
        # electron mass
        e = 0.00054857990946
        # Proton mass
        p = 1.007276466621
        # loss mass
        # L = mass.calculate_mass(formula, absolute=True)
        L = 0
        # loss mass
        n = 1.008664915

        # Calculate fragment mass
        if charge == 0:
            fragment_mass = P + M + SI + EI + (n * isotope) - L
        else:
            fragment_mass = (
                P + M + SI + EI + (p * (charge)) + (n * isotope) - L
            ) / np.abs(charge)

        fragment_mass = round(float(fragment_mass), 10)
        return fragment_mass

    def update_framgent_theoretical_mz_charge(self, mz, charge=None):
        # from an mz of corresponding to charge 0, calculate the mz corresponding to charged fragments
        if charge is None:
            raise ValueError("Charge must be specified")
        elif charge == 0:
            return mz
        else:
            return round((mz + (1.007276466621 * charge)) / charge, 10)

    def update_fragment_theoretical_mz_isotope(self, mz, charge, isotope=None):
        # from an mz of corresponding to isotope 0, calculate the mz corresponding to isotopic fragments
        if isotope is None:
            raise ValueError("Isotope must be specified")
        elif charge == 0:
            return mz + (1.008664915 * isotope)
        else:
            return mz + ((1.008664915 * isotope) / charge)

    def match_fragment(self, mz_theo):
        idx = np.argmin(np.abs(self.mzs - mz_theo))

        # Check if the closest experimental mass is within the tolerance.
        if self.msms_tol_unit == "ppm":
            if np.abs(self.mzs[idx] - mz_theo) / mz_theo * 1e6 <= self.msms_tol:
                return float(self.its[idx]), float(self.mzs[idx])
            else:
                return 0, 0
        elif self.msms_tol_unit == "da":
            if np.abs(self.mzs[idx] - mz_theo) <= self.msms_tol:
                return float(self.its[idx]), float(self.mzs[idx])
            else:
                return 0, 0
        else:
            raise ValueError("Tolerance unit not recognized")

    # ---------------------------------------------------------------------------- #

    def add_intermediate_0_nodes(self, parent, peptidoform):
        # Create list of tuple of possible fragment positions list[(start, end)], filter min length
        n_fragment_range = [
            (1, i)
            for i in range(len(peptidoform.sequence) - 1)
            if i >= self.min_frag_length
        ]
        c_fragment_range = [
            (i, len(peptidoform.sequence))
            for i in range(2, len(peptidoform.sequence))
            if len(peptidoform.sequence) - i >= self.min_frag_length
        ]
        internal_fragment_range = [
            (i, j)
            for i in range(2, len(peptidoform.sequence) - 1)
            for j in range(i + 1, len(peptidoform.sequence))
            if j - i >= self.min_frag_length
        ]

        # Total iterations for the loading bar
        total_iterations = (
            len(n_fragment_range) + len(c_fragment_range) + len(internal_fragment_range)
        )

        with tqdm(total=total_iterations, desc="Adding Nodes") as pbar:
            # Add N fragments
            for pos in n_fragment_range:
                frag_code = self.get_frag_code(pos[0], pos[1])
                frag_dir = self.get_fragment_direction(
                    pos[0], pos[1], peptidoform.sequence
                )
                frag_mz = int(
                    self.get_fragment_theoretical_mz(
                        peptidoform, pos[0], pos[1], "t", "t", [], 0, 0
                    )
                )  # theoretical mz of the fragment

                node_id = str(f"{pos[0]}:{pos[1]}_{frag_mz}")
                peptidoform_index = self.peptidoforms.index(peptidoform)
                # if already exists, add the peptideform to the list of peptidoforms
                if node_id in self.nodes:
                    self.set_node_attributes(
                        node_id,
                        peptidoforms=self.nodes[node_id]["peptidoforms"]
                        + [peptidoform_index],
                    )
                    self.add_edge(parent, node_id)
                else:
                    self.add_node(
                        node_id,
                        parent=parent,
                        frag_codes=[frag_code],
                        peptidoforms=[peptidoform_index],
                        node_type="intermediate_0",
                        color=self.colors["intermediate_0"],
                        frag_dir=frag_dir,
                        start_pos=pos[0],
                        end_pos=pos[1],
                    )

                    self.add_intermediate_1_nodes(node_id, peptidoform)

                    pbar.update(1)

            # Add C fragments
            for pos in c_fragment_range:
                frag_code = self.get_frag_code(pos[0], pos[1])
                frag_dir = self.get_fragment_direction(
                    pos[0], pos[1], peptidoform.sequence
                )
                frag_mz = int(
                    self.get_fragment_theoretical_mz(
                        peptidoform, pos[0], pos[1], "t", "t", [], 0, 0
                    )
                )  # theoretical mz of the fragment

                node_id = str(f"{pos[0]}:{pos[1]}_{frag_mz}")
                # get peptidoform index
                peptidoform_index = self.peptidoforms.index(peptidoform)

                # if already exists, add the peptideform to the list of peptidoforms and edge
                if node_id in self.nodes:
                    self.set_node_attributes(
                        node_id,
                        peptidoforms=self.nodes[node_id]["peptidoforms"]
                        + [peptidoform_index],
                    )
                    self.add_edge(parent, node_id)

                else:
                    self.add_node(
                        node_id,
                        parent=parent,
                        frag_codes=[frag_code],
                        peptidoforms=[peptidoform_index],
                        node_type="intermediate_0",
                        color=self.colors["intermediate_0"],
                        frag_dir=frag_dir,
                        start_pos=pos[0],
                        end_pos=pos[1],
                    )

                    self.add_intermediate_1_nodes(node_id, peptidoform)

                    pbar.update(1)

            # Add internal fragments (if allowed)
            if not self.terminals_only:
                for pos in internal_fragment_range:
                    frag_code = self.get_frag_code(pos[0], pos[1])
                    frag_dir = self.get_fragment_direction(
                        pos[0], pos[1], peptidoform.sequence
                    )
                    frag_mz = int(
                        self.get_fragment_theoretical_mz(
                            peptidoform, pos[0], pos[1], "t", "t", [], 0, 0
                        )
                    )  # theoretical mz of the fragment

                    node_id = str(f"{pos[0]}:{pos[1]}_{frag_mz}")
                    peptidoform_index = self.peptidoforms.index(peptidoform)
                    # if already exists, add the peptideform to the list of peptidoforms
                    if node_id in self.nodes:
                        self.set_node_attributes(
                            node_id,
                            peptidoforms=self.nodes[node_id]["peptidoforms"]
                            + [peptidoform_index],
                        )
                        self.add_edge(parent, node_id)
                    else:
                        self.add_node(
                            node_id,
                            parent=parent,
                            frag_codes=[frag_code],
                            peptidoforms=[peptidoform_index],
                            node_type="intermediate_0",
                            color=self.colors["intermediate_0"],
                            frag_dir=frag_dir,
                            start_pos=pos[0],
                            end_pos=pos[1],
                        )

                        self.add_intermediate_1_nodes(node_id, peptidoform)

                        pbar.update(1)

    # ---------------------------------------------------------------------------- #

    def add_intermediate_1_nodes(self, parent, peptidoform):
        frag_dir = self.nodes[parent]["frag_dir"]

        if frag_dir == "N":
            ioncaps_types = [("t", end_ioncap) for end_ioncap in self.end_ioncaps_types]
        elif frag_dir == "C":
            ioncaps_types = [
                (start_ioncap, "t") for start_ioncap in self.start_ioncaps_types
            ]
        elif frag_dir == "I":
            ioncaps_types = list(
                itertools.product(self.start_ioncaps_types, self.end_ioncaps_types)
            )

        for ioncaps in ioncaps_types:
            frag_mz = self.get_fragment_theoretical_mz(
                peptidoform,
                self.nodes[parent]["start_pos"],
                self.nodes[parent]["end_pos"],
                ioncaps[0],
                ioncaps[1],
                [],
                0,
                0,
            )

            frag_code = self.get_frag_code(
                self.nodes[parent]["start_pos"],
                self.nodes[parent]["end_pos"],
                ioncaps[0],
                ioncaps[1],
            )

            # Is their an already existing intermediate 1 node with an mz within msms_tol:

            # Get the intermediate 1 nodes with the closest mz
            closest_mz = min(self.I_1_nodes_mz, key=lambda x: abs(x - frag_mz))

            # if the mz within the tolerance
            is_within_tol = False
            if self.msms_tol_unit == "ppm":
                # print("ppm error: ", np.abs(closest_mz - frag_mz) / frag_mz * 1e6)
                if np.abs(closest_mz - frag_mz) / frag_mz * 1e6 <= self.msms_tol:
                    is_within_tol = True
            elif self.msms_tol_unit == "da":
                # print("Da error: ", np.abs(closest_mz - frag_mz))
                if np.abs(closest_mz - frag_mz) <= self.msms_tol:
                    is_within_tol = True

            # is whiting the tolerance, add the frag_code and an edge to the existing intermediate 1 node
            if is_within_tol:
                # if internal over terminal is not allowed, check if the intermediate 1 node is terminal
                if self.internal_over_terminal == False:
                    if (
                        frag_dir == "I"
                        and self.nodes[str(closest_mz)]["frag_dir"] != "I"
                    ):
                        continue
                else:
                    print(
                        f"node {frag_code} (mz: {frag_mz}) is within the tolerance of node {closest_mz} (mz: {closest_mz})"
                    )
                    self.set_node_attributes(
                        str(closest_mz),
                        frag_codes=self.nodes[str(closest_mz)]["frag_codes"]
                        + [frag_code],
                    )
                    self.add_edge(parent, str(closest_mz))
            # if not within the tolerance, create a new intermediate 1 node
            else:
                self.add_node(
                    str(frag_mz),
                    parent=parent,
                    frag_codes=[frag_code],
                    node_type="intermediate_1",
                    color=self.colors["intermediate_1_" + frag_dir],
                    start_pos=self.nodes[parent]["start_pos"],
                    end_pos=self.nodes[parent]["end_pos"],
                    frag_dir=frag_dir,
                    mz=frag_mz,
                )

                # add mz to the list of I_1 nodes
                self.I_1_nodes_mz.append(frag_mz)

                self.add_intermediate_2_nodes(str(frag_mz), peptidoform)

    # ---------------------------------------------------------------------------- #

    def add_intermediate_2_nodes(self, parent, peptidoform, charge_prob=None):
        # if no charge prob, use default charge lookup
        if charge_prob is None and self.max_charge == 0:
            charge_lookup = [0]
        elif charge_prob is None:
            charge_lookup = list(range(1, self.max_charge + 1))
        else:
            length = self.nodes[parent]["end_pos"] - self.nodes[parent]["start_pos"]
            charge_lookup = self.get_charge_lookup(charge_prob, length)

        # annotate the intermediate_1 node with the charge lookup
        node_ids = []
        for charge in charge_lookup:
            frag_mz = self.update_framgent_theoretical_mz_charge(
                self.nodes[parent]["mz"], charge
            )
            frag_codes = [
                self.update_frag_code_charge(frag_code, charge)
                for frag_code in self.nodes[parent]["frag_codes"]
            ]

            self.add_node(
                str(frag_mz) + "_2",
                frag_codes=frag_codes,
                parent=parent,
                node_type="intermediate_2",
                color=self.colors["intermediate_2"],
                charge=charge,
                mz=frag_mz,
            )

            # print("-----ADDING INTERMEDIATE 2 NODE: ", str(frag_mz), "WITH CHARGE: ", charge)

            self.add_leaf_nodes_and_match(str(frag_mz) + "_2", peptidoform)

    def get_charge_lookup(self, charge_prob, length):
        if charge_prob is not None and length in charge_prob:
            charge_prob = charge_prob[length]
            charge_lookup = sorted(charge_prob, key=charge_prob.get, reverse=True)
        elif self.max_charge == 0:
            charge_lookup = [0]
        else:
            charge_lookup = range(1, self.max_charge + 1)
        return charge_lookup

    # ---------------------------------------------------------------------------- #

    def add_leaf_nodes_and_match(self, parent, peptidoform, no_overlaping_match=False):
        """
        this method compute the probability of the isotopic pattern of a fragment
        and generate the leaf node from the most likely isotopic peak to the less likely,
        and match the fragment to the peak list
        """
        if self.monoisotopic:
            self.max_isotope = 1
            iso_probs = [1]
            iso_probs_order = [0]

        else:
            iso_probs = self.get_fragment_isotope_probabilities(
                parent, peptidoform, max_isotope=self.max_isotope
            )
            iso_probs_order = np.argsort(iso_probs)[::-1]
            iso_probs_order = iso_probs_order[: self.max_isotope]

        first_isotope = True

        for i in iso_probs_order:
            # Matching the fragment to the peak list
            frag_leaf_mz = self.update_fragment_theoretical_mz_isotope(
                self.nodes[parent]["mz"], self.nodes[parent]["charge"], i
            )
            frag_leaf_id = str(frag_leaf_mz) + "_i"
            frag_codes_leaf = [
                self.update_frag_code_isotope(frag_code, i)
                for frag_code in self.nodes[parent]["frag_codes"]
            ]

            # if theoretical mz in mz range of spectrum range :
            if frag_leaf_mz > self.min_mz and frag_leaf_mz < self.max_mz:
                #print(self.match_fragment(frag_leaf_mz))
                its_match, mz_match = self.match_fragment(frag_leaf_mz)

                # add the leaf node to the graph
                if its_match > 0:
                    if first_isotope:
                        self.add_node(
                            frag_leaf_id,
                            parent=parent,
                            frag_codes=frag_codes_leaf,
                            node_type="leaf",
                            color=self.colors["leaf_first_isotope"],
                            intensity=its_match,
                            mz=frag_leaf_mz,
                            isotope=int(i),
                            isotope_prob=iso_probs[i],
                        )

                        # print("----------matched fragment", frag_codes_leaf, "intensity", its_match, "isotope", i)

                        first_isotope_intensity = its_match
                        first_isotope_probability = iso_probs[i]

                    else:
                        self.add_node(
                            frag_leaf_id,
                            parent=parent,
                            node_type="leaf",
                            color=self.colors["leaf_other_isotopes"],
                            intensity=its_match,
                            mz=frag_leaf_mz,
                            isotope=int(i),
                            isotope_prob=iso_probs[i],
                        )

                        # print("----------matched fragment", frag_codes_leaf, "intensity", its_match, "isotope", i)

                    self.add_edge(frag_leaf_id, str(mz_match) + "_0")

                elif its_match == 0 and not first_isotope:
                    self.add_node(
                        str(frag_leaf_mz) + "_i",
                        parent=parent,
                        node_type="leaf",
                        color=self.colors["leaf_not_matched"],
                        intensity=its_match,
                        mz=frag_leaf_mz,
                        isotope=int(i),
                        isotope_prob=iso_probs[i],
                    )

                    # add fake peak node of half the minimal intensity
                    fake_its = self.min_it / self.fake_peak_factor

                    x0, y0 = self.point_on_circle_mz(
                        frag_leaf_mz, radius=self.spectrum_node_radius
                    )
                    x1, y1 = self.point_on_circle_mz(
                        frag_leaf_mz,
                        radius=self.spectrum_node_radius
                        + self.get_peak_height(fake_its),
                    )

                    self.add_node(
                        str(frag_leaf_mz) + "_0",
                        parent=str(frag_leaf_mz) + "_i",
                        node_type="peak",
                        color=self.colors["peak_fake"],
                        mz=frag_leaf_mz,
                        intensity=self.min_it / self.fake_peak_factor,
                        physics=False,
                        x=x0,
                        y=y0,
                        size=2,
                    )

                    self.add_node(
                        str(frag_leaf_mz) + "_1",
                        parent=str(frag_leaf_mz) + "_0",
                        node_type="peak",
                        color=self.colors["peak_fake"],
                        mz=frag_leaf_mz,
                        intensity=self.min_it / 2,
                        physics=False,
                        x=x1,
                        y=y1,
                        size=2,
                    )

                    self.add_edge(str(frag_leaf_mz) + "_i", str(frag_leaf_mz) + "_0")

                elif its_match == 0 and first_isotope:
                    break

                # add line to matching_success dataframe with concat:

                # self.matching_success = pd.concat([self.matching_success, pd.DataFrame([[frag_code_leaf,its_match]],columns=['fragment_code','intensity'])])

                # if last isotopic peak in list
                if len(iso_probs) - 1 < i + 1:
                    break

                first_isotope = False

                # stop the matching if the theoretical intensity of the next isotope is too low
                try:
                    next_theo_i = first_isotope_intensity * (
                        iso_probs[i + 1] / first_isotope_probability
                    )
                    if next_theo_i < (self.min_it * self.min_it_threshold_factor):
                        break
                except ZeroDivisionError:
                    break

    def remove_outlier_peaks_mz(self, threshold=3):
        """remove outlyig peaks based on their m/z zscore"""
        mean = np.mean(self.mzs)
        std = np.std(self.mzs)
        zscores = [(mz - mean) / std for mz in self.mzs]
        # get index of outliers
        outliers = np.where(np.abs(zscores) > threshold)[0]
        # remove outliers from mz and its
        self.mzs = np.delete(self.mzs, outliers)
        self.its = np.delete(self.its, outliers)

    def point_on_circle_mz(self, mz, radius):
        if mz > self.max_mz or mz < self.min_mz:
            raise ValueError("mz value out of range")

        angle_range = 360  # 250 - 10
        angle_per_unit = angle_range / (self.max_mz - self.min_mz)

        # calculate angle based on mz value
        angle = (mz - self.min_mz) * angle_per_unit

        # calculate coordinates based on angle and radius
        angle_radians = math.radians(angle)
        x = radius * math.cos(angle_radians)
        y = radius * math.sin(angle_radians)

        return x, y

    def point_on_circle_mz_intermediate_node(self, mz, radius, min_mz, max_mz):
        pass

    def point_on_circle_length(self, length, radius):
        angle_range = 360
        angle_per_unit = angle_range / len(self.peptidoforms[0].sequence)

        # calculate angle based on length of the fragment
        angle = length * angle_per_unit

        # calculate coordinates based on angle and radius
        angle_radians = math.radians(angle)
        x = radius * math.cos(angle_radians)
        y = radius * math.sin(angle_radians)
        return round(x, 1), round(y, 1)

    def point_on_circle_length_charge(self, length_charge, radius):
        angle_range = 360
        angle_per_unit = angle_range / (
            len(self.peptidoform.sequence) * self.max_charge + self.max_charge
        )

        # calculate angle based on length of the fragment
        angle = length_charge * angle_per_unit

        # calculate coordinates based on angle and radius
        angle_radians = math.radians(angle)
        x = radius * math.cos(angle_radians)
        y = radius * math.sin(angle_radians)
        return round(x, 1), round(y, 1)

    def get_peak_height(self, intensity):
        if intensity == 0:
            return 0
        else:
            # determine max intensity as 2STD above mean
            upper_limit = np.mean(self.its) + (1 * np.std(self.its))

            return (intensity / upper_limit) * self.peak_scaling_factor + 10

    def set_position_nodes(self):

        with tqdm(total=len(list(self.nodes.keys())), desc="Setting node positions (Vis)") as pbar:
            for node in self.nodes:
                if self.nodes[node]["node_type"] == "leaf":
                    x, y = self.point_on_circle_mz(
                        self.nodes[node]["mz"], self.circle_radius_leaf
                    )
                    self.set_node_attributes(node, x=x, y=y, physics=False)
                if self.nodes[node]["node_type"] == "root":
                    self.set_node_attributes(node, x=0, y=0, physics=False)
                if self.nodes[node]["node_type"] == "intermediate_0":
                    if self.nodes[node]["frag_dir"] == "I":
                        x, y = self.point_on_circle_length(
                            self.nodes[node]["end_pos"] - self.nodes[node]["start_pos"],
                            self.circle_radius_intermediate_0_I,
                        )
                    elif self.nodes[node]["frag_dir"] == "N":
                        x, y = self.point_on_circle_length(
                            self.nodes[node]["end_pos"] - self.nodes[node]["start_pos"],
                            self.circle_radius_intermediate_0_N,
                        )
                    elif self.nodes[node]["frag_dir"] == "C":
                        x, y = self.point_on_circle_length(
                            self.nodes[node]["end_pos"] - self.nodes[node]["start_pos"],
                            self.circle_radius_intermediate_0_C,
                        )

                    self.set_node_attributes(node, x=x, y=y, physics=False)

                pbar.update(1)

    def add_spectrum_peaks_to_graph(self, radius):
        for mz, intensity in zip(self.mzs, self.its):
            # Define the coordinates of the two node
            x0, y0 = self.point_on_circle_mz(mz, radius)
            x1, y1 = self.point_on_circle_mz(
                mz, radius + self.get_peak_height(intensity)
            )

            # add the two node with unique identifiers
            self.add_node(
                str(mz) + "_0",
                parent=None,
                size=2,
                node_type="peak",
                x=x0,
                y=y0,
                physics=False,
                type="peak",
                mz=mz,
                intensity=intensity,
                color=self.colors["peak"],
            )
            self.add_node(
                str(mz) + "_1",
                parent=str(mz) + "_0",
                size=2,
                node_type="peak_viz",
                x=x1,
                y=y1,
                physics=False,
                color=self.colors["peak"],
            )

    def clean_graph(self, node_to_remove=[]):
        visited = {}
        for node in self.nodes:
            if self.nodes[node]["node_type"] == "leaf":
                visited[node] = True

        for node in visited.keys():
            visited.update(self.dfs(node, visited))

        for node in visited.keys():
            self.remove_node(node)

    def remove_non_leaf_branches(self):
        visited = {}  # dictionary to keep track of visited nodes

        # Step 1: Identify all nodes of type "leaf" and mark them as visited
        for node in self.nodes:
            if self.nodes[node]["node_type"] == "leaf":
                visited[node] = True

        # Step 2: Perform DFS from visited nodes, marking all visited nodes as visited
        def dfs(node):
            visited[node] = True  # mark node as visited
            for neighbor in self.predecessors(node):
                if (
                    neighbor not in visited
                    and self.nodes[neighbor]["node_type"] != "leaf"
                ):
                    dfs(neighbor)

        for node in list(visited.keys()):
            dfs(node)

        # Step 3: Remove all unmarked nodes and their branches except node of type intermediate_0
        with tqdm(total=len(list(self.nodes.keys())), desc="removing non annotated fragments") as pbar:
            for node in list(self.nodes.keys()):
                if node not in visited and self.nodes[node]["node_type"] not in [
                    "root",
                    "peak",
                    "peak_viz",
                ]:
                    self.remove_node(node)
                pbar.update(1)

    def get_overlapping_leaf_nodes(self, leaf_node_code):
        # from a leaf node code get all the nodes linked to the same peak node

        overlapping_leaf_nodes = []

        # get peak nodes
        peak_nodes = self.successors(leaf_node_code)

        for peak_node in peak_nodes:
            for node in self.predecessors(peak_node):
                if self.nodes[node]["node_type"] == "leaf":
                    overlapping_leaf_nodes.append(node)

        return overlapping_leaf_nodes

    def remove_overlapping_internal_nodes(self):
        # remove all nodes corresponding to internal ions that are overlapping with other terminal ions (C, N)
        to_remove = []
        for node in self.nodes:
            if self.nodes[node]["node_type"] == "leaf":
                if self.nodes[node]["frag_dir"] == "I":
                    # get the types of the overlapping leaf nodes
                    for overlapping_leaf_node in self.get_overlapping_leaf_nodes(node):
                        if self.nodes[overlapping_leaf_node]["frag_dir"] in ["C", "N"]:
                            to_remove.append(node)
                            break
        # print(len(to_remove), " Internal nodes have been removed")
        for node in to_remove:
            self.remove_node(node)

    def set_node_title_as_attributes(self):
        for node in self.nodes:
            str_attr = ""
            for key, value in self.nodes[node].items():
                str_attr += key + "=" + str(value) + "\n"

            self.set_node_attributes(node, title=str_attr)

    def node_size_from_intensity(self, intensity):
        if intensity == 0:
            return 0
        else:
            return (intensity / max(self.its)) * 30 + 3

    def propagate_intensity(self):
        # propagate intensities from peaks to I0

        # for each peak node
        with tqdm(total=len(self.nodes), desc="Propagating Intensity (peak nodes)") as pbar:
            for node in self.nodes:
                parents = list(self.predecessors(node))
                if len(parents) != 0:
                    if (
                        self.nodes[node]["node_type"] == "peak"
                        and self.nodes[parents[0]]["node_type"] == "leaf"
                    ):
                        I2_weights = []
                        I2_weights = [
                            self.nodes[list(self.predecessors(parent))[0]]["weight"]
                            for parent in parents
                        ]

                        # normalize weights to sum to 1 and set leaf node weight
                        I2_weights = [i / sum(I2_weights) for i in I2_weights]
                        for leaf_node, weight in zip(parents, I2_weights):
                            self.set_node_attributes(leaf_node, weight=weight)
                pbar.update(1)

        # set the intensity size of the leaf nodes
        with tqdm(total=len(self.nodes), desc="Propagating Intensity (leaf nodes)") as pbar:
            for node in self.nodes:
                if self.nodes[node]["node_type"] == "leaf":
                    # set leaf intensity from connected peak node
                    if len(list(self.successors(node))) > 0:
                        intensity = (
                            self.nodes[list(self.successors(node))[0]]["intensity"]
                            * self.nodes[node]["weight"]
                        )
                        self.set_node_attributes(
                            node,
                            size=self.node_size_from_intensity(intensity),
                            intensity=intensity,
                        )
                    else:
                        self.set_node_attributes(
                            node, size=self.node_size_from_intensity(0), intensity=0
                        )
                pbar.update(1)

        # propagate the intensity from leaf nodes to intermediate_2 nodes
        with tqdm(total=len(self.nodes), desc="Propagating Intensity (I2 nodes)") as pbar:
            for node in self.nodes:
                if self.nodes[node]["node_type"] == "leaf":
                    intensity = self.nodes[node]["intensity"]
                    # propagate the intensity to intermediate_2 nodes
                    for parent in self.predecessors(node):
                        intensity = self.nodes[parent]["intensity"] + intensity
                        # add intensity to intermediate_2 node
                        self.set_node_attributes(
                            parent,
                            size=self.node_size_from_intensity(intensity),
                            intensity=intensity,
                        )
                pbar.update(1)

        # propagate the intensity from intermediate_2 nodes to intermediate_1 nodes
        with tqdm(total=len(self.nodes), desc="Propagating Intensity (I1 nodes)") as pbar:
            for node in self.nodes:
                if self.nodes[node]["node_type"] == "intermediate_2":
                    intensity = self.nodes[node]["intensity"]
                    # propagate the intensity to intermediate_1 nodes
                    for parent in self.predecessors(node):
                        intensity = self.nodes[parent]["intensity"] + intensity
                        # add intensity to intermediate_1 node
                        self.set_node_attributes(
                            parent,
                            size=self.node_size_from_intensity(intensity),
                            intensity=intensity,
                        )
                pbar.update(1)

        # propagate the intensity from intermediate_1 nodes to intermediate_0 nodes
        with tqdm(total=len(self.nodes), desc="Propagating Intensity (I0 nodes)") as pbar:
            for node in self.nodes:
                if self.nodes[node]["node_type"] == "intermediate_1":
                    intensity = self.nodes[node]["intensity"]
                    # propagate the intensity to intermediate_0 nodes

                    # get the number of parents
                    parents = list(self.predecessors(node))
                    if len(parents) != 0:
                        for parent in parents:
                            # add intensity to intermediate_0 node
                            intensity = self.nodes[parent]["intensity"] + (
                                intensity / len(parents)
                            )
                            # add intensity to intermediate_0 node
                            self.set_node_attributes(
                                parent,
                                size=self.node_size_from_intensity(intensity),
                                intensity=intensity,
                            )
                pbar.update(1)



    def compute_cosine_similarity(self):
        for node in self.nodes:
            # print("node", node)
            # print(self.nodes[node])
            if self.nodes[node]["node_type"] == "intermediate_2":
                iso_probs = []
                iso_inten = []

                for child in self.successors(node):
                    iso_probs.append(self.nodes[child]["isotope_prob"])
                    iso_inten.append(self.nodes[child]["intensity"])

                if sum(iso_inten) != 0 and sum(iso_probs) != 0:
                    # zero values are replaced with half the detection limit
                    # iso_inten = [self.min_it / 2 if i == 0 else i for i in iso_inten]

                    # convert intensity to probability
                    iso_inten = [i / sum(iso_inten) for i in iso_inten]
                    # sum iso_probs to 1
                    iso_probs = [i / sum(iso_probs) for i in iso_probs]

                    # compute KL divergence
                    kl_val = sum(rel_entr(iso_probs, iso_inten))
                    # compute cosine similarity
                    cs_val = np.dot(iso_probs, iso_inten) / (
                        norm(iso_probs) * norm(iso_inten)
                    )

                    # add val to node attributes
                    self.set_node_attributes(
                        node, kl_divergence=kl_val, cosine_similarity=cs_val
                    )

                else:
                    self.set_node_attributes(
                        node, kl_divergence=-1, cosine_similarity=-1
                    )

    def get_cs_scores(self, iondir):
        cs_scores = []

        for node in self.nodes:
            if (
                self.nodes[node]["node_type"] == "intermediate_2"
                and self.nodes[node]["frag_dir"] in iondir
            ):
                cs_scores.append(self.nodes[node]["cosine_similarity"])

        return cs_scores

    def get_intermediate_2_intensities(self, iondir):
        intensities = []
        for node in self.nodes:
            if (
                self.nodes[node]["node_type"] == "intermediate_2"
                and self.nodes[node]["frag_dir"] in iondir
            ):
                intensities.append(self.nodes[node]["intensity"])

        return intensities

    def get_fragmentation_coverage_terminal(self):
        # return the number of position in the sequence for which fragment are found:
        found_positions_n_term = {
            i: 0 for i in range(len(self.peptidoform.sequence) + 1)
        }
        found_positions_c_term = {
            i: 0 for i in range(len(self.peptidoform.sequence) + 1)
        }

        for node in self.nodes:
            if self.nodes[node]["node_type"] == "intermediate_2":
                if self.nodes[node]["frag_dir"] == "N":
                    found_positions_n_term[self.nodes[node]["end_pos"]] += 1
                elif self.nodes[node]["frag_dir"] == "C":
                    found_positions_c_term[self.nodes[node]["start_pos"]] += 1

        # get percentage of position with at least one fragment
        coverage_n_term = sum(
            [1 for i in found_positions_n_term.values() if i > 0]
        ) / len(found_positions_n_term)
        coverage_c_term = sum(
            [1 for i in found_positions_c_term.values() if i > 0]
        ) / len(found_positions_c_term)

        return coverage_n_term, coverage_c_term

    def get_ion_series_mz(self, peptidoform, ion_cap_type):
        # return the mz of the ion series
        codes = []
        mzs = []

        for p in range(1, len(peptidoform.sequence)):
            for c in range(1, self.max_charge):
                # check ion series direction:
                if constant.ion_direction[ion_cap_type] == "n-term":
                    # n-term ion series
                    mzs.append(
                        self.get_fragment_theoretical_mz(
                            peptidoform=peptidoform,
                            start=0,
                            end=p,
                            start_ioncap="t",
                            end_ioncap=ion_cap_type,
                            formula=[],
                            charge=c,
                            isotope=0,
                        )
                    )
                    codes.append(
                        self.get_frag_code(
                            0, p, start_ioncap="t", end_ioncap=ion_cap_type, charge=c
                        )
                    )
                elif constant.ion_direction[ion_cap_type] == "c-term":
                    # c-term ion series
                    mzs.append(
                        self.get_fragment_theoretical_mz(
                            peptidoform=peptidoform,
                            start=p,
                            end=len(peptidoform.sequence),
                            start_ioncap=ion_cap_type,
                            end_ioncap="t",
                            formula=[],
                            charge=c,
                            isotope=0,
                        )
                    )
                    codes.append(
                        self.get_frag_code(
                            p,
                            len(peptidoform.sequence),
                            start_ioncap=ion_cap_type,
                            end_ioncap="t",
                            charge=c,
                        )
                    )
        # make pandas dataframe
        df = pd.DataFrame({"mz": mzs, "code": codes})

        return df

    def filter_nodes_cosine_similarity(self, threshold):
        """filter intermediate_2 and corresponding leaf nodes based on cosine similarity threshold"""

        to_remove = []
        for node in self.nodes:
            if self.nodes[node]["node_type"] == "intermediate_2":
                if self.nodes[node]["cosine_similarity"] < threshold:
                    # remove intermediate_2 node and corresponding leaf nodes
                    for child in self.successors(node):
                        to_remove.append(child)

                        # remove peak node if the leaf node is the only parent
                        for child_peak in self.successors(child):
                            if len(list(self.predecessors(child_peak))) == 1:
                                to_remove.append(child_peak)
                                # remove the child peak node
                                for child_peak_peak in self.successors(child_peak):
                                    to_remove.append(child_peak_peak)

                    to_remove.append(node)

        # print(len(to_remove), " nodes have been removed based on cosine similarity threshold")

        for node in to_remove:
            if self.has_node(node):
                self.remove_node(node)

    def filter_not_consistent_charge(self, min_occurence=3):
        """This function genrate the set of accepted charge in function of length from the
        terminal framgent and then filter internal framgent if they do not have a consistent charge
        """

        # first extract the charges and length of the terminal fragments
        charges = []
        lengths = []

        for node in self.nodes:
            if self.nodes[node]["node_type"] == "intermediate_2":
                if self.nodes[node]["frag_dir"] in ["N", "C"]:
                    charges.append(self.nodes[node]["charge"])
                    lengths.append(
                        self.nodes[node]["end_pos"] - self.nodes[node]["start_pos"]
                    )

        # #gather in a dataframe
        df = pd.DataFrame({"charge": charges, "length": lengths})

        # plot (plotly)
        fig = px.scatter(df, x="length", y="charge", color="charge")
        fig.show()

        # for each length define a range of accepted charges (90% probability)
        accepted_charge = {}
        for length in set(lengths):
            # get the charges for this length
            charges = df[df["length"] == length]["charge"]
            # get the 90% probability interval
            interval = np.percentile(charges, [5, 95])
            accepted_charge[length] = interval

        # print(accepted_charge)

        # if missing length use average from neighbouring available lengths
        for length in range(1, len(self.peptidoforms[0].sequence)):
            if length not in accepted_charge.keys():
                # get the closest length
                closest_length = min(
                    accepted_charge.keys(), key=lambda x: abs(x - length)
                )
                accepted_charge[length] = accepted_charge[closest_length]

        # print(accepted_charge)

        to_remove = []
        for node in self.nodes:
            if self.nodes[node]["node_type"] == "intermediate_2":
                if self.nodes[node]["frag_dir"] == "I":
                    if (
                        self.nodes[node]["charge"]
                        < accepted_charge[
                            self.nodes[node]["end_pos"] - self.nodes[node]["start_pos"]
                        ][0]
                        or self.nodes[node]["charge"]
                        > accepted_charge[
                            self.nodes[node]["end_pos"] - self.nodes[node]["start_pos"]
                        ][1]
                    ):
                        # remove intermediate_2 node and corresponding leaf nodes
                        for child in self.successors(node):
                            to_remove.append(child)

                            # remove peak node if the leaf node is the only parent
                            for child_peak in self.successors(child):
                                if len(list(self.predecessors(child_peak))) == 1:
                                    to_remove.append(child_peak)
                                    # remove the child peak node
                                    for child_peak_peak in self.successors(child_peak):
                                        to_remove.append(child_peak_peak)

        # print(len(to_remove), " nodes have been removed based on charge consistency")

        for node in to_remove:
            if self.has_node(node):
                self.remove_node(node)

    def negative_cosine_similarity(self, weights, target, candidates):
        weighted_sum = np.sum(weights[:, np.newaxis] * candidates, axis=0)
        return -np.dot(target, weighted_sum) / (
            np.linalg.norm(target) * np.linalg.norm(weighted_sum)
        )

    def find_best_combination(self, target, candidates):
        # standardize target distribution
        target = target / np.sum(target)
        # convert candidates to list of np arrays
        candidates = [np.array(c) for c in candidates]

        num_candidates = len(candidates)
        initial_weights = (
            np.ones(num_candidates) / num_candidates
        )  # Initial weights are equally distributed
        constraints = {
            "type": "eq",
            "fun": lambda w: np.sum(w) - 1,
        }  # Constraint: weights sum to 1
        bounds = [
            (0, 1) for _ in range(num_candidates)
        ]  # Weights must be between 0 and 1

        # Use the minimize function from SciPy to find the optimal weights
        result = minimize(
            self.negative_cosine_similarity,
            initial_weights,
            args=(target, candidates),
            method="SLSQP",
            bounds=bounds,
            constraints=constraints,
        )
        optimal_weights = result.x

        # Apply the minimum weight threshold
        # set weights to 0 if below threshold

        optimal_weights = [
            0 if w < self.min_weight_threshold else w for w in optimal_weights
        ]
        # standardize weights to sum up to 1
        optimal_weights = optimal_weights / np.sum(optimal_weights)

        # if only one value in target distributions, weights as 1/number candidates
        if len(target) == 1:
            optimal_weights = np.ones(num_candidates) / num_candidates

        if round(sum(optimal_weights), 1) != 1:
            # print(sum(optimal_weights))
            # print(optimal_weights)

            print("could not find optimal weights, using 0 weights")
            # return 0 optimal weights if sum is not 1
            optimal_weights = np.zeros(num_candidates)

            # raise ValueError("weights do not sum to 1")

        # Calculate the best combination of candidate distributions using the optimal weights
        best_combination = np.sum(optimal_weights[:, np.newaxis] * candidates, axis=0)

        # print resuls:

        # print("target distribution:", target)
        # print("candidate distributions:", candidates)
        # print("optimal weights:", optimal_weights)
        # print("best combination:", best_combination)

        return optimal_weights

    def get_overlaping_fragment_groups(self):
        """returns groups of intermediate_2 nodes whose leaf are connected to the same peak node"""

        # create a copy of the graph by creating a new networkx object:

        graph = nx.Graph()
        graph.add_nodes_from(self.nodes)
        graph.add_edges_from(self.edges)

        # remove all nodes that are not intermediate_2 or leaf or peak nodes
        for node in self.nodes:
            if self.nodes[node]["node_type"] not in ["intermediate_2", "leaf", "peak"]:
                graph.remove_node(node)

        # get all connected components
        groups = list(nx.connected_components(graph))

        # # keep only intermediate_2 nodes in connected components
        # for g in groups:
        #     to_remove = []
        #     for node in g:
        #         if self.nodes[node]["node_type"] != "intermediate_2":
        #             to_remove.append(node)
        #     for node in to_remove:
        #         g.remove(node)

        #    g = list(g)

        return groups

    def filter_fragment_best_isotopic_fit(self):
        """find the best isotopic fit in each group of overlapping fragments and remove the other ones"""

        groups = self.get_overlaping_fragment_groups()
        total_removed = 0

        # with loading bar
        with tqdm(total=len(groups), desc="Filtering Nodes") as pbar:
            for g in groups:
                if len(g) > 1:
                    # get the target distribution from the peak node
                    peak_nodes = []
                    peak_nodes_mzs = []
                    peak_nodes_its = []

                    for node in g:
                        if self.nodes[node]["node_type"] == "peak":
                            peak_nodes.append(node)
                            peak_nodes_mzs.append(float(self.nodes[node]["mz"]))
                            peak_nodes_its.append(self.nodes[node]["intensity"])

                    if len(peak_nodes) == 0:
                        continue

                    # order the three lists by increasing values in peak_nodes_mzs
                    peak_nodes_mzs, peak_nodes, peak_nodes_its = zip(
                        *sorted(zip(peak_nodes_mzs, peak_nodes, peak_nodes_its))
                    )
                    # convert intensities value to distribution
                    peak_nodes_probs = [i / sum(peak_nodes_its) for i in peak_nodes_its]

                    # print("sorted peak node mzs:", peak_nodes_mzs)
                    # print("target distribution:", peak_nodes_probs)

                    # get the candidate distributions from the intermediate_2 nodes
                    candidate_distributions = []
                    candidate_distributions_parent_nodes = []
                    for node in g:
                        if self.nodes[node]["node_type"] == "intermediate_2":
                            candidate_distributions_parent_nodes.append(node)
                            # get peak node connected to the intermediate_2 node (iterate over leaf and get peaks)
                            connected_peak_nodes = []
                            connected_leaf_nodes = []
                            for leaf_node in self.successors(node):
                                connected_leaf_nodes.append(leaf_node)
                                for peak_node in self.successors(leaf_node):
                                    if peak_node not in connected_peak_nodes:
                                        connected_peak_nodes.append(peak_node)

                            # for each peak nodes get the theoretical probability of the fragment if there is one else appedn 0:
                            candidate_distribution = [0 for i in range(len(peak_nodes))]
                            i = 0
                            for peak_node in peak_nodes:
                                if peak_node in connected_peak_nodes:
                                    # find the connected leaf node
                                    for leaf_from_peak in self.predecessors(peak_node):
                                        if leaf_from_peak in connected_leaf_nodes:
                                            candidate_distribution[i] = self.nodes[
                                                leaf_from_peak
                                            ]["isotope_prob"]
                                i += 1
                            # standardize the distribution to sum to 1
                            if sum(candidate_distribution) != 0:
                                candidate_distribution = [
                                    i / sum(candidate_distribution)
                                    for i in candidate_distribution
                                ]
                            candidate_distributions.append(candidate_distribution)

                    # for i in range(len(candidate_distributions)):
                    #     print(
                    #         candidate_distributions_parent_nodes[i],
                    #         "---->",
                    #         candidate_distributions[i],
                    #     )

                    # find the best combination of candidate distributions
                    # print(len(g), " overlapping fragments")

                    if len(g) < self.max_overlap_group_size:
                        best_weights = self.find_best_combination(
                            peak_nodes_probs, candidate_distributions
                        )

                        # print("best combination (by weights):", best_weights)

                        # remove nodes that are where the weight is 0
                        # add calculated weight to intermediate_2 nodes
                        to_remove = []
                        for i in range(len(best_weights)):
                            if best_weights[i] == 0:
                                to_remove.append(
                                    candidate_distributions_parent_nodes[i]
                                )

                                # remove child nodes (leaf nodes)
                                # print("parent node distribution:", candidate_distributions_parent_nodes[i])
                                for leaf_node in self.successors(
                                    candidate_distributions_parent_nodes[i]
                                ):
                                    # print("removing leaf node:", leaf_node)
                                    to_remove.append(leaf_node)
                                    # remove peak nodes if only one leaf node is connected to it
                                    for peak_node in self.successors(leaf_node):
                                        if len(list(self.predecessors(peak_node))) == 1:
                                            to_remove.append(peak_node)
                                            # remove succesor nodes of peak nodes
                                            to_remove.append(
                                                list(self.successors(peak_node))[0]
                                            )
                            if best_weights[i] != 0:
                                self.set_node_attributes(
                                    candidate_distributions_parent_nodes[i],
                                    weight=best_weights[i],
                                )
                        # print(len(to_remove), " nodes have been removed based on best isotopic fit")

                        # remove the nodes from the graph
                        for node in to_remove:
                            self.remove_node(node)
                            total_removed += 1

                    else:  # the group is too big, use the average of the candidate distributions
                        print(
                            "length of group is too large, removing internal fragments"
                        )

                        # remove candidate distribution and node for internal fragments
                        to_remove = []
                        for i in range(len(candidate_distributions)):
                            if (
                                self.nodes[candidate_distributions_parent_nodes[i]][
                                    "frag_dir"
                                ]
                                == "I"
                            ):
                                to_remove.append(
                                    candidate_distributions_parent_nodes[i]
                                )

                                # remove child nodes (leaf nodes)
                                # print("parent node distribution:", candidate_distributions_parent_nodes[i])
                                for leaf_node in self.successors(
                                    candidate_distributions_parent_nodes[i]
                                ):
                                    # print("removing leaf node:", leaf_node)
                                    to_remove.append(leaf_node)
                                    # remove peak nodes if only one leaf node is connected to it
                                    for peak_node in self.successors(leaf_node):
                                        if len(list(self.predecessors(peak_node))) == 1:
                                            to_remove.append(peak_node)
                                            # remove succesor nodes of peak nodes
                                            to_remove.append(
                                                list(self.successors(peak_node))[0]
                                            )

                        # remove the nodes from the graph
                        for node in to_remove:
                            self.remove_node(node)
                            total_removed += 1

                pbar.update(1)

        print(total_removed, " nodes have been removed based on best isotopic fit")

    def viz_fragment_coverage_and_charge(self):
        """plotly matrix that displays the coverage of terminal fragment of the peptide sequence
        each fragment type is a row and the columns represent the amino acid position in the peptide sequence
        """

        # 1 for each ion cap types generate empty list of length peptide sequence
        intensity_matrix = {}
        row_names = []
        for ion_type in self.start_ioncaps_types + self.end_ioncaps_types:
            intensity_matrix[ion_type] = [
                0 for i in range(len(self.peptidoform.sequence) + 1)
            ]
            row_names.append(ion_type)

        # print(intensity_matrix)

        # charge length matrix
        charge_length_matrix = {}
        row_names = []
        for charge in range(1, self.max_charge + 1):
            charge_length_matrix[charge] = [
                0 for i in range(len(self.peptidoform.sequence) + 1)
            ]
            row_names.append(charge)

        # print(charge_length_matrix)

        # charge length matrix intern

        charge_length_matrix_intern = {}
        row_names = []
        for charge in range(1, self.max_charge + 1):
            charge_length_matrix_intern[charge] = [
                0 for i in range(len(self.peptidoform.sequence) + 1)
            ]
            row_names.append(charge)

        # 2 for each leaf node that are terminal ions, get the amino acid position, the fragment type and the intensity
        for node in self.nodes:
            if self.nodes[node]["node_type"] == "leaf":
                # if intensity is 0, skip

                frag_len = self.nodes[node]["end_pos"] - self.nodes[node]["start_pos"]
                frag_charge = self.nodes[node]["charge"]
                if self.nodes[node]["intensity"] == 0:
                    continue

                # if n-terminal ion, get the end position
                if self.nodes[node]["frag_dir"] == "N":
                    frag_site = self.nodes[node]["end_pos"]
                    ion_type = self.nodes[node]["end_ioncap"]
                    intensity_matrix[ion_type][frag_site] = self.nodes[node][
                        "intensity"
                    ]
                    charge_length_matrix[frag_charge][frag_len] = (
                        charge_length_matrix[charge][frag_len]
                        + self.nodes[node]["intensity"]
                    )

                # if c-terminal ion, get the start position
                if self.nodes[node]["frag_dir"] == "C":
                    frag_site = self.nodes[node]["start_pos"]
                    ion_type = self.nodes[node]["start_ioncap"]
                    intensity_matrix[ion_type][frag_site] = self.nodes[node][
                        "intensity"
                    ]
                    # charge_length_matrix[frag_charge][frag_len] = charge_length_matrix[charge][frag_len]+self.nodes[node]["intensity"]

                if self.nodes[node]["frag_dir"] == "I":
                    charge_length_matrix_intern[frag_charge][frag_len] = (
                        charge_length_matrix_intern[charge][frag_len]
                        + self.nodes[node]["intensity"]
                    )

        #print(intensity_matrix)

        # log transform all amtrices
        for ion_type in intensity_matrix:
            intensity_matrix[ion_type] = [
                math.log(i + 1) for i in intensity_matrix[ion_type]
            ]

        fig = go.Figure(
            data=go.Heatmap(
                z=list(intensity_matrix.values()),
                x=list(range(len(self.peptidoform.sequence))),
                y=row_names,
            )
        )
        fig.show()

        fig = go.Figure(
            data=go.Heatmap(
                z=list(charge_length_matrix.values()),
                x=list(range(len(self.peptidoform.sequence))),
                y=row_names,
            )
        )
        fig.show()

        fig = go.Figure(
            data=go.Heatmap(
                z=list(charge_length_matrix_intern.values()),
                x=list(range(len(self.peptidoform.sequence))),
                y=row_names,
            )
        )

        fig.show()

        charge_length_matrix = np.array(list(charge_length_matrix.values()))
        # standardize the matrix
        charge_length_matrix = (
            charge_length_matrix - np.mean(charge_length_matrix)
        ) / np.std(charge_length_matrix)

        # Sample charge states and peptide lengths
        charge_states = np.arange(1, charge_length_matrix.shape[0] + 1)
        peptide_lengths = np.arange(1, charge_length_matrix.shape[1] + 1)

        # Create a meshgrid for the input space
        X, Y = np.meshgrid(peptide_lengths, charge_states)
        X_flat = X.flatten()
        Y_flat = Y.flatten()

        # Reshape the charge_length_matrix to a 1D array
        Z = charge_length_matrix.flatten()

        # Define the quadratic kernel
        kernel = ConstantKernel(1.0) * RBF(length_scale=1.0) ** 2 + WhiteKernel(
            noise_level=1e-5
        )

        # Create the Gaussian Process Regressor
        model = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10)

        # Fit the model to the data
        model.fit(np.column_stack((X_flat, Y_flat)), Z)

        # Define new data points for prediction
        new_peptide_lengths = np.arange(1, 100)
        new_charge_states = np.arange(1, 10)
        X_new, Y_new = np.meshgrid(new_peptide_lengths, new_charge_states)
        X_new_flat = X_new.flatten()
        Y_new_flat = Y_new.flatten()

        # Predict the charge probabilities
        predictions, _ = model.predict(
            np.column_stack((X_new_flat, Y_new_flat)), return_std=True
        )

        # Assuming a Gaussian distribution, you can calculate the probabilities using the cumulative distribution function
        from scipy.stats import norm

        probabilities = norm.cdf(predictions)

        # Reshape the probabilities to match the input meshgrid shape
        probabilities = probabilities.reshape(X_new.shape)

        #print("Probabilities:")
        #print(probabilities)

    def model_charge_from_length_probability_3(
        self, smoothing_factor=1, apply_weighting=False, column_smoothing_window=0
    ):
        import numpy as np
        import matplotlib.pyplot as plt
        from scipy.optimize import curve_fit

        # get data
        peptide_lengths = []
        charge_states = []
        intensities = []

        for node in self.nodes:
            if self.nodes[node]["node_type"] == "leaf":
                if self.nodes[node]["frag_dir"] in ["N", "C"]:
                    if self.nodes[node]["intensity"] != 0:
                        peptide_lengths.append(
                            self.nodes[node]["end_pos"] - self.nodes[node]["start_pos"]
                        )
                        charge_states.append(self.nodes[node]["charge"])
                        intensities.append(self.nodes[node]["intensity"])

        #print("peptide_lengths", peptide_lengths)
        #print("charge_states", charge_states)
        #print("intensities", intensities)
        # convert to numpy array
        peptide_lengths = np.array(peptide_lengths)
        charge_states = np.array(charge_states)
        intensities = np.array(intensities)

        num_lengths = len(set(peptide_lengths))
        num_states = len(set(charge_states))
        length_state_counts = {}
        length_state_intensity_sums = {}

        # Count occurrences and sum intensities of each (peptide length, charge state) combination
        for length, state, intensity in zip(
            peptide_lengths, charge_states, intensities
        ):
            if length not in length_state_counts:
                length_state_counts[length] = {}
                length_state_intensity_sums[length] = {}
            if state not in length_state_counts[length]:
                length_state_counts[length][state] = 1
                length_state_intensity_sums[length][state] = intensity
            else:
                length_state_counts[length][state] += 1
                length_state_intensity_sums[length][state] += intensity

        # Estimate probabilities with Laplace Smoothing and custom smoothing factor
        probabilities = {}
        for length, state_counts in length_state_counts.items():
            total_counts = sum(state_counts.values())
            total_intensity = sum(length_state_intensity_sums[length].values())
            probabilities[length] = {}
            for state in range(1, num_states + 1):
                state_count = state_counts.get(state, 0)
                intensity_sum = length_state_intensity_sums[length].get(state, 0)

                # Apply intensity weighting if enabled
                if apply_weighting:
                    prob = (
                        (state_count + smoothing_factor)
                        / (total_counts + smoothing_factor * num_states)
                        * (intensity_sum / total_intensity)
                    )
                else:
                    prob = (state_count + smoothing_factor) / (
                        total_counts + smoothing_factor * num_states
                    )

                probabilities[length][state] = prob

        # Normalize probabilities to ensure they sum up to 1
        for length in probabilities:
            prob_sum = sum(probabilities[length].values())
            for state in probabilities[length]:
                probabilities[length][state] /= prob_sum

        # Apply column smoothing based on moving average
        if column_smoothing_window > 1:
            for state in range(1, num_states + 1):
                smoothed_probs = []
                for length_idx, length in enumerate(probabilities.keys()):
                    start_idx = max(0, length_idx - column_smoothing_window // 2)
                    end_idx = min(
                        num_lengths, length_idx + column_smoothing_window // 2 + 1
                    )
                    smoothed_prob = sum(
                        probabilities[length][state]
                        for length in list(probabilities.keys())[start_idx:end_idx]
                    ) / (end_idx - start_idx)
                    smoothed_probs.append(smoothed_prob)
                for length_idx, length in enumerate(probabilities.keys()):
                    probabilities[length][state] = smoothed_probs[length_idx]

        #print("probabilities", probabilities)

        # convert to long dataframe

        probabilities_long = []
        for length in probabilities:
            for state in probabilities[length]:
                probabilities_long.append([length, state, probabilities[length][state]])

        probabilities_long = pd.DataFrame(
            probabilities_long, columns=["peptide_length", "charge", "probability"]
        )

        # plot as a barplot

        fig = px.bar(
            probabilities_long,
            x="peptide_length",
            y="probability",
            color="charge",
            barmode="group",
        )
        fig.show()

        return probabilities

    def plot_charge_over_fragment_length(self, frag_dir="all"):
        charges = []
        lengths = []
        intensities = []

        for node in self.nodes:
            if self.nodes[node]["node_type"] == "intermediate_2":
                if self.nodes[node]["intensity"] != 0:
                    if frag_dir == "all":
                        charges.append(self.nodes[node]["charge"])
                        lengths.append(
                            self.nodes[node]["end_pos"] - self.nodes[node]["start_pos"]
                        )
                        intensities.append(self.nodes[node]["intensity"])
                    elif frag_dir != "all" and self.nodes[node]["frag_dir"] == frag_dir:
                        charges.append(self.nodes[node]["charge"])
                        lengths.append(
                            self.nodes[node]["end_pos"] - self.nodes[node]["start_pos"]
                        )
                        intensities.append(self.nodes[node]["intensity"])
                    else:
                        continue

        # merge data and sum up intensities
        data = pd.DataFrame(
            {"charge": charges, "length": lengths, "intensity": intensities}
        )
        data = data.groupby(["charge", "length"]).sum().reset_index()

        #print(data)

        # plot as scatter intensity is dot size
        fig = px.scatter(
            data,
            x="length",
            y="charge",
            size="intensity",
            color="charge",
            hover_data=["intensity"],
        )

        fig.show()

    def calculate_number_full_overlap_leaf(self):
        """calculate the number of intermediate_2 nodes whose leaf nodes fully overlap with the leaf nodes of another intermediate_2 node"""

        # get all intermediate_2 nodes
        intermediate_2_nodes = [
            node
            for node in self.nodes
            if self.nodes[node]["node_type"] == "intermediate_2"
        ]

        # iterate over all intermediate_2 nodes
        num_full_overlap = 0
        for node in intermediate_2_nodes:
            # get the leaf nodes
            leaf_nodes = [
                leaf_node
                for leaf_node in self.successors(node)
                if self.nodes[leaf_node]["node_type"] == "leaf"
            ]

            # iterate over all other intermediate_2 nodes
            for other_node in intermediate_2_nodes:
                # get the leaf nodes
                other_leaf_nodes = [
                    leaf_node
                    for leaf_node in self.successors(other_node)
                    if self.nodes[leaf_node]["node_type"] == "leaf"
                ]

                # check if the leaf nodes are the same
                if set(leaf_nodes) == set(other_leaf_nodes):
                    num_full_overlap += 1

        print("##############")
        print("number of intermediate_2 nodes with full overlap:", num_full_overlap)
        print("##############")

        return num_full_overlap

    def print_number_each_node_type(self):
        """print the number of each node type"""

        node_types = [self.nodes[node]["node_type"] for node in self.nodes]

        print("##############")
        print("number of each node type:")
        print(Counter(node_types))
        print("##############")

    def intensity_difference_terminal(self, return_intensities=False):
        """for each fragmentation site, retrieve the intensity of all N and C terminal ions at that site and compute their diference, return a list of all differences"""

        # create empty list corresponding to each fragmentation site
        intensity_n_term = [0 for i in range(len(self.peptidoform.sequence) + 1)]
        intensity_c_term = [0 for i in range(len(self.peptidoform.sequence) + 1)]
        intensity_c_term_intern = [0 for i in range(len(self.peptidoform.sequence) + 1)]
        intensity_n_term_intern = [0 for i in range(len(self.peptidoform.sequence) + 1)]

        #print(intensity_c_term)

        # iterate over all leaf nodes
        for node in self.nodes:
            if self.nodes[node]["node_type"] == "leaf":
                # if intensity is 0, skip
                if self.nodes[node]["intensity"] == 0:
                    continue

                # if n-terminal ion, get the end position
                if self.nodes[node]["frag_dir"] == "N":
                    frag_site = self.nodes[node]["end_pos"]
                    intensity_n_term[frag_site] += self.nodes[node]["intensity"]

                # if c-terminal ion, get the start position
                if self.nodes[node]["frag_dir"] == "C":
                    frag_site = self.nodes[node]["start_pos"]
                    intensity_c_term[frag_site] += self.nodes[node]["intensity"]

                if self.nodes[node]["frag_dir"] == "I":
                    frag_site = self.nodes[node]["start_pos"]
                    intensity_c_term_intern[frag_site] += self.nodes[node]["intensity"]

                    frag_site = self.nodes[node]["end_pos"]
                    intensity_n_term_intern[frag_site] += self.nodes[node]["intensity"]

        # plot the intensity in a bar plot, one bar for each fragmentation site, c-terminal as negative and n-terminal as positive

        if return_intensities:
            return intensity_n_term, intensity_c_term

        fig = go.Figure()
        fig.add_trace(
            go.Bar(
                x=list(range(len(self.peptidoform.sequence))),
                y=intensity_n_term,
                name="n-terminal",
                marker_color="red",
            )
        )
        fig.add_trace(
            go.Bar(
                x=list(range(len(self.peptidoform.sequence))),
                y=[-i for i in intensity_c_term],
                name="c-terminal",
                marker_color="purple",
            )
        )

        fig.add_trace(
            go.Bar(
                x=list(range(len(self.peptidoform.sequence))),
                y=[i / 2 for i in intensity_n_term_intern],
                name="n-terminal intern",
                marker_color="blue",
            )
        )

        fig.add_trace(
            go.Bar(
                x=list(range(len(self.peptidoform.sequence))),
                y=[-i / 2 for i in intensity_c_term_intern],
                name="c-terminal intern",
                marker_color="green",
            )
        )

        fig.update_layout(barmode="group", xaxis_tickangle=-45)
        # fig.show()

        # compute the difference between n-terminal and c-terminal intensities
        intensity_difference = [
            intensity_n_term[i] - intensity_c_term[i]
            for i in range(len(self.peptidoform.sequence))
        ]

        # plot the intensity difference in a bar plot

        fig = go.Figure()
        fig.add_trace(
            go.Bar(
                x=list(range(len(self.peptidoform.sequence))),
                y=intensity_difference,
                name="intensity difference",
                marker_color="indianred",
            )
        )

        fig.update_layout(barmode="group", xaxis_tickangle=-45)
        # fig.show()

        #print(intensity_difference)

        # find most likely pairs
        pairs = self.find_ordered_pairs(intensity_difference)

        # for pair in pairs:
        #     print(f"Indexes: {pair[0]}, {pair[1]}; Sum: {pair[2]}")

        return pairs

    def find_ordered_pairs(self, values):
        n = len(values)
        pairs = []

        for i in range(n - 1):
            for j in range(i + 1, n):
                if values[i] > values[j]:
                    pairs.append((i, j, values[i] - values[j]))

        pairs.sort(key=lambda x: x[2], reverse=True)

        return pairs

    def get_fragment_table(self, peptidoform_index=0):
        """return a table with fragment information at the intermediate_2 nodes level"""

        df = pd.DataFrame(
            columns=[
                "frag_code",
                "start_pos",
                "end_pos",
                "frag_dir",
                "charge",
                "intensity",
                "cosine_similarity",
                "match_iso",
                "n_theo_iso",
                "n_over_term",
                "n_over_tot",
            ]
        )

        # add progress bar
        for node in tqdm(self.nodes, desc="Extracting Annotation"):
            if (
                self.nodes[node]["node_type"] == "intermediate_2"
                and peptidoform_index in self.nodes[node]["peptidoforms"]
            ):
                # calculate_number_ child matched
                match_iso = len(
                    [
                        i
                        for i in self.successors(node)
                        if self.nodes[i]["intensity"] != 0
                    ]
                )
                n_theo_iso = len([i for i in self.successors(node)])

                # number of peaks that overlap with terminal ions
                n_over_term = 0
                n_over_tot = 0
                # for leaf_node in self.successors(node):
                #     overlaping_leafs = self.get_overlapping_leaf_nodes(leaf_node)
                #     for o in overlaping_leafs:
                #         if self.nodes[o]["node_type"] == "leaf":
                #             n_over_tot += 1
                #             if self.nodes[o]["frag_dir"] in ["N", "C"]:
                #                 n_over_term += 1

                df = pd.concat(
                    [
                        df,
                        pd.DataFrame.from_dict(
                            {
                                "frag_code": node,
                                "start_pos": [self.nodes[node]["start_pos"]],
                                "end_pos": [self.nodes[node]["end_pos"]],
                                "frag_dir": [self.nodes[node]["frag_dir"]],
                                "charge": [self.nodes[node]["charge"]],
                                "intensity": [self.nodes[node]["intensity"]],
                                "cosine_similarity": [
                                    self.nodes[node]["cosine_similarity"]
                                ],
                                "match_iso": match_iso,
                                "n_theo_iso": n_theo_iso,
                                "n_over_term": n_over_term,
                                "n_over_tot": n_over_tot,
                            }
                        ),
                    ]
                )

        return df

    def get_ms2_sum_its(self):
        summed_its = np.sum(self.its)

        return summed_its

    def get_peptidoform_indexes(self):
        indexes = []
        for node in self.nodes:
            if "peptidoforms" in self.nodes[node]:
                for i in self.nodes[node]["peptidoforms"]:
                    if i not in indexes:
                        indexes.append(i)
        return indexes

    def get_fragment_table_I0(self, peptidoform_index=0):
        """return a table with fragment information at the intermediate_0 nodes level"""

        df = pd.DataFrame(columns=["start_pos", "end_pos", "frag_codes", "intensity"])

        for node in self.nodes:
            if (
                self.nodes[node]["node_type"] == "intermediate_0"
                and peptidoform_index in self.nodes[node]["peptidoforms"]
            ):
                # get the number of child "intermediate_1" nodes
                n_child_I1 = len(
                    list(
                        [
                            i
                            for i in self.successors(node)
                            if self.nodes[i]["node_type"] == "intermediate_1"
                        ]
                    )
                )
                # get second degree child of intermediate_1 nodes
                n_child_I2 = 0
                # charge of I2 nodes
                charges = []
                cosine_similarities = []
                for i in self.successors(node):
                    if self.nodes[i]["node_type"] == "intermediate_1":
                        n_child_I2 += len(list(self.successors(i)))
                        # get charges and CS val of I2 nodes
                        for j in self.successors(i):
                            charges.append(self.nodes[j]["charge"])
                            cosine_similarities.append(
                                self.nodes[j]["cosine_similarity"]
                            )

                # weighted mean (intensity) of charge:
                charges = np.array(charges)
                cosine_similarities = np.array(cosine_similarities)
                # print("charges", charges)
                # print("cosine_similarities", cosine_similarities)
                # print("intensities", intensities)

                # get average successors mz:
                mzs = []
                for i in self.successors(node):
                    if self.nodes[i]["node_type"] == "intermediate_1":
                        mzs.append(self.nodes[i]["mz"])
                mz = np.mean(mzs)

                df = pd.concat(
                    [
                        df,
                        pd.DataFrame.from_dict(
                            {
                                "start_pos": [self.nodes[node]["start_pos"]],
                                "end_pos": [self.nodes[node]["end_pos"]],
                                "frag_codes": [self.nodes[node]["frag_codes"]],
                                "intensity": [self.nodes[node]["intensity"]],
                                "n_child_I1": n_child_I1,
                                "n_child_I2": n_child_I2,
                                "avg_charge": np.mean(charges),
                                "avg_cosine_similarity": np.mean(cosine_similarities),
                                "mz": mz,
                            }
                        ),
                    ]
                )

        # print(df)
        # how many row below self.min_intermediate_2
        # print("number of rows below min_intermediate_2:", len(df[df["n_child_I2"] < self.min_intermediate_2_evidence]))

        #print head of the dataframe

        print(df.head(10))

        # filter n_child_I2
        df = df[df["n_child_I2"] >= self.min_intermediate_2_evidence]

        # filter n_child_I1
        df = df[df["n_child_I1"] >= self.min_intermediate_1_evidence]

        # make sure that start_pos and end_pos are integers
        df["start_pos"] = df["start_pos"].astype(int)
        df["end_pos"] = df["end_pos"].astype(int)

        return df

    def get_all_fragment_table_I0(self):
        """
        This function returns a table with fragment information at the intermediate_0 nodes level for each peptidoforms
        intensity and mz column for each peptidoform
        (add a suffix to the column name to indicate the peptidoform index)
        """

        # get all peptidoform indexes
        indexes = self.get_peptidoform_indexes()

        print("peptidoform indexes:", indexes)

        list_df = []

        for index in indexes:
            df = self.get_fragment_table_I0(index)

            df["peptidoform_index"] = index

            list_df.append(df)

        #concatenate rows

        df = pd.concat(list_df)




        return df



    def get_leaf_nodes_linked_to_peak(self, peak_node_code):
        """
        This function returns all leaf nodes linked to a specific peak node.

        Parameters:
        peak_node_code (str): The code of the peak node.

        Returns:
        list: A list of leaf nodes linked to the peak node.
        """
        linked_nodes = self.predecessors(peak_node_code)
        leaf_nodes = [
            node for node in linked_nodes if self.nodes[node]["node_type"] == "leaf"
        ]
        return leaf_nodes

    def get_peak_table(self):
        # create atabel where each line is an experimental peak in the sectrum
        df = pd.DataFrame(
            columns=[
                "mz",
                "intensity",
                "list_leaf_nodes",
                "number_leaf_nodes",
                "terminal",
                "internal",
            ]
        )

        for node in self.nodes:
            if self.nodes[node]["node_type"] == "peak":
                # get the leaf nodes
                leaf_nodes = self.get_leaf_nodes_linked_to_peak(node)


                # get the mz and intensity
                mz = self.nodes[node]["mz"]
                intensity = self.nodes[node]["intensity"]
                # get the number of leaf nodes
                number_leaf_nodes = len(leaf_nodes)
                # get the list of leaf nodes
                list_leaf_nodes = leaf_nodes
                # get the number of terminal and internal leaf nodes
                terminal = 0
                internal = 0
                for leaf_node in leaf_nodes:
                    if self.nodes[leaf_node]["frag_dir"] in ["N", "C"]:
                        terminal += 1
                    if self.nodes[leaf_node]["frag_dir"] == "I":
                        internal += 1


                df = pd.concat(
                    [
                        df,
                        pd.DataFrame.from_dict(
                            {
                                "mz": [mz],
                                "intensity": [intensity],
                                "list_leaf_nodes": [list_leaf_nodes],
                                "number_leaf_nodes": [number_leaf_nodes],
                                "terminal": [terminal],
                                "internal": [internal],
                            }
                        ),
                    ]
                )

        return df
