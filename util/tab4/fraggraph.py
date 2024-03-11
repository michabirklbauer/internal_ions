#!/usr/bin/env python3

import os
import shutil
import random
from datetime import datetime

import streamlit as st
import streamlit.components.v1 as components

# IMPORTANT: Some imports are function level because they require an active R
# installation which will be checked on function call
from util.redirect import st_stdout
from fraggraph.combine_spectra import combine_spectra
from fraggraph.frag_graph_fast import FragGraph
from fraggraph.frag_graph_viz import draw_graph3 as draw_graph

from util.constants import DIV_COLOR

def single_fraggraph(pep1: str, frag: str, iso: str, verbose: bool = False) -> None:
    with st.status("Generating the fragmention graph...") as fg_gen_status:
        with st_stdout("info"):
            fg = FragGraph(fragmentation_parameters = frag)
            fg.generate_graph([pep1],
                              st.session_state["consensus_spectrum"]["mz_mean"],
                              st.session_state["consensus_spectrum"]["its_mean"])
        fg_gen_successful = st.success("Successfully generated fragmentation graph!")
        fg_gen_status.update(label = f"Successfully generated fragmentation graph for {pep1}.", state = "complete")

    # file storage handling
    tmp_dir_name = "tmp_fragannot_files_653205774"
    if os.path.exists(tmp_dir_name) and os.path.isdir(tmp_dir_name):
        shutil.rmtree(tmp_dir_name)
        os.makedirs(tmp_dir_name)
    else:
        os.makedirs(tmp_dir_name)

    output_name_prefix = tmp_dir_name + "/" + datetime.now().strftime("%b-%d-%Y_%H-%M-%S") + "_" + str(random.randint(10000, 99999))

    # visualizing the fragmentation graph
    draw_graph(fg, output_filename = output_name_prefix + "pyvis.html")
    with open(output_name_prefix + "pyvis.html", "r", encoding = "utf-8") as f:
        graph = f.read()
        f.close()

    graph_vis_header = st.subheader("Visualization of the fragmentation graph")
    graph_vis_desc = st.markdown("Description text")

    components.html(graph,
                    width = 900,
                    height = 900)

    graph_vis_caption = st.markdown(f"**Figure 1:** caption.")

    # check for R installation
    if "R_HOME" in os.environ:
        has_R_installed = True
    else:
        with open("R_HOME.config", "r", encoding = "utf-8") as f:
            R_path = f.read().strip()
            f.close()
        if os.path.exists(R_path):
            os.environ["R_HOME"] = R_path
            has_R_installed = True
        else:
            has_R_installed = False
            st.error("No valid R installation found! Please refer to the documentation for troubleshooting!", icon = "🚨")

    if has_R_installed:
        # other plots
        plots_header = st.subheader("Other Plots")
        plots_desc = st.markdown("Description text")

        # plotting the spectrum
        from fraggraph.frag_graph_report import plot_spectrum

        plot1_title = st.markdown("**Experimental Spectrum:**")
        plot1 = st.pyplot(plot_spectrum(fg), use_container_width = True)
        plot1_caption = st.markdown(f"**Figure 2:** caption.")

        # formatting
        plot_col1, plot_col2 = st.columns(2)

        # plotting distribution
        from plotnine import ggplot
        from fraggraph.frag_graph_report import percentage_intensity_distribution

        with plot_col1:
            plot2_title = st.markdown("**Intensity Distribution:**")
            plot2 = st.pyplot(ggplot.draw(percentage_intensity_distribution(fg, filename = None)), use_container_width = True)
            plot2_caption = st.markdown(f"**Figure 3:** caption.")

        # fragment coverage matrix
        from fraggraph.frag_graph_report import draw_fragment_coverage_matrix

        with plot_col2:
            plot3_title = st.markdown("**Fragment Coverage Matrix:**")
            plot3 = st.pyplot(ggplot.draw(draw_fragment_coverage_matrix(fg,
                                                                        x = "avg_cosine_similarity",
                                                                        x_min = 0.7,
                                                                        x_max = 1)),
                              use_container_width = True)
            plot3_caption = st.markdown(f"**Figure 4:** caption.")

    # file storage clean up
    try:
        os.remove(output_name_prefix + "pyvis.html")
    except Exception as e:
        if verbose:
            print("Could not remove file: " + output_name_prefix + "pyvis.html")

    return

def double_fraggraph(pep1: str, pep2: str, frag: str, iso: str, verbose: bool = False) -> None:
    with st.status("Generating the fragmention graph...") as fg_gen_status:
        with st_stdout("info"):
            fg = FragGraph(fragmentation_parameters = frag)
            fg.generate_graph([pep1, pep2],
                              st.session_state["consensus_spectrum"]["mz_mean"],
                              st.session_state["consensus_spectrum"]["its_mean"])
        fg_gen_successful = st.success("Successfully generated fragmentation graph!")
        fg_gen_status.update(label = f"Successfully generated fragmentation graph for {pep1} and {pep2}.", state = "complete")

    # file storage handling
    tmp_dir_name = "tmp_fragannot_files_653205774"
    if os.path.exists(tmp_dir_name) and os.path.isdir(tmp_dir_name):
        shutil.rmtree(tmp_dir_name)
        os.makedirs(tmp_dir_name)
    else:
        os.makedirs(tmp_dir_name)

    output_name_prefix = tmp_dir_name + "/" + datetime.now().strftime("%b-%d-%Y_%H-%M-%S") + "_" + str(random.randint(10000, 99999))

    # visualizing the fragmentation graph
    draw_graph(fg, output_filename = output_name_prefix + "pyvis.html")
    with open(output_name_prefix + "pyvis.html", "r", encoding = "utf-8") as f:
        graph = f.read()
        f.close()

    graph_vis_header = st.subheader("Visualization of the fragmentation graph")
    graph_vis_desc = st.markdown("Description text")

    components.html(graph,
                    width = 900,
                    height = 900)

    graph_vis_caption = st.markdown(f"**Figure 1:** caption.")

    # check for R installation
    if "R_HOME" in os.environ:
        has_R_installed = True
    else:
        with open("R_HOME.config", "r", encoding = "utf-8") as f:
            R_path = f.read().strip()
            f.close()
        if os.path.exists(R_path):
            os.environ["R_HOME"] = R_path
            has_R_installed = True
        else:
            has_R_installed = False
            st.error("No valid R installation found! Please refer to the documentation for troubleshooting!", icon = "🚨")

    if has_R_installed:
        # other plots
        plots_header = st.subheader("Other Plots")
        plots_desc = st.markdown("Description text")

        # plotting the spectrum
        from fraggraph.frag_graph_report import plot_spectrum

        plot1_title = st.markdown("**Experimental Spectrum:**")
        plot1 = st.pyplot(plot_spectrum(fg), use_container_width = True)
        plot1_caption = st.markdown(f"**Figure 2:** caption.")

        # formatting
        plot_col1, plot_col2 = st.columns(2)

        # fragment coverage matrix pep1
        from plotnine import ggplot
        from fraggraph.frag_graph_report import draw_fragment_coverage_matrix

        with plot_col1:
            plot2_title = st.markdown("**Fragment Coverage Matrix Peptidoform 1:**")
            plot2 = st.pyplot(ggplot.draw(draw_fragment_coverage_matrix(fg,
                                                                        x = "intensity",
                                                                        x_min = 1,
                                                                        x_max = 12,
                                                                        peptidoform_index = 0)),
                              use_container_width = True)
            plot2_caption = st.markdown(f"**Figure 3:** caption.")

        # fragment coverage matrix pep2
        from plotnine import ggplot
        from fraggraph.frag_graph_report import draw_fragment_coverage_matrix

        with plot_col2:
            plot3_title = st.markdown("**Fragment Coverage Matrix Peptidoform 2:**")
            plot3 = st.pyplot(ggplot.draw(draw_fragment_coverage_matrix(fg,
                                                                        x = "intensity",
                                                                        x_min = 1,
                                                                        x_max = 12,
                                                                        peptidoform_index = 1)),
                              use_container_width = True)
            plot3_caption = st.markdown(f"**Figure 4:** caption.")

        #### TODO ####
        # the remaining plots would need a function that can determine modification sites
        # based on the peptidoform str
        st.error("SiteDeterminingPlots: NotImplementedException", icon = "🚨")

    # file storage clean up
    try:
        os.remove(output_name_prefix + "pyvis.html")
    except Exception as e:
        if verbose:
            print("Could not remove file: " + output_name_prefix + "pyvis.html")

    return

def main(argv = None) -> None:

    ############################################################################
    header = st.subheader("Fraggraph Results", divider = DIV_COLOR)
    tab3_desc = st.markdown("Description of results.")

    params = argv
    params_keys = ["mzd", "cov", "pep1", "pep2", "frag", "iso"]
    for key in params_keys:
        if key not in params:
            # this actually should never happen!
            st.error(f"Insufficient parameters for running Fraggraph!", icon = "🚨")

    with st.status("Combining spectra to consensus spectrum...") as consensus_status:
        with st_stdout("info"):
            # merge spectra into a single consensus spectrum
            consensus_spectrum = combine_spectra(spectra_list = st.session_state["filtered_spectra"]["spectra"],
                                                 mzd = params["mzd"])
            # filter out peaks that are found in less than cov of the spectra
            consensus_spectrum = consensus_spectrum[consensus_spectrum["cov_spectra"] >= params["cov"]]
            consensus_spectrum = consensus_spectrum.reset_index(drop = True)
            print("Number of peaks in consensus spectrum after filtering: ", len(consensus_spectrum))
            st.session_state["consensus_spectrum"] = consensus_spectrum
        consensus_spectrum_success = st.success("Successfully created consensus spectrum!")
        consensus_status.update(label = f"Successfully created consensus spectrum from " +
                                        f"filtered spectra from file {st.session_state['filtered_spectra']['name']}!",
                                state = "complete")

    if "consensus_spectrum" in st.session_state:
        # empty strings are handled as None values
        if params["pep1"] == "":
            params["pep1"] = None
        if params["pep2"] == "":
            params["pep2"] = None
        # check cases
        if params["pep1"] is None and params["pep2"] is None:
            st.error(f"Please select at least one peptide/peptidoform to run Fraggraph!", icon = "🚨")
        elif params["pep1"] is None and params["pep2"] is not None:
            single_fraggraph(params["pep2"], params["frag"], params["iso"])
        elif params["pep2"] is None and params["pep1"] is not None:
            single_fraggraph(params["pep1"], params["frag"], params["iso"])
        else:
            double_fraggraph(params["pep1"], params["pep2"], params["frag"], params["iso"])
