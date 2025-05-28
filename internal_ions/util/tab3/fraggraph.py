import os
import shutil
import random
from datetime import datetime

import streamlit as st
import streamlit.components.v1 as components

# IMPORTANT: Some imports are function level because they require an active R
# installation which will be checked on function call
from ...fraggraph.combine_spectra import combine_spectra
from ...fraggraph.frag_graph_fast import FragGraph
from ...fraggraph.frag_graph_viz import draw_graph3 as draw_graph

from ...util.redirect import st_stdout
from .plots import draw_fragment_coverage_matrix_plotly
from .plots import draw_fragment_coverage_matrix_difference_plotly
from .plots import draw_barplot_intensity_SDI

from ...util.constants import DIV_COLOR


def create_fraggraph(peptidoforms):
    with st.status("Generating the fragmention graph...") as fg_gen_status:
        with st_stdout("info"):
            start_ioncaps_types = st.session_state["selected_ions_cterm"]
            end_ioncaps_types = st.session_state["selected_ions_nterm"]
            msms_tol = st.session_state["tolerance"]
            msms_tol_unit = "da"

            if st.session_state["deconvoluted_spectra"]:
                monoisotopic = True
                max_prec_charge = 1
                max_isotope = 1
                charge_loss = False
            else:
                monoisotopic = False
                charge_loss = st.session_state["charge_reduction"]
                max_prec_charge = "auto" if st.session_state["max_charge_auto"] else st.session_state["max_charge"]
                max_isotope = "auto" if st.session_state["max_isotope_auto"] else st.session_state["max_isotope"]

            print(f"Parameters: {start_ioncaps_types}, {end_ioncaps_types}, {msms_tol}, {msms_tol_unit}, {monoisotopic}, {max_prec_charge}, {max_isotope}, {charge_loss}")

            fg = FragGraph(start_ioncaps_types=start_ioncaps_types,
                           end_ioncaps_types=end_ioncaps_types,
                           msms_tol=msms_tol,
                           msms_tol_unit=msms_tol_unit,
                           monoisotopic=monoisotopic,
                           max_prec_charge=max_prec_charge,
                           max_isotope=max_isotope,
                           charge_loss=charge_loss)

            fg.generate_graph(peptidoforms,
                              st.session_state["consensus_spectrum"]["mz_mean"],
                              st.session_state["consensus_spectrum"]["its_mean"])

        st.success("Successfully generated fragmentation graph!")
        fg_gen_status.update(label=f"Successfully generated fragmentation graph for {', '.join(peptidoforms)}.", state="complete")
    return fg


def single_or_double_fraggraph(peptidoforms: list[str], verbose: bool = False) -> None:
    fg = create_fraggraph(peptidoforms)

    # file storage handling
    tmp_dir_name = "tmp_fragannot_files_653205774"
    if os.path.exists(tmp_dir_name) and os.path.isdir(tmp_dir_name):
        shutil.rmtree(tmp_dir_name)
    os.makedirs(tmp_dir_name)

    output_name_prefix = tmp_dir_name + "/" + datetime.now().strftime("%b-%d-%Y_%H-%M-%S") + "_" + str(random.randint(10000, 99999))

    # visualizing the fragmentation graph
    draw_graph(fg, output_filename=output_name_prefix + "pyvis.html")
    with open(output_name_prefix + "pyvis.html", "r", encoding="utf-8") as f:
        graph = f.read()

    st.subheader("Visualization of the fragmentation graph", divider=DIV_COLOR)
    # st.markdown("Description text")
    components.html(graph, width=900, height=900)

    cols = st.columns(len(peptidoforms))

    for i, col in enumerate(cols):
        with col:
            st.markdown(f"**Fragment Coverage Matrix Peptidoform {i+1}:**")
            frag_mat_plot = draw_fragment_coverage_matrix_plotly(fg, x="intensity", peptidoform_index=i)
            st.plotly_chart(frag_mat_plot, use_container_width=False, height=500, width=500)

    if len(peptidoforms) > 1:
        diff_plot = draw_fragment_coverage_matrix_difference_plotly(fg)
        st.markdown("**Fragment Coverage Matrix Difference:**")
        st.plotly_chart(diff_plot, use_container_width=False, height=600, width=600)

        sdi_plot = draw_barplot_intensity_SDI(fg)
        st.markdown("**Barplot of Intensity SDI:**")
        st.plotly_chart(sdi_plot, use_container_width=False, height=600, width=600)

    try:
        os.remove(output_name_prefix + "pyvis.html")
    except Exception:
        if verbose:
            print("Could not remove file: " + output_name_prefix + "pyvis.html")


def main(argv=None) -> None:

    ############################################################################
    st.subheader("Fraggraph Results", divider=DIV_COLOR)
    # st.markdown("Description of results.")

    params = argv
    params_keys = ["mzd", "cov", "pep1", "pep2"]
    for key in params_keys:
        if key not in params:
            # this actually should never happen!
            st.error("Insufficient parameters for running Fraggraph!", icon="🚨")

    if "consensus_spectrum" not in st.session_state:
        with st.status("Combining spectra to consensus spectrum...") as consensus_status:
            with st_stdout("info"):
                # merge spectra into a single consensus spectrum
                consensus_spectrum = combine_spectra(spectra_list=st.session_state["filtered_spectra"]["spectra"],
                                                     mzd=params["mzd"])
                # filter out peaks that are found in less than cov of the spectra
                consensus_spectrum = consensus_spectrum[consensus_spectrum["cov_spectra"] >= params["cov"]]
                consensus_spectrum = consensus_spectrum.reset_index(drop=True)
                print("Number of peaks in consensus spectrum after filtering: ", len(consensus_spectrum))
                st.session_state["consensus_spectrum"] = consensus_spectrum
            st.success("Successfully created consensus spectrum!")
            consensus_status.update(label="Successfully created consensus spectrum from "
                                          f"filtered spectra from file {st.session_state['filtered_spectra']['name']}!",
                                    state="complete")

    if "consensus_spectrum" in st.session_state:
        # empty strings are handled as None values
        if params["pep1"] == "":
            params["pep1"] = None
        if params["pep2"] == "":
            params["pep2"] = None
        # check cases
        if params["pep1"] is None and params["pep2"] is None:
            st.error("Please select or provide at least one peptidoform/proteoform to run Fraggraph!", icon="🚨")
        elif params["pep1"] is None and params["pep2"] is not None:
            single_or_double_fraggraph([params["pep2"]])
        elif params["pep2"] is None and params["pep1"] is not None:
            single_or_double_fraggraph([params["pep1"]])
        else:
            single_or_double_fraggraph([params["pep1"], params["pep2"]])
