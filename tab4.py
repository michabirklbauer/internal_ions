#!/usr/bin/env python3

import streamlit as st

from util.redirect import st_stdout
from util.spectrumio import read_filtered_spectra
from util.tab4.fraggraph import main as fraggraph_main

from util.constants import DIV_COLOR

def reset_filtered_spectra() -> None:
    st.session_state["rerun_filtered_spectra_reading"] = True

def main(argv = None) -> None:

    ############################################################################
    header = st.subheader("Fraggraph", divider = DIV_COLOR)
    tab4_desc = st.markdown("Description text of tab4.")

    ############################################################################
    data_import_tab4_header = st.subheader("Data Import", divider = DIV_COLOR)

    if "filtered_spectra" in st.session_state:
        st.success(f"Filtered spectra from file \"{st.session_state['filtered_spectra']['name']}\" were successfully loaded!")
        with st.expander("Display filtering parameters:"):
            scans_from_protein_str = ""
            if st.session_state["filtered_spectra"]["filter_params"]["selected_protein"] is None:
                scans_from_protein_str = None
            else:
                if len(st.session_state["filtered_spectra"]["filter_params"]["scans_from_protein"]) > 10:
                    scans_from_protein_str = "\"omitted due to size (>10)\""
                else:
                    scans_from_protein_str = "[" + ", ".join([str(x) for x in st.session_state["filtered_spectra"]["filter_params"]["scans_from_protein"]]) + "]"

            scans_from_peptide_str = ""
            if st.session_state["filtered_spectra"]["filter_params"]["selected_peptide"] is None:
                scans_from_peptide_str = None
            else:
                if len(st.session_state["filtered_spectra"]["filter_params"]["scans_from_peptide"]) > 10:
                    scans_from_peptide_str = "\"omitted due to size (>10)\""
                else:
                    scans_from_peptide_str = "[" + ", ".join([str(x) for x in st.session_state["filtered_spectra"]["filter_params"]["scans_from_peptide"]]) + "]"

            params_str = "\t{\"source_filename\": " + f"{st.session_state['filtered_spectra']['name']}\n" + \
                         "\t \"filter_params\": \n\t\t{\"first_scan\": " + f"{st.session_state['filtered_spectra']['filter_params']['first_scan']}\n" + \
                         "\t\t \"last_scan\": " + f"{st.session_state['filtered_spectra']['filter_params']['last_scan']}\n" + \
                         "\t\t \"min_mz\": " + f"{st.session_state['filtered_spectra']['filter_params']['min_mz']}\n" + \
                         "\t\t \"max_mz\": " + f"{st.session_state['filtered_spectra']['filter_params']['max_mz']}\n" + \
                         "\t\t \"min_rt\": " + f"{st.session_state['filtered_spectra']['filter_params']['min_rt']}\n" + \
                         "\t\t \"max_rt\": " + f"{st.session_state['filtered_spectra']['filter_params']['max_rt']}\n" + \
                         "\t\t \"max_charge\": " + f"{st.session_state['filtered_spectra']['filter_params']['max_charge']}\n" + \
                         "\t\t \"max_isotope\": " + f"{st.session_state['filtered_spectra']['filter_params']['max_isotope']}\n" + \
                         "\t\t \"selected_protein\": " + f"{st.session_state['filtered_spectra']['filter_params']['selected_protein']}\n" + \
                         "\t\t \"scans_from_protein\": " + f"{scans_from_protein_str}\n" + \
                         "\t\t \"selected_peptide\": " + f"{st.session_state['filtered_spectra']['filter_params']['selected_peptide']}\n" + \
                         "\t\t \"scans_from_peptide\": " + f"{scans_from_peptide_str}\n" + \
                         "\t\t}\n\t}"

            filter_dl_desc_text1t4 = st.markdown("**Filter Parameters:**")
            filter_dl_desc_text1t4 = st.text(params_str)

    ############################################################################
        fraggraph_params_header = st.subheader("Fraggraph Parameters", divider = DIV_COLOR)

        fraggraph_params_text = st.markdown("Please specify the parameters used for running Fraggraph.")

        fg_input_col1, fg_input_col2 = st.columns(2)

        with fg_input_col1:
            fg_mzd = st.number_input("Select mzd for combining spectra:",
                                     min_value = 0.0,
                                     max_value = 1.0,
                                     value = 0.01,
                                     step = 0.001,
                                     format = "%0.3f",
                                     help = "mzd.")

        with fg_input_col2:
            fg_cov = st.number_input("Select coverage for peak filtering:",
                                     min_value = 0.0,
                                     max_value = 1.0,
                                     value = 0.1,
                                     step = 0.001,
                                     format = "%0.3f",
                                     help = "cov.")

        fg_peptidoform1 = st.text_input("Specify a peptidoform to consider:",
                                        value = None,
                                        help = "Specify a peptidoform to consider for generating the fragment graph, e.g. " +
                                               "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALRE.",
                                        placeholder = "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALRE")

        fg_peptidoform2 = st.text_input("Optionally, specify a peptidoform to compare to:",
                                        value = None,
                                        help = "Specify a peptidoform to compare the fragment graph of the first peptidoform to, e.g. " +
                                               "ARTKQTARKSTGGKAPRKQLATKAARKSAPAT[-79.966331]GGV[+79.966331]KKPHRYRPGTVALRE.",
                                        placeholder = "ARTKQTARKSTGGKAPRKQLATKAARKSAPAT[-79.966331]GGV[+79.966331]KKPHRYRPGTVALRE")

        fg_fragmentation = st.selectbox("Specify the fragmentation method:",
                                        options = ["ETD", "ETD_deconv", "HCD", "CID",
                                                   "ETHCD", "ETHCD_deconv", "DECOY_high_res", "DECOY_high_res_deconv"],
                                        index = 2,
                                        help = "Specify the fragmention method used.")

        fg_run_l, fg_run_center, fg_run_r = st.columns(3)

        with fg_run_center:
            run_fraggraph = st.button("Run Fraggraph!", use_container_width = True)

        if run_fraggraph:
            fraggraph_main({"mzd": fg_mzd,
                            "cov": fg_cov,
                            "pep1": fg_peptidoform1,
                            "pep2": fg_peptidoform2,
                            "frag": fg_fragmentation})

    else:
        st.info("No spectra selected! Please upload a file in the \"Annotation\" tab and select the " +
                "desired spectra in the \"Spectrum\" tab or upload a file here.")

        filtered_spectrum_file = st.file_uploader("Upload a spectrum file (we assume all spectra are MS2-level):",
                                                  key = "filtered_spectrum_file",
                                                  type = ["mgf"],
                                                  on_change = reset_filtered_spectra,
                                                  help = "Upload a spectrum file to be analyzed in .mgf format.")

        if filtered_spectrum_file is not None:
            trigger_rerun = False
            with st.status("Reading spectra...") as filtered_spectra_reading_status:
                with st_stdout("info"):
                    if "filtered_spectra" not in st.session_state:
                        st.session_state["filtered_spectra"] = read_filtered_spectra(filtered_spectrum_file, filtered_spectrum_file.name)
                        st.session_state["rerun_filtered_spectra_reading"] = False
                        trigger_rerun = True
                    if "filtered_spectra" in st.session_state:
                        if st.session_state["rerun_filtered_spectra_reading"]:
                            st.session_state["filtered_spectra"] = read_filtered_spectra(filtered_spectrum_file, filtered_spectrum_file.name)
                            st.session_state["rerun_filtered_spectra_reading"] = False
                            trigger_rerun = True
                read_filtered_spectra_successfully = st.success("Read all spectra successfully!")
                filtered_spectra_reading_status.update(label = f"Read all spectra from file {st.session_state.filtered_spectrum_file.name} successfully!", state = "complete")
                if trigger_rerun:
                    st.rerun()
