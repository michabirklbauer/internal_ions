#!/usr/bin/env python3

import json

import streamlit as st

from util.redirect import st_stdout
from util.spectrumio import filter_spectra
from util.streamlit_utils import spectra_to_mgf_stream
from util.tab3.plots import plot_spectrum

from util.constants import DIV_COLOR

def main(argv = None) -> None:

    ############################################################################
    header = st.subheader("Spectrum Tab", divider = DIV_COLOR)
    tab3_desc = st.markdown("Description text of tab3.")

    ############################################################################
    data_import_tab3_header = st.subheader("Data Import", divider = DIV_COLOR)

    if "spectra" in st.session_state:
        st.success(f"Spectra from file \"{st.session_state['spectra']['name']}\" were successfully loaded!")

    ############################################################################
        spectrum_viewer_header = st.subheader("Spectrum Viewer", divider = DIV_COLOR)

        scan_nr = st.selectbox("Select a scan number to display the corresponding spectrum:",
                               st.session_state["spectra"]["spectra"].keys(),
                               index = None,
                               help = "Select a scan number that is available in the uploaded .mgf from the drop down to display the corresponding spectrum.")

        if scan_nr is not None:
            spectrum_plot = plot_spectrum(st.session_state["spectra"]["spectra"][scan_nr]["peaks"].values(),
                                          st.session_state["spectra"]["spectra"][scan_nr]["peaks"].keys())
            st.plotly_chart(spectrum_plot)

    ############################################################################
        spectrum_selection_header = st.subheader("Spectrum Selector", divider = DIV_COLOR)

        spectrum_selection_text = st.markdown("Using the following fields you can filter your spectra for further analyis.")

        spec_sel_col1, spec_sel_col2 = st.columns(2)

        with spec_sel_col1:
            first_scan = st.selectbox("Select the first scan number to analyse:",
                                      st.session_state["spectra"]["spectra"].keys(),
                                      index = 0,
                                      help = "Select the scan number where analysis should start.")

            min_mz = st.number_input("Select the minimum m/z for a spectrum:",
                                     min_value = 0.0,
                                     max_value = 10000.0,
                                     value = 0.0,
                                     step = 0.01,
                                     help = "The minimum m/z.")

            min_rt = st.number_input("Select the minimum retention time for a spectrum:",
                                     min_value = min([float(s["rt"]) for s in st.session_state["spectra"]["spectra"].values()]),
                                     max_value = max([float(s["rt"]) for s in st.session_state["spectra"]["spectra"].values()]),
                                     value = min([float(s["rt"]) for s in st.session_state["spectra"]["spectra"].values()]),
                                     step = 0.01,
                                     help = "The minimum rt.")

            max_charge = st.number_input("Select the maximum charge for a spectrum:",
                                         min_value = 1,
                                         max_value = 20,
                                         value = 4,
                                         step = 1,
                                         help = "The maximum charge.")

        with spec_sel_col2:
            last_scan = st.selectbox("Select the last scan number to analyse:",
                                     st.session_state["spectra"]["spectra"].keys(),
                                     index = len(st.session_state["spectra"]["spectra"].keys()) - 1,
                                     help = "Select the scan number where analysis should end.")

            max_mz = st.number_input("Select the maximum m/z for a spectrum:",
                                     min_value = 0.0,
                                     max_value = 10000.0,
                                     value = 10000.0,
                                     step = 0.01,
                                     help = "The maximum m/z.")

            max_rt = st.number_input("Select the maximum retention time for a spectrum:",
                                     min_value = min([float(s["rt"]) for s in st.session_state["spectra"]["spectra"].values()]),
                                     max_value = max([float(s["rt"]) for s in st.session_state["spectra"]["spectra"].values()]),
                                     value = max([float(s["rt"]) for s in st.session_state["spectra"]["spectra"].values()]),
                                     step = 0.01,
                                     help = "The maximum rt.")

            max_isotope = st.number_input("Select the maximum isotope for a spectrum:",
                                          min_value = 1,
                                          max_value = 20,
                                          value = 4,
                                          step = 1,
                                          help = "The maximum isotope.")

        if "identifications" in st.session_state:
            st.success(f"Identifications from file \"{st.session_state['identifications']['name']}\" were successfully loaded!")

            protein_col, peptide_col = st.columns(2)

            with protein_col:
                scans_from_protein = st.selectbox("Select scans from a specific protein:",
                                                  st.session_state["identifications"]["proteins"].keys(),
                                                  key = "selected_protein_scans",
                                                  index = None,
                                                  help = "Select a protein of interest.")

            with peptide_col:
                scans_from_peptide = st.selectbox("Select scans from a specific peptide:",
                                                  st.session_state["identifications"]["peptides"].keys(),
                                                  key = "selected_peptide_scans",
                                                  index = None,
                                                  help = "Select a peptide of interest.")

        else:
            st.info("No identifications file was provided! Filtering based on proteins/peptides not available unless a identifications file is uploaded in the \"Annotation\" tab!")

    ############################################################################
        filter_spectra_header = st.subheader("Filter Spectra", divider = DIV_COLOR)

        l1, l2, center_button, r1, r2 = st.columns(5)

        with center_button:
            run_filter = st.button("Run filter!", use_container_width = True)

        if run_filter:
            if st.session_state.spectrum_file is not None:
                with st.status("Filtering spectra...") as filter_status:
                    scans_from_protein_val = None
                    scans_from_protein_list = [i for i in range(int(first_scan), int(last_scan) + 1)]
                    scans_from_peptide_val = None
                    scans_from_peptide_list = [i for i in range(int(first_scan), int(last_scan) + 1)]
                    if "selected_protein_scans" in st.session_state:
                        if st.session_state["selected_protein_scans"] is not None:
                            scans_from_protein_val = st.session_state["selected_protein_scans"]
                            scans_from_protein_list = st.session_state["identifications"]["proteins"][scans_from_protein_val]
                    if "selected_peptide_scans" in st.session_state:
                        if st.session_state["selected_peptide_scans"] is not None:
                            scans_from_peptide_val = st.session_state["selected_peptide_scans"]
                            scans_from_peptide_list = st.session_state["identifications"]["peptides"][scans_from_peptide_val]
                    filter_params = {"first_scan": first_scan,
                                     "last_scan": last_scan,
                                     "min_mz": min_mz,
                                     "max_mz": max_mz,
                                     "min_rt": min_rt,
                                     "max_rt": max_rt,
                                     "max_charge": max_charge,
                                     "max_isotope": max_isotope,
                                     "selected_protein": scans_from_protein_val,
                                     "scans_from_protein": scans_from_protein_list,
                                     "selected_peptide": scans_from_peptide_val,
                                     "scans_from_peptide": scans_from_peptide_list}
                    with st_stdout("info"):
                        st.session_state["filtered_spectra"] = filter_spectra(st.session_state.spectrum_file,
                                                                              filter_params,
                                                                              st.session_state.spectrum_file.name)
                        filter_spectra_successfully = st.success("Filtered all spectra successfully!")
                    filter_status.update(label = "Successfully finished filtering spectra.", state = "complete")
            else:
                st.error("Error reading spectra file! Please re-upload your file in the \"Annotation\" tab!", icon = "🚨")

    else:
        st.error("No spectra file uploaded! Please upload a file in the \"Annotation\" tab!", icon = "🚨")

    ############################################################################
    if "filtered_spectra" in st.session_state:
        download_header = st.subheader("Download Filtered Spectra", divider = DIV_COLOR)

        filter_dl_desc_text1 = st.markdown("**Filtered spectra based on following criteria:**")
        filter_dl_desc_text1 = st.markdown("**Input Filename:**")
        filter_dl_desc_text1 = st.text(f"{st.session_state['filtered_spectra']['name']}")

        scans_from_protein_str = ""
        if st.session_state["filtered_spectra"]["filter_params"]["selected_protein"] is None:
            scans_from_protein_str = None
        else:
            if len(st.session_state["filtered_spectra"]["filter_params"]["scans_from_protein"]) > 10:
                scans_from_protein_str = "\"omitted due to size (>10)\""
            else:
                scans_from_protein_str = "[" + ", ".join(st.session_state["filtered_spectra"]["filter_params"]["scans_from_protein"]) + "]"

        scans_from_peptide_str = ""
        if st.session_state["filtered_spectra"]["filter_params"]["selected_peptide"] is None:
            scans_from_peptide_str = None
        else:
            if len(st.session_state["filtered_spectra"]["filter_params"]["scans_from_peptide"]) > 10:
                scans_from_peptide_str = "\"omitted due to size (>10)\""
            else:
                scans_from_peptide_str = "[" + ", ".join(st.session_state["filtered_spectra"]["filter_params"]["scans_from_peptide"]) + "]"

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

        filter_dl_desc_text1 = st.markdown("**Filter Parameters:**")
        filter_dl_desc_text1 = st.text(params_str)

        dl_l1, dl_r1 = st.columns(2)

        with dl_l1:
            mgf_spectra = st.download_button(label = "Download spectra in .mgf format!",
                                             data = spectra_to_mgf_stream(st.session_state["filtered_spectra"]),
                                             file_name = "filtered_spectra.mgf.txt",
                                             mime = "text/plain",
                                             help = "Download the filtered spectra in .mgf file format.",
                                             type = "primary",
                                             use_container_width = True)

        with dl_r1:
            json_meta = st.download_button(label = "Download meta-data in .json format!",
                                           data = json.dumps({"source_filename": st.session_state["filtered_spectra"]["name"],
                                                              "filter_params": st.session_state["filtered_spectra"]["filter_params"]}),
                                           file_name = "filtered_spectra_meta.json",
                                           mime = "text/json",
                                           help = "Download meta-data of the filtered spectra in .json file format.",
                                           type = "primary",
                                           use_container_width = True)

cmt = \
"""
    if "dataframes" in st.session_state:
        if st.session_state["dataframes"] is not None:

            st.write("Files are available.")
            result_copy = st.session_state["dataframes"][0]
            result_copy["id_extracted"] = [int(item.split("=")[-1]) for item in result_copy["spectrum_id"]] # put in parser??

            unique_ids = list(set(result_copy["id_extracted"]))

            selected_spectrum = st.slider("Select a spectrum", min_value=min(unique_ids), max_value=max(unique_ids))

            ints = result_copy[result_copy['id_extracted'] == selected_spectrum]["frag_intensity"]
            mzs = result_copy[result_copy['id_extracted'] == selected_spectrum]["frag_mz"]
            # Choose from list of avalible spectra indices?
            spectrum_plot = plot_spectrum(ints, mzs) # input selected_spectrum
            st.plotly_chart(spectrum_plot)


    else:
        st.write("File is not available. Upload your files in the annotation tab.")
"""
