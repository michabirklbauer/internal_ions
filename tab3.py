#!/usr/bin/env python3

import streamlit as st

from util.tab3.plots import plot_spectrum

def main(argv = None) -> None:

    ############################################################################
    header = st.subheader("Spectrum Tab", divider = "rainbow")

    if "spectra" in st.session_state:
        st.success(f"Spectra from file \"{st.session_state['spectra']['name']}\" were successfully loaded!")

    ############################################################################
        spectrum_viewer_header = st.subheader("Spectrum Viewer", divider = "rainbow")

        scan_nr = st.selectbox("Select a scan number to display the corresponding spectrum:",
                               st.session_state["spectra"]["spectra"].keys(),
                               index = None,
                               help = "Select a scan number that is available in the uploaded .mgf from the drop down to display the corresponding spectrum.")

        if scan_nr is not None:
            spectrum_plot = plot_spectrum(st.session_state["spectra"]["spectra"][scan_nr]["peaks"].values(),
                                          st.session_state["spectra"]["spectra"][scan_nr]["peaks"].keys())
            st.plotly_chart(spectrum_plot)

    ############################################################################
        spectrum_selection_header = st.subheader("Spectrum Selector", divider = "rainbow")

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
        filter_spectra_header = st.subheader("Filter Spectra", divider = "rainbow")

        l1, l2, center_button, r1, r2 = st.columns(5)

        with center_button:
            run_filter = st.button("Run filter!", use_container_width = True)

        if run_filter:
            #### TODO ####
            # get all params
            # implement filter function
            # do filtering
            st.error("NotImplementedException", icon = "🚨")

    else:
        st.error("No spectra file uploaded! Please upload a file in the \"Annotation\" tab!", icon = "🚨")

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
