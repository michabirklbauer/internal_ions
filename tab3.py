#!/usr/bin/env python3

import json

import streamlit as st

from streamlit_plotly_events import plotly_events

from util.redirect import st_stdout
from util.spectrumio import filter_spectra
from util.streamlit_utils import spectra_to_mgf_stream
from util.tab3.plots import plot_spectrum
from util.tab3.plots import plot_consensus_spectrum
from util.tab3.plots import plot_spectra_chromatogram
from fraggraph.combine_spectra import combine_spectra
from util.tab4.fraggraph import main as fraggraph_main

from util.constants import DIV_COLOR

def main(argv = None) -> None:

    ############################################################################
    
    tab3_desc = st.markdown("Reannotation of fragments on a given spectrum (or merged spectra) using fragmentation graph.")

    ############################################################################
    data_import_tab3_header = st.subheader("Data Import", divider = DIV_COLOR)

    if "spectra" in st.session_state:
        st.success(f"Spectra from file \"{st.session_state['spectra']['name']}\" were successfully loaded!")

    

    ############################################################################
        spectrum_selection_header = st.subheader("Spectrum Selector", divider = DIV_COLOR)
        spectrum_selection_text = st.markdown("Using the following fields you can filter your spectra for further analyis.")

        
        
        # Plot chromatogram
        spectra_chromatogram = plot_spectra_chromatogram(st.session_state["spectra"]["spectra"])
        spectra_chromatogram_selection = plotly_events(spectra_chromatogram, select_event=True)
        
        if spectra_chromatogram_selection is not None and len(spectra_chromatogram_selection) > 0: 
            #get min and max rt
            min_rt_select = min([rt["x"] for rt in spectra_chromatogram_selection])
            max_rt_select = max([rt["x"] for rt in spectra_chromatogram_selection])
            #get min and max mz
            min_mz_select = min([mz["y"] for mz in spectra_chromatogram_selection])
            max_mz_select = max([mz["y"] for mz in spectra_chromatogram_selection])
            
            #change value in the corresponding number_input fields
            st.session_state["min_rt_filter"] = min_rt_select
            st.session_state["max_rt_filter"] = max_rt_select
            st.session_state["min_mz_filter"] = min_mz_select
            st.session_state["max_mz_filter"] = max_mz_select
            
            
            
            
        
        spec_sel_col1, spec_sel_col2 = st.columns(2)
        
        with spec_sel_col1:
            first_scan = st.selectbox("Select the first scan number to analyse:",
                                      st.session_state["spectra"]["spectra"].keys(),
                                      index = 0,
                                      help = "Select the scan number where analysis should start.",
                                      key="first_scan_filter")

            min_mz = st.number_input("Select the minimum m/z for a spectrum:",
                                     min_value = 0.0,
                                     max_value = 10000.0,
                                     value = 0.0,
                                     step = 0.01,
                                     help = "The minimum m/z.", 
                                     key="min_mz_filter")

            min_rt = st.number_input("Select the minimum retention time for a spectrum:",
                                     min_value = min([float(s["rt"]) for s in st.session_state["spectra"]["spectra"].values()]),
                                     max_value = max([float(s["rt"]) for s in st.session_state["spectra"]["spectra"].values()]),
                                     value = min([float(s["rt"]) for s in st.session_state["spectra"]["spectra"].values()]),
                                     step = 0.01,
                                     help = "The minimum rt.",
                                     key="min_rt_filter")

            # max_charge = st.number_input("Select the maximum charge for a spectrum:",
            #                              min_value = 1,
            #                              max_value = 20,
            #                              value = 4,
            #                              step = 1,
            #                              help = "The maximum charge.",
            #                              key="max_charge_filter")

        with spec_sel_col2:
            last_scan = st.selectbox("Select the last scan number to analyse:",
                                     st.session_state["spectra"]["spectra"].keys(),
                                     index = len(st.session_state["spectra"]["spectra"].keys()) - 1,
                                     help = "Select the scan number where analysis should end.",
                                     key="last_scan_filter")

            max_mz = st.number_input("Select the maximum m/z for a spectrum:",
                                     min_value = 0.0,
                                     max_value = 10000.0,
                                     value = 10000.0,
                                     step = 0.01,
                                     help = "The maximum m/z.",
                                     key="max_mz_filter")

            max_rt = st.number_input("Select the maximum retention time for a spectrum:",
                                     min_value = min([float(s["rt"]) for s in st.session_state["spectra"]["spectra"].values()]),
                                     max_value = max([float(s["rt"]) for s in st.session_state["spectra"]["spectra"].values()]),
                                     value = max([float(s["rt"]) for s in st.session_state["spectra"]["spectra"].values()]),
                                     step = 0.01,
                                     help = "The maximum rt.",
                                     key="max_rt_filter")

            # max_isotope = st.number_input("Select the maximum isotope for a spectrum:",
            #                               min_value = 1,
            #                               max_value = 20,
            #                               value = 4,
            #                               step = 1,
            #                               help = "The maximum isotope.",
            #                               key="max_isotope_filter")

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
            run_filter = st.button("Filter spectra and create consensus spectrum", use_container_width = True)
            

        with r1:
            fg_mzd = st.number_input("Select mzd for combining spectra:",
                                     min_value = 0.0,
                                     max_value = 1.0,
                                     value = 0.01,
                                     step = 0.001,
                                     format = "%0.3f",
                                     help = "mzd.")

        with r2:
            fg_cov = st.number_input("Select coverage for peak filtering:",
                                     min_value = 0.0,
                                     max_value = 1.0,
                                     value = 0.25,
                                     step = 0.001,
                                     format = "%0.3f",
                                     help = "cov.")

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
                                     #"max_charge": max_charge,
                                     #"max_isotope": max_isotope,
                                     "selected_protein": scans_from_protein_val,
                                     "scans_from_protein": list(scans_from_protein_list),
                                     "selected_peptide": scans_from_peptide_val,
                                     "scans_from_peptide": list(scans_from_peptide_list)}
                    with st_stdout("info"):
                        st.session_state["filtered_spectra"] = filter_spectra(st.session_state.spectrum_file,
                                                                              filter_params,
                                                                              st.session_state.spectrum_file.name)
                        

                    filter_status.update(label = "Successfully finished filtering spectra.", state = "complete")
                    
            else:
                st.error("Error reading spectra file! Please re-upload your file in the \"Annotation\" tab!", icon = "🚨")

            # build consensus spectrum
            if "filtered_spectra" in st.session_state:         
                with st.status("Combining spectra to consensus spectrum...") as consensus_status:
                    with st_stdout("info"):
                        # merge spectra into a single consensus spectrum
                        consensus_spectrum = combine_spectra(spectra_list = st.session_state["filtered_spectra"]["spectra"],
                                                            mzd = fg_mzd)
                        print("Number of peaks in consensus spectrum before filtering: ", len(consensus_spectrum), "\n")
                        # filter out peaks that are found in less than cov of the spectra
                        consensus_spectrum = consensus_spectrum[consensus_spectrum["cov_spectra"] >=  fg_cov]
                        consensus_spectrum = consensus_spectrum.reset_index(drop = True)
                        print("Number of peaks in consensus spectrum after filtering: ", len(consensus_spectrum))
                        st.session_state["consensus_spectrum"] = consensus_spectrum
                    consensus_spectrum_success = st.success("Successfully created consensus spectrum!")
                    consensus_status.update(label = f"Successfully created consensus spectrum from " +
                                                    f"filtered spectra from file {st.session_state['filtered_spectra']['name']}!",
                                            state = "complete")
                                    

    else:
        st.error("No spectra file uploaded! Please upload a file in the \"Annotation\" tab!", icon = "🚨")
        
        
    ############################################################################
    
    spectrum_viewer_header = st.subheader("Spectrum Viewer", divider = DIV_COLOR)

    if "filtered_spectra" in st.session_state:
        spectrum_plot = plot_consensus_spectrum(st.session_state["consensus_spectrum"]["its_mean"],
                                        st.session_state["consensus_spectrum"]["mz_mean"],
                                        st.session_state["consensus_spectrum"]["cov_spectra"]
                                        )
        st_spectrum_plot = st.plotly_chart(spectrum_plot, use_container_width = True)
        st_spectrum_plot_caption = st.markdown(f"**Figure 1:** Displaying consensus spetrum.")


    ############################################################################
    # if "filtered_spectra" in st.session_state:
    #     download_header = st.subheader("Download Filtered Spectra", divider = DIV_COLOR)

    #     filter_dl_desc_text1 = st.markdown("**Filtered spectra based on following criteria:**")
    #     filter_dl_desc_text1 = st.markdown("**Input Filename:**")
    #     filter_dl_desc_text1 = st.text(f"{st.session_state['filtered_spectra']['name']}")

    #     scans_from_protein_str = ""
    #     if st.session_state["filtered_spectra"]["filter_params"]["selected_protein"] is None:
    #         scans_from_protein_str = None
    #     else:
    #         if len(st.session_state["filtered_spectra"]["filter_params"]["scans_from_protein"]) > 10:
    #             scans_from_protein_str = "\"omitted due to size (>10)\""
    #         else:
    #             scans_from_protein_str = "[" + ", ".join([str(x) for x in st.session_state["filtered_spectra"]["filter_params"]["scans_from_protein"]]) + "]"

    #     scans_from_peptide_str = ""
    #     if st.session_state["filtered_spectra"]["filter_params"]["selected_peptide"] is None:
    #         scans_from_peptide_str = None
    #     else:
    #         if len(st.session_state["filtered_spectra"]["filter_params"]["scans_from_peptide"]) > 10:
    #             scans_from_peptide_str = "\"omitted due to size (>10)\""
    #         else:
    #             scans_from_peptide_str = "[" + ", ".join([str(x) for x in st.session_state["filtered_spectra"]["filter_params"]["scans_from_peptide"]]) + "]"

    #     params_str = "\t{\"source_filename\": " + f"{st.session_state['filtered_spectra']['name']}\n" + \
    #                  "\t \"filter_params\": \n\t\t{\"first_scan\": " + f"{st.session_state['filtered_spectra']['filter_params']['first_scan']}\n" + \
    #                  "\t\t \"last_scan\": " + f"{st.session_state['filtered_spectra']['filter_params']['last_scan']}\n" + \
    #                  "\t\t \"min_mz\": " + f"{st.session_state['filtered_spectra']['filter_params']['min_mz']}\n" + \
    #                  "\t\t \"max_mz\": " + f"{st.session_state['filtered_spectra']['filter_params']['max_mz']}\n" + \
    #                  "\t\t \"min_rt\": " + f"{st.session_state['filtered_spectra']['filter_params']['min_rt']}\n" + \
    #                  "\t\t \"max_rt\": " + f"{st.session_state['filtered_spectra']['filter_params']['max_rt']}\n" + \
    #                  "\t\t \"max_charge\": " + f"{st.session_state['filtered_spectra']['filter_params']['max_charge']}\n" + \
    #                  "\t\t \"max_isotope\": " + f"{st.session_state['filtered_spectra']['filter_params']['max_isotope']}\n" + \
    #                  "\t\t \"selected_protein\": " + f"{st.session_state['filtered_spectra']['filter_params']['selected_protein']}\n" + \
    #                  "\t\t \"scans_from_protein\": " + f"{scans_from_protein_str}\n" + \
    #                  "\t\t \"selected_peptide\": " + f"{st.session_state['filtered_spectra']['filter_params']['selected_peptide']}\n" + \
    #                  "\t\t \"scans_from_peptide\": " + f"{scans_from_peptide_str}\n" + \
    #                  "\t\t}\n\t}"

    #     filter_dl_desc_text1 = st.markdown("**Filter Parameters:**")
    #     filter_dl_desc_text1 = st.text(params_str)

    #     dl_l1, dl_r1 = st.columns(2)

    #     with dl_l1:
    #         mgf_spectra = st.download_button(label = "Download spectra in .mgf format!",
    #                                          data = spectra_to_mgf_stream(st.session_state["filtered_spectra"]),
    #                                          file_name = "filtered_spectra.mgf.txt",
    #                                          mime = "text/plain",
    #                                          help = "Download the filtered spectra in .mgf file format.",
    #                                          type = "primary",
    #                                          use_container_width = True)

    #     with dl_r1:
    #         json_meta = st.download_button(label = "Download meta-data in .json format!",
    #                                        data = json.dumps({"source_filename": st.session_state["filtered_spectra"]["name"],
    #                                                           "filter_params": st.session_state["filtered_spectra"]["filter_params"]}),
    #                                        file_name = "filtered_spectra_meta.json",
    #                                        mime = "text/json",
    #                                        help = "Download meta-data of the filtered spectra in .json file format.",
    #                                        type = "primary",
    #                                        use_container_width = True)
            
            
            
        ############################################################################
        
    # if "filtered_spectra" in st.session_state:
    #     st.success(f"Filtered spectra from file \"{st.session_state['filtered_spectra']['name']}\" were successfully loaded!")
    #     with st.expander("Display filtering parameters:"):
    #         scans_from_protein_str = ""
    #         if st.session_state["filtered_spectra"]["filter_params"]["selected_protein"] is None:
    #             scans_from_protein_str = None
    #         else:
    #             if len(st.session_state["filtered_spectra"]["filter_params"]["scans_from_protein"]) > 10:
    #                 scans_from_protein_str = "\"omitted due to size (>10)\""
    #             else:
    #                 scans_from_protein_str = "[" + ", ".join([str(x) for x in st.session_state["filtered_spectra"]["filter_params"]["scans_from_protein"]]) + "]"

    #         scans_from_peptide_str = ""
    #         if st.session_state["filtered_spectra"]["filter_params"]["selected_peptide"] is None:
    #             scans_from_peptide_str = None
    #         else:
    #             if len(st.session_state["filtered_spectra"]["filter_params"]["scans_from_peptide"]) > 10:
    #                 scans_from_peptide_str = "\"omitted due to size (>10)\""
    #             else:
    #                 scans_from_peptide_str = "[" + ", ".join([str(x) for x in st.session_state["filtered_spectra"]["filter_params"]["scans_from_peptide"]]) + "]"

    #         params_str = "\t{\"source_filename\": " + f"{st.session_state['filtered_spectra']['name']}\n" + \
    #                      "\t \"filter_params\": \n\t\t{\"first_scan\": " + f"{st.session_state['filtered_spectra']['filter_params']['first_scan']}\n" + \
    #                      "\t\t \"last_scan\": " + f"{st.session_state['filtered_spectra']['filter_params']['last_scan']}\n" + \
    #                      "\t\t \"min_mz\": " + f"{st.session_state['filtered_spectra']['filter_params']['min_mz']}\n" + \
    #                      "\t\t \"max_mz\": " + f"{st.session_state['filtered_spectra']['filter_params']['max_mz']}\n" + \
    #                      "\t\t \"min_rt\": " + f"{st.session_state['filtered_spectra']['filter_params']['min_rt']}\n" + \
    #                      "\t\t \"max_rt\": " + f"{st.session_state['filtered_spectra']['filter_params']['max_rt']}\n" + \
    #                      "\t\t \"max_charge\": " + f"{st.session_state['filtered_spectra']['filter_params']['max_charge']}\n" + \
    #                      "\t\t \"max_isotope\": " + f"{st.session_state['filtered_spectra']['filter_params']['max_isotope']}\n" + \
    #                      "\t\t \"selected_protein\": " + f"{st.session_state['filtered_spectra']['filter_params']['selected_protein']}\n" + \
    #                      "\t\t \"scans_from_protein\": " + f"{scans_from_protein_str}\n" + \
    #                      "\t\t \"selected_peptide\": " + f"{st.session_state['filtered_spectra']['filter_params']['selected_peptide']}\n" + \
    #                      "\t\t \"scans_from_peptide\": " + f"{scans_from_peptide_str}\n" + \
    #                      "\t\t}\n\t}"

    #         filter_dl_desc_text1t4 = st.markdown("**Filter Parameters:**")
    #         filter_dl_desc_text1t4 = st.text(params_str)

    ############################################################################
    
    if "consensus_spectrum" in st.session_state:
        
        
        fraggraph_params_header = st.subheader("Fraggraph Parameters", divider = DIV_COLOR)

        fraggraph_params_text = st.markdown("Please specify the parameters used for running Fraggraph.")


        fg_peptidoform1 = st.text_input("Specify a peptidoform to consider:",
                                        value = "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALRE",
                                        help = "Specify a peptidoform to consider for generating the fragment graph, e.g. " +
                                               "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALRE.",
                                        placeholder = "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALRE")

        fg_peptidoform2 = st.text_input("Optionally, specify a peptidoform to compare to:",
                                        value = None,
                                        help = "Specify a peptidoform to compare the fragment graph of the first peptidoform to, e.g. " +
                                               "ARTKQTARKSTGGKAPRKQLATKAARKSAPAT[-79.966331]GGV[+79.966331]KKPHRYRPGTVALRE.",
                                        placeholder = "ARTKQTARKSTGGKAPRKQLATKAARKSAPAT[-79.966331]GGV[+79.966331]KKPHRYRPGTVALRE")

        
        fg_run_l, fg_run_center, fg_run_r = st.columns(3)

        with fg_run_center:
            run_fraggraph = st.button("Run Fraggraph!", use_container_width = True)

        if run_fraggraph:
            fraggraph_main({"mzd": fg_mzd,
                            "cov": fg_cov,
                            "pep1": fg_peptidoform1,
                            "pep2": fg_peptidoform2})
        
    else:
        st.error("No consensus spectrum available. Please check spectra selection and create a consensus spectrum first", icon = "🚨")


























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
