#!/usr/bin/env python3

import pandas as pd

import streamlit as st

from util.tab2.filter import filter_dataframes

# don't use star imports here!
from util.tab2.plots import common_type_hist
from util.tab2.plots import common_type_pie
from util.tab2.plots import mz_dist_ion_type
from util.tab2.plots import rel_ion_intens_perc
from util.tab2.plots import rel_ion_intens_ridge
from util.tab2.plots import per_spec_ion_type
from util.tab2.plots import per_spec_ion_intens
from util.tab2.plots import logo_of_fraction

from util.constants import DIV_COLOR

def main(argv = None) -> None:

    ############################################################################
    header = st.subheader("Statistics Tab", divider = DIV_COLOR)
    tab2_desc = st.markdown("Description text of tab2.")

    ############################################################################
    file_upload_header_fs_csvs = st.subheader("Data Import", divider = DIV_COLOR)

    if "dataframes" not in st.session_state:
        st.info("No results loaded! Please upload them here!")

        fragments_csv = st.file_uploader("Upload the fragment-centric .csv output of the \"Annotation\" tab:",
                                         type = ["csv"],
                                         help = "Upload the fragment-centric .csv output of Fragannot in the \"Annotation\" tab.")

        spectra_csv = st.file_uploader("Upload the spectrum-centric .csv output of the \"Annotation\" tab:",
                                       type = ["csv"],
                                       help = "Upload the spectrum-centric .csv output of Fragannot in the \"Annotation\" tab.")

        if fragments_csv is not None and spectra_csv is not None:
            fragments_csv_df = pd.read_csv(fragments_csv)
            spectra_csv_df = pd.read_csv(spectra_csv)
            if "frag_code" in list(fragments_csv_df.columns) and "perc_internal" in list(spectra_csv_df.columns):
                st.session_state["dataframes"] = [fragments_csv_df, spectra_csv_df]
                st.session_state["dataframes_source"] = {"spectrum_file": None,
                                                         "identifications_file": None,
                                                         "fragment_centric_csv": fragments_csv.name,
                                                         "spectrum_centric_csv": spectra_csv.name}
                st.rerun()
            else:
                st.error("Uploaded .csv files are not in the right format! Did you switch up fragment- and spectra-centric .csv files?", icon = "🚨")

    if "dataframes" in st.session_state:
        if st.session_state["dataframes"] is not None:

            if st.session_state["dataframes_source"]["fragment_centric_csv"] is not None:
                loaded_from_which_files = "uploaded " + str(st.session_state["dataframes_source"]["fragment_centric_csv"]) + \
                                          ", " + str(st.session_state["dataframes_source"]["spectrum_centric_csv"])
            else:
                loaded_from_which_files = "Fragannot run with files " + str(st.session_state["dataframes_source"]["spectrum_file"]) + \
                                          ", " + str(st.session_state["dataframes_source"]["identifications_file"])
            st.success(f"Loaded results from {loaded_from_which_files}!")

    ############################################################################
            show_previews_header = st.subheader("Results Viewer", divider = DIV_COLOR)
            frag_preview_text = st.markdown("Preview of the fragment-centric dataframe:")
            frag_preview = st.dataframe(st.session_state["dataframes"][0], height = 200)
            spec_preview_text = st.markdown("Preview of the spectrum-centric dataframe:")
            spec_preview = st.dataframe(st.session_state["dataframes"][1], height = 200)

    ############################################################################
            filter_dfs_header = st.subheader("Filter Results", divider = DIV_COLOR)
            filter_dfs_desc_text = st.markdown("Results can be further filtered with the filtering options below (optional).")

            # START filtering options
            start_seq_length, end_seq_length = st.select_slider("Peptide sequence lenght:",
                                                            options = range(0, 1001),
                                                            value = (0, 1000),
                                                            help = "For \"Spectrum-centric Statistics\" only spectra with peptides withing range are considered.")

            start_frag_len, end_frag_len = st.select_slider("Fragment ion sequence lenght:",
                                                            options = range(0, 1001),
                                                            value = (0, 1000),
                                                            help = "For \"Fragment-centric Statistics\" only fragments within in range are considered.")

            start_mz, end_mz = st.select_slider("Select m/z range for internal fragments:",
                                                options = range(0, 10001),
                                                value = (0, 10000),
                                                help = "For \"Fragment-centric Statistics\" only fragments within in range are considered.")

            start_int, end_int = st.select_slider("Select intensity range for internal fragments:",
                                                  options = range(0, 1000001),
                                                  value =(0, 1000000),
                                                  help = "For \"Fragment-centric Statistics\" only fragments within in range are considered.")

            ion_selection_text = st.markdown("Select which ions to analyse:")

            ions_col1, ions_col2, ions_col3 = st.columns(3)

            with ions_col1:

                ions_checkbox_cterm = st.markdown("**C-terminal ions:**")

                X_ion = st.checkbox("X ions", value = True)
                Y_ion = st.checkbox("Y ions", value = True)
                Z_ion = st.checkbox("Z ions", value = True)
                Zdot_ion = st.checkbox("Zdot ions", value = True)
                Zp1_ion = st.checkbox("Z+1 ions", value = True)
                Zp2_ion = st.checkbox("Z+2 ions", value = True)
                Zp3_ion = st.checkbox("Z+3 ions", value = True)

            with ions_col2:

                ions_checkbox_nterm = st.markdown("**N-terminal ions:**")

                A_ion = st.checkbox("A ions", value = True)
                B_ion = st.checkbox("B ions", value = True)
                C_ion = st.checkbox("C ions", value = True)
                Cdot_ion = st.checkbox("Cdot ions", value = True)
                Cm1_ion = st.checkbox("C-1 ions", value = True)
                Cp1_ion = st.checkbox("C+1 ions", value = True)

            with ions_col3:

                ions_checkbox_notannotated = st.markdown("**Non-annotated ions:**")

                N_ion = st.checkbox("Non-annotated ions", value = True)

            ion_filter_param = [N_ion, A_ion, B_ion, C_ion, Cdot_ion, Cm1_ion, Cp1_ion, X_ion, Y_ion, Z_ion, Zdot_ion, Zp1_ion, Zp2_ion, Zp3_ion]
            # END filtering options

            filtered_fragments_df, filtered_spectra_df = filter_dataframes(st.session_state["dataframes"],
                                                                           start_seq_length,
                                                                           end_seq_length,
                                                                           start_frag_len,
                                                                           end_frag_len,
                                                                           start_mz,
                                                                           end_mz,
                                                                           start_int,
                                                                           end_int,
                                                                           ion_filter_param)

    ############################################################################
            frag_center_stats_header = st.subheader("Fragment-centric Statistics", divider = DIV_COLOR)
            frag_center_stats_desc = st.markdown("Fragment-centric stats description.")

            frag_center_plot_col1_1, frag_center_plot_col1_2 = st.columns(2)

            with frag_center_plot_col1_1:
                plot1_title = st.markdown("**Histogram of ion types:**")
                plot1 = st.plotly_chart(common_type_hist(filtered_fragments_df), use_container_width = True)
                plot1_caption = st.markdown(f"**Figure 1:** caption.")

            with frag_center_plot_col1_2:
                plot2_title = st.markdown("**Pie chart of ion types:**")
                plot2 = st.plotly_chart(common_type_pie(filtered_fragments_df), use_container_width = True)
                plot2_caption = st.markdown(f"**Figure 2:** caption.")

            plot3_title = st.markdown("**Histogram of m/z per ion type:**")
            plot3 = st.plotly_chart(mz_dist_ion_type(filtered_fragments_df), use_container_width = True)
            plot3_caption = st.markdown(f"**Figure 3:** caption.")

            frag_center_plot_col3_1, frag_center_plot_col3_2 = st.columns(2)

            with frag_center_plot_col3_1:
                plot4_title = st.markdown("**Log Intensities Distribution:**")
                plot4 = st.plotly_chart(rel_ion_intens_perc(filtered_fragments_df), use_container_width = True)
                plot4_caption = st.markdown(f"**Figure 4:** caption.")

            with frag_center_plot_col3_2:
                plot5_title = st.markdown("**Relative Log Intensities:**")
                plot5 = st.plotly_chart(rel_ion_intens_ridge(filtered_fragments_df), use_container_width = True)
                plot5_caption = st.markdown(f"**Figure 5:** caption.")

    ############################################################################
            spec_center_stats_header = st.subheader("Spectrum-centric Statistics", divider = DIV_COLOR)
            spec_center_stats_desc = st.markdown("Spectrum-centric stats description.")

            spec_center_plot_col1_1, spec_center_plot_col1_2 = st.columns(2)

            with spec_center_plot_col1_1:
                plot6_title = st.markdown("**Ion type per spectrum:**")
                plot6 = st.plotly_chart(per_spec_ion_type(filtered_spectra_df), use_container_width = True)
                plot6_caption = st.markdown(f"**Figure 6:** caption.")

            with spec_center_plot_col1_2:
                plot7_title = st.markdown("**Log Intensities:**")
                plot7 = st.plotly_chart(per_spec_ion_intens(filtered_spectra_df), use_container_width = True)
                plot7_caption = st.markdown(f"**Figure 7:** caption.")

            plot8_title = st.markdown("**Logo view of internal fragments:**")
            with st.expander("Expand for changing logo view options."):
                plot8_text_1 = st.markdown("Select number of top spectra with the highest number of internal ions to inlcude in logo. " +
                                           "Choose 0 for making a logo of all spectra.")
                plot8_topn = st.number_input("Choose number of spectra in logo:",
                                             min_value = 0,
                                             max_value = 3,
                                             step = 1)
                p8_min_length_logo, p8_max_length_logo = st.select_slider("Peptide sequence lenght:",
                                                                          options = range(0, 21),
                                                                          value = (0, 20))
            plot8 = st.pyplot(logo_of_fraction(filtered_spectra_df,
                                               filtered_fragments_df,
                                               plot8_topn,
                                               p8_max_length_logo,
                                               p8_min_length_logo))
            plot8_caption = st.markdown(f"**Figure 8:** caption.")
