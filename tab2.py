#!/usr/bin/env python3

import pandas as pd

import streamlit as st

from util.tab2.filter import filter_dataframes
from util.tab2.filter2 import filter_dataframe

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
            show_previews_header = st.subheader("View Result", divider = DIV_COLOR)
            frag_preview_text = st.markdown("Data of the fragment-centric dataframe:")
            N_ion = st.checkbox("Show non-annotated fragments",
                                value=False,
                                help="Whether or not to show non-annotated fragments in the output.")
            if N_ion:
                frag_preview = st.dataframe(st.session_state["dataframes"][0], height = 400, use_container_width = True)
            else:
                frag_preview = st.dataframe(st.session_state["dataframes"][0].dropna(subset = "frag_seq"), height = 400, use_container_width = True)

            spec_preview_text = st.markdown("Data of the spectrum-centric dataframe:")
            spec_preview = st.dataframe(st.session_state["dataframes"][1], height = 400, use_container_width = True)
            modify = st.checkbox("Filter data", value=False, help="Display filter options")
            if modify:
                st.session_state["frag_center_filtered"] = filter_dataframe(st.session_state["dataframes"][0], 'spectrum')
                st.session_state["spec_center_filtered"] = filter_dataframe(st.session_state["dataframes"][1], 'fragment')
            else:
                st.session_state["frag_center_filtered"] = st.session_state["dataframes"][0]
                st.session_state["spec_center_filtered"] = st.session_state["dataframes"][1]

            ion_filter_param = [N_ion,
                                "a" in st.session_state["selected_ions_nterm"],
                                "b" in st.session_state["selected_ions_nterm"],
                                "c" in st.session_state["selected_ions_nterm"],
                                "cdot" in st.session_state["selected_ions_nterm"],
                                "c-1" in st.session_state["selected_ions_nterm"],
                                "c+1" in st.session_state["selected_ions_nterm"],
                                "x" in st.session_state["selected_ions_cterm"],
                                "y" in st.session_state["selected_ions_cterm"],
                                False,  # z ions, this is a relict
                                "zdot" in st.session_state["selected_ions_cterm"],
                                "z+1" in st.session_state["selected_ions_cterm"],
                                "z+2" in st.session_state["selected_ions_cterm"],
                                "z+3" in st.session_state["selected_ions_cterm"]]
            st.session_state["frag_center_filtered"], st.session_state["spec_center_filtered"] = filter_dataframes(
                [st.session_state["frag_center_filtered"], st.session_state["spec_center_filtered"]], ion_filter_param)

    ############################################################################

            frag_center_stats_header = st.subheader("Fragment-centric Statistics", divider = DIV_COLOR)
            frag_center_stats_desc = st.markdown("Fragment-centric stats description.")

            frag_center_plot_col1_1, frag_center_plot_col1_2 = st.columns(2)

            with frag_center_plot_col1_1:
                plot1_title = st.markdown("**Distribution of Ion Types**")
                plot1 = st.plotly_chart(common_type_hist(st.session_state["frag_center_filtered"]), use_container_width = True)
                plot1_caption = st.markdown(f"**Figure 1:** Histogram illustrating the frequency distribution of ion types present in the dataset.")

            with frag_center_plot_col1_2:
                plot2_title = st.markdown("**Ion Type Proportions**")
                plot2 = st.plotly_chart(common_type_pie(st.session_state["frag_center_filtered"]), use_container_width = True)
                plot2_caption = st.markdown(f"**Figure 2:** Pie chart displaying the proportional composition of ion types within the dataset.")

            plot3_title = st.markdown("**Distribution of m/z Values Across Ion Types**")
            plot3 = st.plotly_chart(mz_dist_ion_type(st.session_state["frag_center_filtered"]), use_container_width = True)
            plot3_caption = st.markdown(f"**Figure 3:** Histogram depicting the distribution of mass-to-charge ratio (m/z) values across different ion types.")

            frag_center_plot_col3_1, frag_center_plot_col3_2 = st.columns(2)

            with frag_center_plot_col3_1:
                plot4_title = st.markdown("**Distribution of the Log Intensity**")
                plot4 = st.plotly_chart(rel_ion_intens_perc(st.session_state["frag_center_filtered"]), use_container_width = True)
                plot4_caption = st.markdown(f"**Figure 4:** Distribution of log-transformed intensity values.")

            with frag_center_plot_col3_2:
                plot5_title = st.markdown("**Relative Log Intensities**")
                plot5 = st.plotly_chart(rel_ion_intens_ridge(st.session_state["frag_center_filtered"]), use_container_width = True)
                plot5_caption = st.markdown(f"**Figure 5:** Distribution of log-transformed intensities relative to their respective scales.")

    ############################################################################
            spec_center_stats_header = st.subheader("Spectrum-centric Statistics", divider = DIV_COLOR)
            spec_center_stats_desc = st.markdown("Spectrum-centric stats description.")

            spec_center_plot_col1_1, spec_center_plot_col1_2 = st.columns(2)

            with spec_center_plot_col1_1:
                plot6_title = st.markdown("**Ion Type Distribution Per Spectra**")
                plot6 = st.plotly_chart(per_spec_ion_type(st.session_state["spec_center_filtered"]), use_container_width = True)
                plot6_caption = st.markdown(f"**Figure 6:** Distribution of ion types within each spectrum, providing insights into the diversity and abundance of ions detected across the dataset.")

            with spec_center_plot_col1_2:
                plot7_title = st.markdown("**Log Intensities:**")
                plot7 = st.plotly_chart(per_spec_ion_intens(st.session_state["spec_center_filtered"]), use_container_width = True)
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
            plot8 = st.pyplot(logo_of_fraction(st.session_state["spec_center_filtered"],
                                               st.session_state["frag_center_filtered"],
                                               plot8_topn,
                                               p8_max_length_logo,
                                               p8_min_length_logo))
            plot8_caption = st.markdown(f"**Figure 8:** caption.")
