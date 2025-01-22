#!/usr/bin/env python3

import io
import json
import numpy as np
import pandas as pd

import streamlit as st

from fragannot.fragannot_call import fragannot_call

from util.converter import JSONConverter
from util.redirect import st_stdout
from util.spectrumio import read_spectra
from util.psmio import read_identifications
from util.streamlit_utils import dataframe_to_csv_stream

from util.constants import FRAGANNOT_ION_NAMES
from util.constants import DIV_COLOR
from util.constants import SUPPORTED_FILETYPES

def reset_spectra() -> None:
    st.session_state["rerun_spectra_reading"] = True

def reset_identifications() -> None:
    st.session_state["rerun_identifications_reading"] = True

def main(argv = None) -> None:

    ############################################################################
    header = st.subheader("Annotation Tab", divider = DIV_COLOR)
    tab1_desc = st.markdown("Description text of tab1.")

    ############################################################################
    data_import_tab1_header = st.subheader("Data Import", divider = DIV_COLOR)

    spectrum_file = st.file_uploader("Upload a spectrum file (we assume all spectra are MS2-level):",
                                     key = "spectrum_file",
                                     type = ["mgf"],
                                     on_change = reset_spectra,
                                     help = "Upload a spectrum file to be analyzed in .mgf format.")

    if spectrum_file is not None:
        with st.status("Reading spectra...") as spectra_reading_status:
            with st_stdout("info"):
                if "spectra" not in st.session_state:
                    st.session_state["spectra"] = read_spectra(spectrum_file,
                                                               spectrum_file.name,
                                                               st.session_state["mgf_parser_pattern"])
                    st.session_state["rerun_spectra_reading"] = False
                if "spectra" in st.session_state:
                    if st.session_state["rerun_spectra_reading"]:
                        st.session_state["spectra"] = read_spectra(spectrum_file,
                                                                   spectrum_file.name,
                                                                   st.session_state["mgf_parser_pattern"])
                        st.session_state["rerun_spectra_reading"] = False
            read_spectra_successfully = st.success("Read all spectra successfully!")
            spectra_reading_status.update(label = f"Read all spectra from file {st.session_state.spectrum_file.name} successfully!", state = "complete")


    identifications_file = st.file_uploader("Upload an identification file:",
                                            key = "identifications_file",
                                            type = None, #["mzid"],
                                            on_change = reset_identifications,
                                            help = "Upload a identification file that contains PSMs of the spectrum file in .mzid format.")
                                            
    identifications_file_format = st.selectbox("Select the file format of the identifications file:",
                                               key = "identifications_file_format",
                                               options = SUPPORTED_FILETYPES,
                                               index = None,
                                               on_change = reset_identifications,
                                               placeholder = "None selected!",
                                               help = "Select the file format of the identifications file, supported options are based on psm_utils.")

    if identifications_file is not None and identifications_file_format is not None:
        with st.status("Reading identifications...") as identifications_reading_status:
            with st_stdout("info"):
                if "identifications" not in st.session_state:
                    st.session_state["identifications"] = read_identifications(identifications_file,
                                                                               identifications_file_format,
                                                                               identifications_file.name,
                                                                               st.session_state["mgf_parser_pattern"])
                    st.session_state["rerun_identifications_reading"] = False
                if "identifications" in st.session_state:
                    if st.session_state["rerun_identifications_reading"]:
                        st.session_state["identifications"] = read_identifications(identifications_file,
                                                                                   identifications_file_format,
                                                                                   identifications_file.name,
                                                                                   st.session_state["mgf_parser_pattern"])
                        st.session_state["rerun_identifications_reading"] = False
            read_identifications_successfully = st.success("Read all identifications successfully!")
            identifications_reading_status.update(label = f"Read all identifications from file {st.session_state.identifications_file.name} successfully!", state = "complete")

    st.session_state["fragannot_call_ion_selection"] = st.session_state["selected_ions_nterm"] + st.session_state["selected_ions_cterm"]

    charges_str = st.text_input("Charges to consider [comma delimited]:",
                                value = "-1, +1",
                                help = "The charges to consider for fragment ions. Multiple entries should be delimited by commas!")
    st.session_state["charges"] = [charge.strip() for charge in charges_str.split(",")]

    losses_str = st.text_input("Neutral losses to consider [comma delimited]",
                               value = "H2O",
                               help = "Neutral losses to consider for fragment ions. Multiple entries should be delimited by commas!")
    st.session_state["losses"] = [loss.strip() for loss in losses_str.split(",")]

    deisotope = st.checkbox("Deisotope spectra",
                            key = "deisotope",
                            value = True,
                            help = "Deisotope uploaded spectra or not.")

    ############################################################################
    st.subheader("Annotation", divider = DIV_COLOR)

    l1, center_button, r1 = st.columns(3)

    with center_button:
        run_analysis = st.button("Load files and run Fragannot!",
                                 type = "primary",
                                 use_container_width = True)

    if run_analysis:
        if st.session_state.spectrum_file is not None and st.session_state.identifications_file is not None:
            run_info_title = st.markdown("**Parameters:**")
            run_info_str = f"\tSpectrum file name: {st.session_state.spectrum_file.name}\n" + \
                           f"\tIdentifications file name: {st.session_state.identifications_file.name}\n" + \
                           f"\tTolerance: {st.session_state.tolerance}\n" + \
                           f"\tSelected ions: {', '.join(st.session_state['fragannot_call_ion_selection'])}\n" + \
                           f"\tCharges: {', '.join(st.session_state['charges'])}\n" + \
                           f"\tLosses: {', '.join(st.session_state['losses'])}\n" + \
                           f"\tDeisotope: {st.session_state.deisotope}"
            run_info = st.text(run_info_str)
            with st.status("Fragannot is running! Show logging info:") as st_status:
                with st_stdout("info"):
                    try:
                        result = fragannot_call(st.session_state.spectrum_file,
                                                st.session_state.identifications_file,
                                                float(st.session_state.tolerance),
                                                st.session_state["fragannot_call_ion_selection"],
                                                st.session_state["charges"],
                                                st.session_state["losses"],
                                                st.session_state.deisotope)
                        converter = JSONConverter()
                        st.session_state["result"] = result
                        st.session_state["dataframes"] = converter.to_dataframes(data = result)
                        st.session_state["dataframes_source"] = {"spectrum_file": st.session_state.spectrum_file.name,
                                                                 "identifications_file": st.session_state.identifications_file.name,
                                                                 "fragment_centric_csv": None,
                                                                 "spectrum_centric_csv": None}
                        status_1 = 0
                        st_status.update(label = "Fragannot finished successfully!", state = "complete")
                    except Exception as e:
                        this_e = st.exception(e)
                        status_1 = 1
            if status_1 == 0:
                res_status_1 = st.success("Fragannot finished successfully!")
            else:
                res_status_1 = st.error("Fragannot stopped prematurely! See log for more information!")
        else:
            res_status_1 = st.error("You need to specify a spectrum AND identifications file!")

    ############################################################################
    if "dataframes" in st.session_state:
        results_preview_header = st.subheader("Results Preview", divider = DIV_COLOR)
        preview_csv_1_desc = st.markdown("Fragment-centric")
        preview_csv_1 = st.dataframe(st.session_state["dataframes"][0].head(10), use_container_width = True)
        preview_csv_1_desc = st.markdown("Spectrum-centric")
        preview_csv_1 = st.dataframe(st.session_state["dataframes"][1].head(10), use_container_width = True)

        results_header = st.subheader("Download Results", divider = DIV_COLOR)

        dl_l1, dl_center, dl_r1 = st.columns(3)

        with dl_l1:
            csv_1 = st.download_button(label = "Download Fragment-centric data!",
                                       data = dataframe_to_csv_stream(st.session_state["dataframes"][0]),
                                       file_name = "fragment_centric.csv",
                                       mime = "text/csv",
                                       help = "Download fragment-centric Fragannot results in .csv format.",
                                       type = "primary",
                                       use_container_width = True)
        with dl_center:
            csv_2 = st.download_button(label = "Download Spectrum-centric data!",
                                       data = dataframe_to_csv_stream(st.session_state["dataframes"][1]),
                                       file_name = "spectrum_centric.csv",
                                       mime = "text/csv",
                                       help = "Download spectrum-centric Fragannot results in .csv format.",
                                       type = "primary",
                                       use_container_width = True)

        if "result" in st.session_state:
            with dl_r1:
                json_output = st.download_button(label = "Download raw result in .json format!",
                                                 data = json.dumps(st.session_state["result"]),
                                                 file_name = "result.json",
                                                 mime = "text/json",
                                                 help = "Download raw Fragannot results in .json file format.",
                                                 type = "primary",
                                                 use_container_width = True)
