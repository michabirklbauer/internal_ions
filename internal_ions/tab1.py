import json

import streamlit as st

from .fragannot.fragannot_call import fragannot_call

from .util.converter import JSONConverter
from .util.redirect import st_stdout
from .util.spectrumio import read_spectrum_file
from .util.psmio import read_identifications, read_id_file
from .util.streamlit_utils import dataframe_to_csv_stream

from .util.constants import DIV_COLOR, SUPPORTED_FILETYPES


def main(argv=None) -> None:

    ############################################################################
    st.subheader("Annotation Tab", divider=DIV_COLOR)
    st.markdown("Upload data and set parameters for data annotation.")

    ############################################################################
    st.subheader("Data Import", divider=DIV_COLOR)

    uploaded_spectrum_file = st.file_uploader(
        "Upload a spectrum file (we assume all spectra are MS2-level):",
        key="uploaded_spectrum_file",
        type=["mgf"],
        help="Upload a spectrum file to be analyzed in .mgf format.",
        disabled=st.session_state["sidebar_disabled"])

    spectrum_file = None
    if uploaded_spectrum_file is not None:
        with st.status("Reading spectra...") as spectra_reading_status:
            with st_stdout("info"):
                spectrum_file = read_spectrum_file(uploaded_spectrum_file)
            st.success("Read all spectra successfully!")
            spectra_reading_status.update(label=f"Read all spectra from file {uploaded_spectrum_file.name} successfully!", state="complete")

    identifications_file = st.file_uploader("Upload an identification file:",
                                            key="identifications_file",
                                            type=None,  # ["mzid"],
                                            help="Upload a identification file that contains PSMs of the spectrum file in .mzid format.",
                                            disabled=st.session_state["sidebar_disabled"])

    identifications_file_format = st.selectbox("Select the file format of the identifications file:",
                                               key="identifications_file_format",
                                               options=['infer'] + SUPPORTED_FILETYPES,
                                               index=0,
                                               placeholder="None selected!",
                                               help="Select the file format of the identifications file, supported options are based on psm_utils.",
                                               disabled=st.session_state["sidebar_disabled"])

    if identifications_file is not None and identifications_file_format is not None:
        if spectrum_file is not None:
            try:
                with st.status("Reading identifications...") as identifications_reading_status:
                    psm_list = read_id_file(identifications_file, identifications_file_format)
                    with st_stdout("info"):
                        if "identifications" not in st.session_state:
                            st.session_state["identifications"] = read_identifications(psm_list, identifications_file.name, spectrum_file)
                    st.success("Read all identifications successfully!")
                    identifications_reading_status.update(label=f"Read all identifications from file {st.session_state.identifications_file.name} successfully!", state="complete")
            except Exception as e:
                st.warning("Something went wrong! If 'file format = infer' was selected please manually specify file format!", icon="⚠️")
                with st.expander("Show exception"):
                    st.exception(e)
        else:
            st.warning("Please upload a spectrum file before the identifications file!", icon="⚠️")

    st.session_state["fragannot_call_ion_selection"] = st.session_state["selected_ions_nterm"] + st.session_state["selected_ions_cterm"]

    charges_str = st.text_input("Charges to consider [comma delimited]:",
                                value="+1",
                                help="The charges to consider for fragment ions. Multiple entries should be delimited by commas!",
                                disabled=st.session_state["sidebar_disabled"])
    st.session_state["charges"] = [charge.strip() for charge in charges_str.split(",")]

    losses_str = st.text_input("Neutral losses to consider [comma delimited]",
                               value="H2O",
                               help="Neutral losses to consider for fragment ions. Multiple entries should be delimited by commas!",
                               disabled=st.session_state["sidebar_disabled"])
    st.session_state["losses"] = [loss.strip() for loss in losses_str.split(",")]

    st.checkbox("Deisotope spectra", key="deisotope", value=True, help="Deisotope uploaded spectra or not.", disabled=st.session_state["sidebar_disabled"])

    ############################################################################
    st.subheader("Annotation", divider=DIV_COLOR)

    l1, center_button, r1 = st.columns(3)

    with center_button:
        if st.session_state["sidebar_disabled"]:
            enable_sidebar = st.button("Re-enable upload and parameter selection!",
                                       type="secondary",
                                       use_container_width=True)
            if enable_sidebar:
                st.session_state["sidebar_disabled"] = False
                st.rerun()
        run_analysis = st.button("Load files and run Fragannot!",
                                 type="primary",
                                 use_container_width=True)

    if run_analysis:
        st.session_state["sidebar_disabled"] = True
        st.session_state["run_fragannot"] = True
        st.rerun()

    if "run_fragannot" in st.session_state and st.session_state["run_fragannot"]:
        cond1 = spectrum_file is not None
        cond2 = st.session_state.identifications_file is not None
        cond3 = st.session_state.identifications_file_format is not None
        if cond1 and cond2 and cond3:
            st.markdown("**Parameters:**")
            run_info_str = f"\tSpectrum file name: {spectrum_file.name}\n" + \
                           f"\tIdentifications file name: {st.session_state.identifications_file.name}\n" + \
                           f"\tTolerance: {st.session_state.tolerance}\n" + \
                           f"\tSelected ions: {', '.join(st.session_state['fragannot_call_ion_selection'])}\n" + \
                           f"\tCharges: {', '.join(st.session_state['charges'])}\n" + \
                           f"\tLosses: {', '.join(st.session_state['losses'])}\n" + \
                           f"\tDeisotope: {st.session_state.deisotope}"
            st.text(run_info_str)
            with st.status("Fragannot is running! Show logging info:") as st_status:
                with st_stdout("info"):
                    try:
                        result = fragannot_call(spectrum_file,
                                                psm_list,
                                                float(st.session_state.tolerance),
                                                st.session_state["fragannot_call_ion_selection"],
                                                st.session_state["charges"],
                                                st.session_state["losses"],
                                                st.session_state.deisotope)

                        converter = JSONConverter()
                        st.session_state["result"] = result
                        st.session_state["dataframes"] = converter.to_dataframes(data=result)
                        st.session_state["dataframes_source"] = {"spectrum_file": spectrum_file.name,
                                                                 "identifications_file": st.session_state.identifications_file.name,
                                                                 "fragment_centric_csv": None,
                                                                 "spectrum_centric_csv": None}
                        status_1 = 0
                        st_status.update(label="Fragannot finished successfully!", state="complete")
                    except Exception as e:
                        st.exception(e)
                        status_1 = 1
            if status_1 == 0:
                st.success("Fragannot finished successfully!")
                st.session_state["run_fragannot"] = False
            else:
                st.error("Fragannot stopped prematurely! See log for more information!")
                st.session_state["run_fragannot"] = False
        else:
            st.error("You need to specify a spectrum AND identifications file AND the file format!")
            st.session_state["run_fragannot"] = False

    ############################################################################
    if "dataframes" in st.session_state:
        st.subheader("Results Preview", divider=DIV_COLOR)
        st.markdown("Fragment-centric")
        st.dataframe(st.session_state["dataframes"][0].head(10), use_container_width=True)
        st.markdown("Spectrum-centric")
        st.dataframe(st.session_state["dataframes"][1].head(10), use_container_width=True)

        st.subheader("Download Results", divider=DIV_COLOR)

        dl_l1, dl_center, dl_r1 = st.columns(3)

        with dl_l1:
            st.download_button(label="Download Fragment-centric data!",
                               data=dataframe_to_csv_stream(st.session_state["dataframes"][0]),
                               file_name="fragment_centric.csv",
                               mime="text/csv",
                               help="Download fragment-centric Fragannot results in .csv format.",
                               type="primary",
                               use_container_width=True)
        with dl_center:
            st.download_button(label="Download Spectrum-centric data!",
                               data=dataframe_to_csv_stream(st.session_state["dataframes"][1]),
                               file_name="spectrum_centric.csv",
                               mime="text/csv",
                               help="Download spectrum-centric Fragannot results in .csv format.",
                               type="primary",
                               use_container_width=True)

        if "result" in st.session_state:
            with dl_r1:
                st.download_button(label="Download raw result in .json format!",
                                   data=json.dumps(st.session_state["result"]),
                                   file_name="result.json",
                                   mime="text/json",
                                   help="Download raw Fragannot results in .json file format.",
                                   type="primary",
                                   use_container_width=True)
