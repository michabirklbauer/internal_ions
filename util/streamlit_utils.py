import io
import pandas as pd
from pyteomics import mgf

import streamlit as st

from util.capture import CaptureStdOut

from typing import Any, Dict


@st.cache_data
def dataframe_to_csv_stream(dataframe: pd.DataFrame):
    text = dataframe.to_csv(index=False).encode("utf-8")

    return text


# currently not used
@st.cache_data
def dataframe_to_xlsx_stream(dataframe: pd.DataFrame, sheet_name: str):
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        dataframe.to_excel(writer, sheet_name=sheet_name)
        writer.close()

    return buffer


# currently not used
@st.cache_data
def spectra_to_mgf_stream(_spectra: Dict[str, Any]):
    std_out = []
    with CaptureStdOut(std_out) as std_out:
        # https://pyteomics.readthedocs.io/en/latest/api/mgf.html#pyteomics.mgf.write
        mgf.write(_spectra["spectra"], write_ions=True)
    return "\n".join(std_out)
