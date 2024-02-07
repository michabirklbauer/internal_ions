#!/usr/bin/env python3

import io
import pandas as pd
import streamlit as st

@st.cache_data
def dataframe_to_csv_stream(dataframe: pd.DataFrame):
    text = dataframe.to_csv(index = False).encode("utf-8")

    return text

@st.cache_data
def dataframe_to_xlsx_stream(dataframe: pd.DataFrame, sheet_name: str):
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine = "openpyxl") as writer:
        dataframe.to_excel(writer, sheet_name = sheet_name)
        writer.close()

    return buffer
