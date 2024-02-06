#!/usr/bin/env python3

import streamlit as st

def main(argv = None) -> None:
    header = st.subheader("header tab 3", divider = "rainbow")

    if "dataframes" in st.session_state:
        if st.session_state["dataframes"] is not None:
            result = st.dataframe(st.session_state["dataframes"][0])
