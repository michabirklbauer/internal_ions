#!/usr/bin/env python3

import streamlit as st

from util.constants import DIV_COLOR

def main(argv = None) -> None:
    header = st.subheader("header tab 4", divider = DIV_COLOR)

    if "dataframes" in st.session_state:
        if st.session_state["dataframes"] is not None:
            result = st.dataframe(st.session_state["dataframes"][0])
