#!/usr/bin/env python3

import streamlit as st

from util.tab3.example_plot import plot_spectrum
from fragannot.spectrumfile import *

def main(argv = None) -> None:
    header = st.subheader("Spectrum view", divider = "rainbow")
    

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
            
            