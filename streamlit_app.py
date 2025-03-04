#!/usr/bin/env python3

"""
#####################################################
##                                                 ##
##            -- STREAMLIT MAIN APP --             ##
##                                                 ##
#####################################################
"""

__version__ = "0.0.3"

import pandas as pd
import streamlit as st
import logging

# import tabs
from tab1 import main as tab1_main
from tab2 import main as tab2_main
from tab3 import main as tab3_main
from tab4 import main as tab4_main

# import constants
from util.constants import REPO_NAME
from util.constants import DIV_COLOR


# main page content
def main_page():

    st.title("Internal Ions Explorer")

    general_description = \
    """
    Description text.
    """
    st.markdown(general_description)

    # set tab names here
    tab1, tab2, tab3, tab4 = st.tabs(["Data Import & Annotation", "Statistics", "Fraggraph", "Documentation"])

    div = \
    """
    #####################################################
    ##                                                 ##
    ##                   -- TAB 1 --                   ##
    ##                                                 ##
    #####################################################
    """

    with tab1:
        tab1_main()

    div = \
    """
    #####################################################
    ##                                                 ##
    ##                   -- TAB 2 --                   ##
    ##                                                 ##
    #####################################################
    """

    with tab2:
        tab2_main()

    div =\
    """
    #####################################################
    ##                                                 ##
    ##                   -- TAB 3 --                   ##
    ##                                                 ##
    #####################################################
    """

    with tab3:
        tab3_main()

    div =\
    """
    #####################################################
    ##                                                 ##
    ##                   -- TAB 3 --                   ##
    ##                                                 ##
    #####################################################
    """

    with tab4:
        tab4_main()


# side bar and main page loader
def main(argv=None) -> None:

    about_str = \
    """
    Description text.
    """ + \
    f"\nThe server is running Fragment Explorer version {__version__}!"

    st.set_page_config(page_title="Fragannot",
                       page_icon=":test_tube:",
                       layout="wide",
                       initial_sidebar_state="expanded",
                       menu_items={"Get Help": f"https://github.com/michabirklbauer/{REPO_NAME}/discussions",
                                   "Report a bug": f"https://github.com/michabirklbauer/{REPO_NAME}/issues",
                                   "About": about_str}
                       )

    with st.sidebar:

        st.title("Internal Ions")

        st.image("img/logo_1.png", caption="Logo")

        st.markdown(about_str)

    ############################################################################
        st.subheader("Parameters", divider=DIV_COLOR)
        st.markdown("These parameters will be used globally across the internal ions explorer.")

        # mgf parsing
        st.markdown("**MGF Parsing**")
        st.text_input("Regex pattern used for parsing scan numbers of spectra from the MGF file.",
                      key="mgf_parser_pattern",
                      value="\\.\\d+\\.",
                      help="The regex pattern for parsing scan numbers of spectra from the MGF file.")

        # tolerance
        st.markdown("**Tolerance**")
        st.number_input("Tolerance in Da:",
                        key="tolerance",
                        value=0.02,
                        help="Fragment ion mass tolerance in Dalton.")

        # ions
        st.markdown("**Ions**")
        fions_nterm_choices = ["a", "b", "c", "cdot", "c-1", "c+1"]
        fions_nterm_mapping = {"a": "A ions", "b": "B ions", "c": "C ions",
                               "cdot": "Cdot ions", "c-1": "C-1 ions", "c+1": "C+1 ions"}
        st.multiselect("Select which n-terminal ions to consider:",
                       options=fions_nterm_choices,
                       default="b",
                       format_func=lambda x: fions_nterm_mapping[x],
                       key="selected_ions_nterm",
                       help="The c-terminal ions considered by Fragannot and Fraggraph.")
        fions_cterm_choices = ["x", "y", "zdot", "z+1", "z+2", "z+3"]
        fions_cterm_mapping = {"x": "X ions", "y": "Y ions", "zdot": "Zdot ions",
                               "z+1": "Z+1 ions", "z+2": "Z+2 ions", "z+3": "Z+3 ions"}
        st.multiselect("Select which c-terminal ions to consider:",
                       options=fions_cterm_choices,
                       default="y",
                       format_func=lambda x: fions_cterm_mapping[x],
                       key="selected_ions_cterm",
                       help="The c-terminal ions considered by Fragannot and Fraggraph.")
        st.markdown("Here are some common selections for HCD/ETD/etc...")
        st.dataframe(pd.DataFrame({"Method": ["HCD", "ETD"], "Ions": ["b, y", "c, z"]}),
                     hide_index=True,
                     use_container_width=True)

        # Non deconvoluted spectra
        st.markdown("----")
        st.markdown("**Annotation of non-deconvoluted spectra**")

        deconvoluted_spectra = st.checkbox("Deconvoluted spectra",
                                           key="deconvoluted_spectra",
                                           value=False,
                                           help="If spectra are already deconvoluted.")

        if not deconvoluted_spectra:
            # Add additional parameters here
            max_charge_auto = st.checkbox("Auto-select Max Charge",
                                          key="max_charge_auto",
                                          value=False,
                                          help="Automatically determine the maximum charge state of the precursor to consider.")
            if not max_charge_auto:
                st.number_input("Max Charge:",
                                key="max_charge",
                                value=3,
                                format="%d",
                                help="Maximum charge state of the precursor to consider.")

            max_isotope_auto = st.checkbox("Auto-select Max Isotope",
                                           key="max_isotope_auto",
                                           value=False,
                                           help="Automatically determine the maximum isotope to consider.")
            if not max_isotope_auto:
                st.number_input("Max Isotope:",
                                key="max_isotope",
                                value=5,
                                format="%d",
                                help="Maximum isotope to consider.")

            st.checkbox("Charge Reduction",
                        key="charge_reduction",
                        value=False,
                        help="Charge reduction implies that a charge is lost upon fragmentation event. This is typically the case for electron-based fragmentation (e.g. ETD ECD).")

    ############################################################################
        st.subheader("About the Project", divider=DIV_COLOR)

        contact_str = \
            "**Contact:**\n- [Arthur Grimaud](mailto:agrimaud@bmb.sdu.dk)\n- [Veit Schwämmle](mailto:veits@bmb.sdu.dk)\n" + \
            "- [Caroline Lennartsson](mailto:caroline.lennartsson@cpr.ku.dk)\n- [Louise Buur](mailto:louise.buur@fh-hagenberg.at)\n" + \
            "- [Micha Birklbauer](mailto:micha.birklbauer@gmail.com)\n- [Vladimir Gorshkov](mailto:homer2k@gmail.com)\n" + \
            "- [Zoltan Udvardy](mailto:zoltan.udvardy.ipbs@gmail.com)"
        st.markdown(contact_str)

        license_str = "**License:** [???]()"
        st.markdown(license_str)

        project_str = f"**Project Page:** [GitHub](https://github.com/michabirklbauer/{REPO_NAME}/)"
        st.markdown(project_str)

    logging.basicConfig(level=logging.INFO)

    main_page()


if __name__ == "__main__":
    main()
