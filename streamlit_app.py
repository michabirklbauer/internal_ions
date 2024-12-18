#!/usr/bin/env python3

"""
#####################################################
##                                                 ##
##            -- STREAMLIT MAIN APP --             ##
##                                                 ##
#####################################################
"""

import pandas as pd
import streamlit as st

#import tabs
from tab1 import main as tab1_main
from tab2 import main as tab2_main
from tab3 import main as tab3_main
from tab4 import main as tab4_main

# import constants
from util.constants import REPO_NAME
from util.constants import DIV_COLOR

# main page content
def main_page():

    title = st.title("Internal Ions Explorer")

    general_description = \
    """
    Description text.
    """
    description = st.markdown(general_description)

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
def main(argv = None) -> None:

    about_str = \
    """
    Description text.
    """

    st.set_page_config(page_title = "Fragannot",
                       page_icon = ":test_tube:",
                       layout = "wide",
                       initial_sidebar_state = "expanded",
                       menu_items = {"Get Help": f"https://github.com/michabirklbauer/{REPO_NAME}/discussions",
                                     "Report a bug": f"https://github.com/michabirklbauer/{REPO_NAME}/issues",
                                     "About": about_str}
                       )

    with st.sidebar:

        title = st.title("Internal Ions")

        logo = st.image("img/logo_1.png", caption = "Logo")

        doc = st.markdown(about_str)

    ############################################################################
        input_header = st.subheader("Parameters", divider = DIV_COLOR)
        input_desc = st.markdown("These parameters will be used globally across the internal ions explorer.")

        # mgf parsing
        mgf_parser_pattern_desc = st.markdown("**MGF Parsing**")
        mgf_parser_pattern = st.text_input("Regex pattern used for parsing scan numbers of spectra from the MGF file.",
                                           key = "mgf_parser_pattern",
                                           value = "\\.\\d+\\.",
                                           help = "The regex pattern for parsing scan numbers of spectra from the MGF file.")

        # tolerance
        ttolerance_desc = st.markdown("**Tolerance**")
        ttolerance = st.number_input("Tolerance in Da:",
                                     key = "tolerance",
                                     value = 0.02,
                                     help = "Fragment ion mass tolerance in Dalton.")

        # ions
        fions_selectbox_desc = st.markdown("**Ions**")
        fions_nterm_choices = ["a", "b", "c", "cdot", "c-1", "c+1"]
        fions_nterm_mapping = {"a": "A ions", "b": "B ions", "c": "C ions",
                               "cdot": "Cdot ions", "c-1": "C-1 ions", "c+1": "C+1 ions"}
        fions_selectbox_nterm = st.multiselect("Select which n-terminal ions to consider:",
                                               options = fions_nterm_choices,
                                               default = "b",
                                               format_func = lambda x: fions_nterm_mapping[x],
                                               key = "selected_ions_nterm",
                                               help = "The c-terminal ions considered by Fragannot and Fraggraph.")
        fions_cterm_choices = ["x", "y", "zdot", "z+1", "z+2", "z+3"]
        fions_cterm_mapping = {"x": "X ions", "y": "Y ions", "zdot": "Zdot ions",
                               "z+1": "Z+1 ions", "z+2": "Z+2 ions", "z+3": "Z+3 ions"}
        fions_selectbox_cterm = st.multiselect("Select which c-terminal ions to consider:",
                                               options = fions_cterm_choices,
                                               default = "y",
                                               format_func = lambda x: fions_cterm_mapping[x],
                                               key = "selected_ions_cterm",
                                               help = "The c-terminal ions considered by Fragannot and Fraggraph.")
        fions_desc_text = st.markdown("Here are some common selections for HCD/ETD/etc...")
        fions_desc_table = st.dataframe(pd.DataFrame({"Method": ["HCD", "ETD"],
                                                      "Ions": ["b, y", "c, z"]}),
                                        hide_index = True,
                                        use_container_width = True)

        # Non deconvoluted spectra
        st.markdown("----")
        fions_selectbox_desc = st.markdown("**Annotation of non-deconvoluted spectra**")

        deconvoluted_spectra = st.checkbox("Deconvoluted spectra",
                                           key = "deconvoluted_spectra",
                                           value = False,
                                           help = "If spectra are already deconvoluted.")

        if not deconvoluted_spectra:
            # Add additional parameters here
            monoisotopic = False
            max_charge_auto = st.checkbox("Auto-select Max Charge",
                                          key = "max_charge_auto",
                                          value = False,
                                          help = "Automatically determine the maximum charge state of the precursor to consider.")
            if not max_charge_auto:
                max_charge = st.number_input("Max Charge:",
                                             key = "max_charge",
                                             value = 3,
                                             format = "%d",
                                             help = "Maximum charge state of the precursor to consider.")

            max_isotope_auto = st.checkbox("Auto-select Max Isotope",
                                           key = "max_isotope_auto",
                                           value = False,
                                           help = "Automatically determine the maximum isotope to consider.")
            if not max_isotope_auto:
                max_isotope = st.number_input("Max Isotope:",
                                              key = "max_isotope",
                                              value = 5,
                                              format = "%d",
                                              help = "Maximum isotope to consider.")

            charge_reduction = st.checkbox("Charge Reduction",
                                           key = "charge_reduction",
                                           value = False,
                                           help = "Charge reduction implies that a charge is lost upon fragmentation event. This is typically the case for electron-based fragmentation (e.g. ETD ECD).")
        else:
            # Set default parameters here
            max_charge = 1
            monoisotopic = True
            max_isotope_auto = True
            max_isotope = 5
            charge_reduction = False

    ############################################################################
        contact_header = st.subheader("About the Project", divider = DIV_COLOR)

        contact_str = \
            "**Contact:**\n- [Arthur Grimaud](mailto:agrimaud@bmb.sdu.dk)\n- [Veit Schwämmle](mailto:veits@bmb.sdu.dk)\n" + \
            "- [Caroline Lennartsson](mailto:caroline.lennartsson@cpr.ku.dk)\n- [Louise Buur](mailto:louise.buur@fh-hagenberg.at)\n" + \
            "- [Micha Birklbauer](mailto:micha.birklbauer@gmail.com)\n- [Vladimir Gorshkov](mailto:homer2k@gmail.com)\n" + \
            "- [Zoltan Udvardy](mailto:zoltan.udvardy.ipbs@gmail.com)"
        contact = st.markdown(contact_str)

        license_str = "**License:** [???]()"
        license = st.markdown(license_str)

        project_str = f"**Project Page:** [GitHub](https://github.com/michabirklbauer/{REPO_NAME}/)"
        project = st.markdown(project_str)

    main_page()

if __name__ == "__main__":
    main()
