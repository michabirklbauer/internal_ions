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
    tab1, tab2, tab3, tab4 = st.tabs(["Annotation", "Statistics", "Spectrum", "Fraggraph"])

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

    div = \
    """
    #####################################################
    ##                                                 ##
    ##                   -- TAB 4 --                   ##
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
                       menu_items = {"Get Help": "https://github.com/michabirklbauer/internal-ions/discussions",
                                     "Report a bug": "https://github.com/michabirklbauer/internal-ions/issues",
                                     "About": about_str}
                       )

    with st.sidebar:

        title = st.title("Internal Ions")

        logo = st.image("img/logo_1.png", caption = "Logo")

        doc = st.markdown(about_str)

    ############################################################################
        input_header = st.subheader("Parameters", divider = DIV_COLOR)
        input_desc = st.markdown("These parameters will be used globally across the internal ions explorer.")

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

    ############################################################################
        contact_header = st.subheader("About the Project", divider = DIV_COLOR)

        contact_str = \
            "**Contact:**\n- [Arthur Grimaud](mailto:agrimaud@bmb.sdu.dk)\n- [Veit Schwämmle](veits@bmb.sdu.dk)\n" + \
            "- [Caroline Lennartsson](mailto:caroline.lennartsson@cpr.ku.dk)\n- [Louise Buur](louise.buur@fh-hagenberg.at)\n" + \
            "- [Micha Birklbauer](mailto:micha.birklbauer@gmail.com)\n- [Vladimir Gorshkov](homer2k@gmail.com)\n" + \
            "- [Zoltan Udvardy](zoltan.udvardy.ipbs@gmail.com)"
        contact = st.markdown(contact_str)

        license_str = "**License:** [???]()"
        license = st.markdown(license_str)

        project_str = "**Project Page:** [GitHub](https://github.com/michabirklbauer/internal-ions/)"
        project = st.markdown(project_str)

    main_page()

if __name__ == "__main__":
    main()
