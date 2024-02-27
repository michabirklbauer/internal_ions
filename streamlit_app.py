#!/usr/bin/env python3

"""
#####################################################
##                                                 ##
##            -- STREAMLIT MAIN APP --             ##
##                                                 ##
#####################################################
"""

import streamlit as st

#import tabs
from tab1 import main as tab1_main
from tab2 import main as tab2_main
from tab3 import main as tab3_main
from tab4 import main as tab4_main

# main page content
def main_page():

    title = st.title("Internal Ions Explorer")

    general_description = \
    """
    Description text.
    """
    description = st.markdown(general_description)

    # set tab names here
    tab1, tab2, tab3, tab4 = st.tabs(["Annotation", "Statistics", "Spectrum", "Tab4"])

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

    title = st.sidebar.title("Internal Ions")

    logo = st.sidebar.image("img/logo_1.png", caption = "Logo")

    doc = st.sidebar.markdown(about_str)

    contact_str = \
        "**Contact:** [Arthur Grimaud](mailto:agrimaud@bmb.sdu.dk), [Veit Schwämmle](veits@bmb.sdu.dk), " + \
        "[Caroline Lennartsson](mailto:caroline.lennartsson@cpr.ku.dk), [Louise Buur](louise.buur@fh-hagenberg.at), " + \
        "[Micha Birklbauer](mailto:micha.birklbauer@gmail.com), [Vladimir Gorshkov](homer2k@gmail.com), " + \
        "[Zoltan Udvardy](zoltan.udvardy.ipbs@gmail.com)"
    contact = st.sidebar.markdown(contact_str)

    license_str = "**License:** ???"
    license = st.sidebar.markdown(license_str)

    project_str = "**Project Page:** [GitHub](https://github.com/michabirklbauer/internal-ions/)"
    project = st.sidebar.markdown(project_str)

    main_page()

if __name__ == "__main__":
    main()
