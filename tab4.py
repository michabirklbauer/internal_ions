import streamlit as st

from util.constants import DIV_COLOR

def main(argv = None) -> None:

    readme = ""

    with open("README.md", "r") as f:
        readme = f.read()
        f.close()

    documentation = st.markdown(readme)

    divider = st.subheader("", divider = DIV_COLOR)

    col_l, col_m, col_r = st.columns(3)

    with col_m:
        ext_docs = st.link_button("Read full documentation!",
                                  url = "https://internal-ions.vercel.app/",
                                  type = "primary",
                                  use_container_width = True)
