# STREAMLIT REDIRECT STD OUT / STD ERR
# 2021 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import sys
import streamlit as st
from io import StringIO
from threading import current_thread
from contextlib import contextmanager
from streamlit.runtime.scriptrunner_utils.script_run_context import SCRIPT_RUN_CONTEXT_ATTR_NAME


# redirect sys.stdout / sys.stderr
@contextmanager
def st_redirect(src, dst):
    placeholder = st.empty()
    output_func = getattr(placeholder, dst)

    with StringIO() as buffer:
        old_write = src.write

        def new_write(b):
            if getattr(current_thread(), SCRIPT_RUN_CONTEXT_ATTR_NAME, None):
                buffer.write(b)
                output_func(buffer.getvalue())
            else:
                old_write(b)

        try:
            src.write = new_write
            yield
        finally:
            src.write = old_write


# write sys.stdout to streamlit dst
@contextmanager
def st_stdout(dst):
    with st_redirect(sys.stdout, dst):
        yield


# write sys.stderr to streamlit dst
@contextmanager
def st_stderr(dst):
    with st_redirect(sys.stderr, dst):
        yield
