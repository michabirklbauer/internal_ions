# internal_ions

## Setup

- Requirements: **This only runs on python 3.11 or lower!**
- Recommended: python 3.11.7
- Packages: `pip install -r requirements.txt`
- Packages (alternative, exact version numbers): `pip install -r python3117.txt`
- Packages (Ubuntu 22.04, python 3.10.13): `pip install -r Docker.config`
  - Downgraded rpy2 package because of missmatch with setuptools.
- You also need to have R installed!
  - Additionally the environmental variable `R_HOME` needs to be set and point to your R installation
  - Alternatively you can specify the path to your R installation in `R_HOME.config` and this will be done for you at runtime

## Development Notes

- Implement plot functions in `util/tab*` either in the `example_plot.py` script in there or preferably write your
  own plotting function in a seperate script and throw it in the folder (to avoid merge conflicts).
- Frontend functions should be implemented directly in `tab*.py`.
- Frontend functions that are used accross several tabs go into `utils`.
- Getting data from other tabs: Access via `st.session_state` (examples given) in the `tab*.py` scripts.

## Example Data

- Can be found in `data`
- For plotting we use `fragment_centric.csv`, `spectrum_centric.csv` (and optionally `result.json`).

## Running the App

- Install requirements!
- Run `streamlit run streamlit_app.py`

## Running the App via Docker

- Running this app via Docker is possible with `docker run -p 8501:8501 michabirklbauer/internalionsv2:latest`

## Server -> [89.58.32.151](http://89.58.32.151/)
