# internal_ions

## Setup

- Requirements: **python 3.11+!**
- Packages: `pip install -r requirements.txt`
- Packages (alternative, exact version numbers): `pip install -r env.txt`

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

- Running this app via Docker is possible with:
  - `docker build . -f Dockerfile -t internalions`
  - `docker run -p 8501:8501 internalions`

## Server

- [Online](https://computproteomics.bmb.sdu.dk/app_direct/internal_ions/)
