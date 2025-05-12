from pyteomics import mass

from psm_utils.io import FILETYPES
SUPPORTED_FILETYPES = list(FILETYPES)

REPO_NAME = "internal_ions"
DIV_COLOR = "rainbow"
FRAGANNOT_ION_NAMES = ["a", "b", "c", "cdot", "c-1", "c+1", "x", "y", "zdot", "z+1", "z+2", "z+3"]

ion_comp = mass.std_ion_comp.copy()
for k in list(ion_comp):
    if k.endswith('dot'):
        ion_comp[k.replace('-', '')] = ion_comp.pop(k)
ion_comp['t'] = {}

ion_cap_delta_mass = {
    name: mass.calculate_mass(composition=comp, absolute=False) for (name, comp) in ion_comp.items()
}

ion_direction = {
    "a": "n-term",
    "b": "n-term",
    "x": "c-term",
    "y": "c-term",
    "cdot": "n-term",
    "c": "n-term",
    "c-1": "n-term",
    "c+1": "n-term",
    "zdot": "c-term",
    "z+1": "c-term",
    "z+2": "c-term",
    "z+3": "c-term",
}


# ---------------------------------------------------------------------------- #
#                               For visualization                              #
# ---------------------------------------------------------------------------- #


colors = [
    "#9b2226",
    "#005f73",
    "#ee9b00",
    "#0a9396",
    "#94d2bd",
    "#ca6702",
    "#e9d8a6",
    "#bb3e03",
    "#001219",
    "#006BA6",
    "#35A7FF",
    "#EFA8B8",
    "#BFACC8",
    "#476A6F",
    "#7067CF",
    "#364156",
    "#98FB98",
    "#8A2BE2",
    "#35682D",
    "#252850",
    "#7E7B52",
]

colors = colors + (["#808080"] * 1000)
