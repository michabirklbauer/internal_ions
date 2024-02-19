#!/usr/bin/env python3

from pyteomics import mgf

from typing import Dict
from typing import Tuple
from typing import BinaryIO

# parsing the scan nr
def parse_scannr(params: Dict, i: int) -> Tuple[int, int]:
    """
    Returns (0, scan nr) if scan number was successfully parsed.
    Otherwise returns (1, i) if scan number could not be parsed.
    """

    # prefer scans attr over title attr
    if "scans" in params:
        try:
            return (0, int(params["scans"]))
        except:
            pass

    # try parse title
    if "title" in params:

        # if there is a scan token in the title, try parse scan_nr
        if "scan" in params["title"]:
            try:
                return (0, int(params["title"].split("scan=")[1].strip("\"")))
            except:
                pass

        # else try parse whole title
        try:
            return (0, int(params["title"]))
        except:
            pass

    # return insuccessful parse
    return (1, i)

# reading spectra
def read_spectra(filename: str | BinaryIO, name: str) -> Dict[int, Dict]:
    """
    Returns a dictionary that maps scan numbers to spectra:
    Dict["name": name,
         "spectra": Dict[int -> Dict["precursor"        -> float
                                     "charge"           -> int
                                     "rt"               -> float
                                     "max_intensity"    -> float
                                     "peaks"            -> Dict[m/z -> intensity]]
    """

    result_dict = dict()

    with mgf.read(filename, use_index = True) as reader:
        for s, spectrum in enumerate(reader):
            scan_nr = parse_scannr(spectrum["params"], -s)[1]
            spectrum_dict = dict()
            spectrum_dict["precursor"] = spectrum["params"]["pepmass"]
            spectrum_dict["charge"] = spectrum["params"]["charge"]
            spectrum_dict["rt"] = spectrum["params"]["rtinseconds"]
            spectrum_dict["max_intensity"] = float(max(spectrum["intensity array"]))
            peaks = dict()
            for i, mz in enumerate(spectrum["m/z array"]):
                peaks[mz] = spectrum["intensity array"][i]
            spectrum_dict["peaks"] = peaks
            result_dict[scan_nr] = spectrum_dict
        reader.close()

    return {"name": name, "spectra": result_dict}
