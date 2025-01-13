#!/usr/bin/env python3

import re
from pyteomics import mgf

from typing import Any
from typing import List
from typing import Dict
from typing import Tuple
from typing import BinaryIO

# parse scan number from pyteomics mgf params
def parse_scannr(obj: dict | str | int,  i: int, pattern: str = "\\.\\d+\\.") -> Tuple[int, int]:
    """Parses the scan number from the params dictionary of the pyteomics mgf
    spectrum.

    Parameters
    ----------
    object : dict, str, or int
        The "params" dictionary of the pyteomics mgf spectrum, or a `psm_utils.spectrum_id`.

    i : int
        The scan number to be returned in case of failure.

    pattern : str
        Regex pattern to use for parsing the scan number from the title if it
        can't be infered otherwise.

    Returns
    -------
    (exit_code, scan_nr) : Tuple
        A tuple with the exit code (0 if successful, 1 if parsing failed) at the
        first position [0] and the scan number at the second position [1].
    """
    
    # if input is already a number, assume it's a scan number
    if type(obj) == int:
        return obj
    
    # if input is a string, assume it's spectrum title or spectrum id
    if type(obj) == str:
        # if there is a scan token in the title, try parse scan_nr
        if "scan" in obj:
            try:
                return (0, int(obj.split("scan=")[1].strip("\"")))
            except:
                pass
    
        # else try to parse by pattern
        try:
            scan_nr = re.findall(pattern, obj)[0]
            scan_nr = re.sub(r"[^0-9]", "", scan_nr)
            if len(scan_nr) > 0:
                return (0, int(scan_nr))
        except:
            pass
    
        # else try parse whole title
        try:
            return (0, int(float(obj)))
        except:
            pass

    # if input is dictionary, assume it's mgf params
    if type(obj) == dict:
        params = obj
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

            # else try to parse by pattern
            try:
                scan_nr = re.findall(pattern, params["title"])[0]
                scan_nr = re.sub(r"[^0-9]", "", scan_nr)
                if len(scan_nr) > 0:
                    return (0, int(scan_nr))
            except:
                pass

            # else try parse whole title
            try:
                return (0, int(params["title"]))
            except:
                pass

    # return unsuccessful parse
    return (1, i)

# reading spectra
def read_spectra(filename: str | BinaryIO, name: str, pattern: str = "\\.\\d+\\.") -> Dict[int, Dict]:
    """
    Returns a dictionary that maps scan numbers to spectra:
    Dict["name": name,
         "spectra": Dict[int -> Dict["spectrum"         -> pyteomics mgf spectrum
                                     "precursor"        -> float
                                     "charge"           -> int
                                     "rt"               -> float
                                     "max_intensity"    -> float
                                     "peaks"            -> Dict[m/z -> intensity]]
    """

    result_dict = dict()

    print("Read spectra in total:")

    with mgf.read(filename, use_index = True) as reader:
        for s, spectrum in enumerate(reader):

            if (s + 1) % 1000 == 0:
                print(f"\t{s + 1}")

            scan_nr = parse_scannr(spectrum["params"], -s, pattern)[1]
            spectrum_dict = dict()
            spectrum_dict["spectrum"] = spectrum
            spectrum_dict["precursor"] = spectrum["params"]["pepmass"]
            spectrum_dict["charge"] = spectrum["params"]["charge"]
            spectrum_dict["rt"] = 0.0 if "rtinseconds" not in spectrum["params"] else spectrum["params"]["rtinseconds"]
            spectrum_dict["max_intensity"] = float(max(spectrum["intensity array"]))
            peaks = dict()
            for i, mz in enumerate(spectrum["m/z array"]):
                peaks[mz] = spectrum["intensity array"][i]
            spectrum_dict["peaks"] = peaks
            result_dict[scan_nr] = spectrum_dict
        reader.close()

    print(f"\nFinished reading {s + 1} spectra!")

    return {"name": name, "spectra": result_dict}

# TODO this can be optimized
def filter_spectra(mass_spectra: Dict[int, Any], filter_params: Dict[str, Any], name: str, pattern: str = "\\.\\d+\\.") -> Dict[str, Any]:
    """
    Returns a Dict including a list of spectra from pyteomics.mgf based on the given filter criteria:
    Dict["name": name,
         "filter_params": filter_params,
         "spectra": List[Dict]]
    """

    spectra = []

    print("Filtered spectra in total:")

    for s, key in enumerate(mass_spectra.keys()):
        spectrum = mass_spectra[key]["spectrum"]
        
        if (s + 1) % 1000 == 0:
            print(f"\t{s + 1}")

        scan_nr = key

        # check spectrum > first scan
        if "first_scan" in filter_params:
            if scan_nr < int(filter_params["first_scan"]):
                continue

        if "last_scan" in filter_params:
            if scan_nr > int(filter_params["last_scan"]):
                continue

        if "min_mz" in filter_params:
            if float(spectrum["params"]["pepmass"][0]) < float(filter_params["min_mz"]):
                continue

        if "max_mz" in filter_params:
            if float(spectrum["params"]["pepmass"][0]) > float(filter_params["max_mz"]):
                continue

        if "min_rt" in filter_params and "rtinseconds" in spectrum["params"]:
            if float(spectrum["params"]["rtinseconds"]) < float(filter_params["min_rt"]):
                continue

        if "max_rt" in filter_params and "rtinseconds" in spectrum["params"]:
            if float(spectrum["params"]["rtinseconds"]) > float(filter_params["max_rt"]):
                continue

        if "max_charge" in filter_params:
            # how to handle mutiple charges?
            for charge in spectrum["params"]["charge"]:
                if int(charge) > int(filter_params["max_charge"]):
                    continue

        if "max_isotope" in filter_params:
            pass
            # todo implement

        if "scans" in filter_params:
            if scan_nr not in filter_params["scans"]:
                continue

        spectra.append(spectrum)

    print(f"\nFinished filtering {s + 1} spectra in total!")

    return {"name": name, "filter_params": filter_params, "spectra": spectra}
