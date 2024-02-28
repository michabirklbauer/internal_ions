#!/usr/bin/env python3

from pyteomics import mgf

from typing import Any
from typing import List
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

    print("Read spectra in total:")

    with mgf.read(filename, use_index = True) as reader:
        for s, spectrum in enumerate(reader):

            if (s + 1) % 1000 == 0:
                print(f"\t{s + 1}")

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

    print(f"\nFinished reading {s + 1} spectra!")

    return {"name": name, "spectra": result_dict}

def filter_spectra(filename: str | BinaryIO, filter_params: Dict[str, Any], name: str) -> Dict[str, Any]:
    """
    Returns a Dict including a list of spectra from pyteomics.mgf based on the given filter criteria:
    Dict["name": name,
         "filter_params": filter_params,
         "spectra": List[Dict]]
    """

    spectra = []

    print("Filtered spectra in total:")

    with mgf.read(filename, use_index = True) as reader:
        for s, spectrum in enumerate(reader):

            if (s + 1) % 1000 == 0:
                print(f"\t{s + 1}")

            scan_nr = parse_scannr(spectrum["params"], -s)[1]

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

            if "min_rt" in filter_params:
                if float(spectrum["params"]["rtinseconds"]) < float(filter_params["min_rt"]):
                    continue

            if "max_rt" in filter_params:
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

            if "scans_from_protein" in filter_params:
                if scan_nr not in filter_params["scans_from_protein"]:
                    continue

            if "scans_from_peptide" in filter_params:
                if scan_nr not in filter_params["scans_from_peptide"]:
                    continue

            spectra.append(spectrum)
        reader.close()

    print(f"\nFinished filtering {s + 1} spectra in total!")

    return {"name": name, "filter_params": filter_params, "spectra": spectra}

def read_filtered_spectra(filename: str | BinaryIO, name: str) -> Dict[str, Any]:

    filter_params = {"first_scan": None,
                     "last_scan": None,
                     "min_mz": None,
                     "max_mz": None,
                     "min_rt": None,
                     "max_rt": None,
                     "max_charge": None,
                     "max_isotope": None,
                     "selected_protein": None,
                     "scans_from_protein": [],
                     "selected_peptide": None,
                     "scans_from_peptide": []}

    spectra = []

    print("Read spectra in total:")

    with mgf.read(filename, use_index = True) as reader:
        for s, spectrum in enumerate(reader):

            if (s + 1) % 1000 == 0:
                print(f"\t{s + 1}")

            spectra.append(spectrum)
        reader.close()

    print(f"\nFinished reading {s + 1} spectra!")

    return {"name": name, "filter_params": filter_params, "spectra": spectra}
