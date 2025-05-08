from pyteomics import mgf
from typing import Any, Dict, BinaryIO


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

    result_dict = {}

    print("Read spectra in total:")

    with mgf.read(filename, use_index=True) as reader:
        for s, spectrum in enumerate(reader):

            if (s + 1) % 1000 == 0:
                print(f"\t{s + 1}")

            scan_nr = s  # parse_scannr(spectrum["params"], -s, pattern)[1]
            spectrum_dict = {}
            spectrum_dict["spectrum"] = spectrum
            spectrum_dict["precursor"] = spectrum["params"]["pepmass"]
            spectrum_dict["charge"] = spectrum["params"]["charge"]
            spectrum_dict["rt"] = spectrum["params"].get("rtinseconds", 0.0)
            spectrum_dict["max_intensity"] = float(max(spectrum["intensity array"]))
            peaks = {}
            for i, mz in enumerate(spectrum["m/z array"]):
                peaks[mz] = spectrum["intensity array"][i]
            spectrum_dict["peaks"] = peaks
            result_dict[scan_nr] = spectrum_dict

    print(f"\nFinished reading {s + 1} spectra!")

    return {"name": name, "spectra": result_dict}


# TODO this can be optimized
def filter_spectra(mass_spectra: Dict[int, Any], filter_params: Dict[str, Any], name: str) -> Dict[str, Any]:
    """
    Returns a Dict including a list of spectra from pyteomics.mgf based on the given filter criteria:
    Dict["name": name,
         "filter_params": filter_params,
         "spectra": List[Dict]]
    """

    spectra = []

    print("Filtered spectra in total:")

    for s, key in enumerate(mass_spectra):
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
                break  # scans are in order, so no need to look anymore

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
