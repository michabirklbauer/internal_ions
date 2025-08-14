import streamlit as st
from pyteomics import mgf
import logging
from collections import Counter
from io import StringIO
import os
from typing import Any, Dict, BinaryIO, TextIO


class SpectrumFile:
    def __init__(self, uploaded_file: BinaryIO):
        self.logger = logging.getLogger(__name__)
        self.file_format = None
        self.name = uploaded_file.name
        self._load(uploaded_file)

    def get_mgf_index_by_scans(self, file: str | TextIO):
        """
        Determine if an MGF file should be indexed by scans or by title.

        Parameters
        ----------
        file : MGF file

        Returns
        -------
        out : bool
            True if scans should be used and False if title should be used

        """
        scans = Counter()
        titles = Counter()
        with mgf.MGF(file) as f:
            for i, spec in enumerate(f):
                if 'scans' in spec['params']:
                    scans[spec['params']['scans']] += 1
                if 'title' in spec['params']:
                    titles[spec['params']['title']] += 1
        use_scans = len(scans) >= len(titles)
        # TODO: indexing by both scans and titles is probably necessary
        if max(len(scans), len(titles)) < i + 1:
            self.logger.warning("Not all spectra can be indexed in %s: %d scans, %d titles, %d spectra", file, len(scans), len(titles), i + 1)
        return use_scans

    def _load(self, uploaded_file):
        uploaded_file.seek(0)
        extension = os.path.splitext(uploaded_file.name)[1]

        if extension.lower() == ".mgf":
            self.logger.info(f"Inferred MGF format from {uploaded_file.name}")
            index_by_scans = self.get_mgf_index_by_scans(StringIO(uploaded_file.getvalue().decode('utf-8')))
            self.spectra_source = mgf.IndexedMGF(uploaded_file, index_by_scans=index_by_scans)
            self.file_format = 'mgf'
            self.index = self.spectra_source.index

        else:
            self.logger.error("Cannot infer format from %s, only MGF format is supported", uploaded_file.name)
            raise Exception("Unsupported spectra file format")

        self.logger.info(f"Loaded a SpectrumFile with {len(self.spectra_source)} spectra.")

    def get_by_id(self, i):
        return self.spectra_source[i]

    def __iter__(self):
        return iter(self.spectra_source)

    def __next__(self):
        return next(self.spectra_source)


@st.cache_data
def read_spectrum_file(uploaded_file: BinaryIO):
    return SpectrumFile(uploaded_file)


@st.cache_data(hash_funcs={SpectrumFile: lambda x: x.name})
def read_spectra(file: SpectrumFile) -> Dict[int, Dict]:
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

    # print("Read spectra in total:")

    for s, spectrum in enumerate(file):

        # if (s + 1) % 1000 == 0:
        #     print(f"\t{s + 1}")

        scan_nr = s
        spectrum_dict = {}
        spectrum_dict["spectrum"] = spectrum
        spectrum_dict["precursor"] = spectrum["params"]["pepmass"]
        spectrum_dict["charge"] = spectrum["params"].get("charge", 0)
        spectrum_dict["rt"] = spectrum["params"].get("rtinseconds", 0.0)
        spectrum_dict["max_intensity"] = float(max(spectrum["intensity array"]))
        peaks = {}
        for i, mz in enumerate(spectrum["m/z array"]):
            peaks[mz] = spectrum["intensity array"][i]
        spectrum_dict["peaks"] = peaks
        result_dict[scan_nr] = spectrum_dict

    # print(f"\nFinished reading {s + 1} spectra!")

    return {"name": file.name, "spectra": result_dict}


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
