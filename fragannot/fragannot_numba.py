# Dependencies
import json
import numpy as np
import re
from pyteomics import mass
from itertools import tee
from random import uniform
from pyteomics import parser
import os
import ms_deisotope

# type hinting
from typing import List
from typing import Dict
from typing import Any

from . import constant
from .parser import Parser as Parser

from numba import jit
from numba.typed import List as numbaList

import multiprocessing
from joblib import Parallel, delayed
from tqdm import tqdm
from stqdm import stqdm

class FragannotNumba:
    # set CPU cores here
    def __init__(self, reserved_cores: int = 2):
        self.nr_used_cores = multiprocessing.cpu_count() - reserved_cores

    def fragment_annotation(
        self,
        ident_file: str,
        spectra_file: str,
        tolerance: float,
        fragment_types: List[str],
        charges: List[str],
        losses: List[str],
        file_format: str,
        deisotope: bool,
        write_file: bool = True) -> List[Dict[str, Any]]:

        return fragment_annotation(ident_file, spectra_file, tolerance,
                                   fragment_types, charges, losses, file_format,
                                   deisotope, write_file, self.nr_used_cores)

# set micro batching and batch params here
def fragment_annotation(
    ident_file: str,
    spectra_file: str,
    tolerance: float,
    fragment_types: List[str],
    charges: List[str] | str,
    losses: List[str],
    file_format: str,
    deisotope: bool,
    write_file: bool = True,
    nr_used_cores: int = 1,
    micro_batch: bool = True,
    batch_size: int = 100) -> List[Dict[str, Any]]:
    """
    Annotate theoretical and observed fragment ions in a spectra file.

    Parameters:
    ----------
    ident_file : str
        Filename of an identification file
    spectra_file : str
        Filename of a spectra file
    tolerance : float
        Tolerance value in ppm for fragment matching
    fragment_types : list
        List of fragment types (fragment type must be defined in constant.py ion_cap_formula)
    charges : list
        List of charges (e.g ["+1", "-2"])
    losses : list
        List of neutral losses molecular formula (e.g ["H2O"])
    file_format : str
        String indicating the file format of the input files

    Returns:
    -------
    None
    """

    print("Fragannot running using " + str(nr_used_cores) + " logical cores.\n")

    P = Parser(is_streamlit = True)

    psms = P.read(spectra_file, ident_file, file_format = file_format)

    print("\nAnnotating spectra in parallel...\n")

    if micro_batch:
        i = 0
        still_spectra_available = True
        psms_json = []
        print("Annotated spectra in total:")
        while still_spectra_available:
            print(f" {i}\t")
            if i + batch_size < len(psms):
                current_batch = psms[i:i+batch_size]
            else:
                current_batch = psms[i:len(psms)]
                still_spectra_available = False
            p_result = Parallel(n_jobs = nr_used_cores)(delayed(calculate_ions_for_psms)(psm, tolerance, fragment_types, charges, losses, deisotope) for psm in current_batch)
            psms_json += list(p_result)
            i += batch_size
    else:
        p_psms = tqdm(psms) # tqdm is good for cli but bad for streamlit progress
        p_result = Parallel(n_jobs = nr_used_cores)(delayed(calculate_ions_for_psms)(psm, tolerance, fragment_types, charges, losses, deisotope) for psm in p_psms)
        psms_json = list(p_result)

    if write_file:
        with open(P.output_fname, "w", encoding = "utf8") as f:
            json.dump(psms_json, f)

    print("\nFinished spectrum annotation.")

    return psms_json

def calculate_ions_for_psms(psm,
                            tolerance: float,
                            fragment_types: List[str],
                            charges: List[str] | str,
                            losses: List[str],
                            deisotope: bool) -> Dict[str, Any]:

    if charges == "auto":  # if charges to consider not specified: use precursor charge as max charge
        charges_used = range(1, abs(psm.get_precursor_charge()), 1)
    else:
        charges_used = charges

    theoretical_fragment_code = compute_theoretical_fragments(
        sequence_length = len(psm.peptidoform.sequence),
        fragment_types = numbaList(fragment_types),
        charges = numbaList([int(c) for c in charges_used]),
        neutral_losses = numbaList(losses),
        internal = True
    )

    theoretical_fragment_dict = {
        f: theoretical_mass_to_charge(f, psm.peptidoform) for f in theoretical_fragment_code
    }

    if deisotope:  # deisotoping TODO check desotoping method to optimize
        mzs, intensities = deisotope_peak_list(
            psm.spectrum["mz"].tolist(), psm.spectrum["intensity"].tolist()
        )
    else:
        mzs = psm.spectrum["mz"]
        intensities = psm.spectrum["intensity"].tolist()

    annotation_mz, annotation_code, annotation_count = match_fragments(
        mzs, theoretical_fragment_dict, tolerance = tolerance
    )

    psm.spectrum["intensity"] = intensities
    psm.spectrum["mz"] = mzs
    psm.spectrum["theoretical_mz"] = annotation_mz
    psm.spectrum["theoretical_code"] = annotation_code
    psm.spectrum["matches_count"] = annotation_count

    return {"sequence": psm.peptidoform.sequence,
            "proforma": psm.peptidoform.proforma,
            "annotation": psm.spectrum,
            "spectrum_id": psm.spectrum_id,
            "identification_score": psm.score,
            "rank": psm.rank,
            # "precursor_charge": int(psm.get_precursor_charge()),
            "precursor_intensity": 666}

def deisotope_peak_list(mzs: List[float], intensities: List[float]) -> List[List[float]]:
    peaks = ms_deisotope.deconvolution.utils.prepare_peaklist(zip(mzs, intensities))
    deconvoluted_peaks, targeted = ms_deisotope.deconvolute_peaks(
        peaks, averagine = ms_deisotope.peptide, scorer = ms_deisotope.MSDeconVFitter(10.0), verbose = True
    )
    mzs = [p.mz for p in deconvoluted_peaks.peaks]
    intensities = [p.intensity for p in deconvoluted_peaks.peaks]

    return mzs, intensities

@jit(nopython = True, cache = True)
def compute_theoretical_fragments(
    sequence_length: int,
    fragment_types: List[str],
    charges: List[int] = [-1],
    neutral_losses: List[str] = [],
    internal: bool = True) -> List[str]:

    #ion_directions = constant.ion_direction
    ion_directions = {
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

    n_term_ions = [ion_type for ion_type in fragment_types if ion_directions[ion_type] == "n-term"]
    c_term_ions = [ion_type for ion_type in fragment_types if ion_directions[ion_type] == "c-term"]

    n_term = ["t:" + ion_type for ion_type in n_term_ions]
    c_term = [ion_type + ":t" for ion_type in c_term_ions]

    # terminal fragments
    n_term_frags = [
        n_term_frag + "@1:" + str(i + 1) for n_term_frag in n_term for i in range(sequence_length - 1)
    ]
    c_term_frags = [
        c_term_frag + "@" + str(i) + ":" + str(sequence_length)
        for c_term_frag in c_term
        for i in range(2, sequence_length + 1)
    ]

    charges_str = [
        "(" + str(int(charge)) + ")" if int(charge) < 0 else "(+" + str(int(charge)) + ")"
        for charge in charges
    ]

    n_term_frags_with_charges = [
        n_term_frag + charge for n_term_frag in n_term_frags for charge in charges_str
    ]

    c_term_frags_with_charges = [
        c_term_frag + charge for c_term_frag in c_term_frags for charge in charges_str
    ]

    neutral_losses_str = ["[" + nl + "]" for nl in neutral_losses]
    neutral_losses_str.append("")
    n_term_frags_with_nl = [
        n_term_frag + nl for n_term_frag in n_term_frags_with_charges for nl in neutral_losses_str
    ]
    c_term_frags_with_nl = [
        c_term_frag + nl for c_term_frag in c_term_frags_with_charges for nl in neutral_losses_str
    ]

    internal_frags_with_nl = ["" for x in range(0)]

    if internal:
        # internal fragments
        internal = [
            n_term_ion + ":" + c_term_ion for n_term_ion in n_term_ions for c_term_ion in c_term_ions
        ]
        internal_pos = [
            str(i) + ":" + str(j)
            for i in range(2, sequence_length)
            for j in range(2, sequence_length)
            if i <= j
        ]
        internal_frags = [
            internal_ions + "@" + internal_positions
            for internal_ions in internal
            for internal_positions in internal_pos
        ]

        internal_frags_with_charges = [
            internal_frag + charge for internal_frag in internal_frags for charge in charges_str
        ]

        internal_frags_with_nl = [
            internal_frag + nl
            for internal_frag in internal_frags_with_charges
            for nl in neutral_losses_str
        ]

    return n_term_frags_with_nl + c_term_frags_with_nl + internal_frags_with_nl

def theoretical_mass_to_charge(fragment_code: str, peptidoform) -> float:

    start, end, ion_cap_start, ion_cap_end, charge, formula = parse_fragment_code(fragment_code)

    # peptide and modification mass
    sequence = []
    mods = []
    for aa, mod in peptidoform.parsed_sequence[start - 1 : end]:
        sequence.append(aa)
        if not mod is None:
            mods.extend([m.mass for m in mod])

    # mass AA sequence
    ps = parser.parse("".join(sequence), show_unmodified_termini=True)
    P = mass.calculate_mass(parsed_sequence=ps)
    # mass modifications
    M = sum(mods)
    # mass start ion cap
    SI = constant.ion_cap_delta_mass[ion_cap_start]
    # mass end ion cap
    EI = constant.ion_cap_delta_mass[ion_cap_end]
    # hydrogen mass
    H = 1.00784
    # loss mass
    L = mass.calculate_mass(formula, absolute=True)

    # Calculate fragment mass
    fragment_mass = (P + M + SI + EI + (H * charge) - L) / abs(charge)

    return fragment_mass

def parse_fragment_code(fragment_code: str):

    # test if fragment code format is valid*
    fragment_code_pattern = re.compile(r".+(:).+(@)[0-9]+(:)[0-9]+(\()(\+|\-)[0-9](\))(\[(.*?)\])?")
    if bool(fragment_code_pattern.match(fragment_code)) == False:
        raise RuntimeError("Incorrect fragment code format: {0}".format(fragment_code))

    ## Parse fragment code

    start, end = [
        int(i) for i in re.search(r"(?<=\@)(.*?)(?=\()", fragment_code).group(1).split(":")
    ]  # Get start and end amino acid indexes
    ion_cap_start, ion_cap_end = [
        str(i) for i in re.search(r"^(.*?)(?=\@)", fragment_code).group(1).split(":")
    ]  # Get start and end ion caps name
    charge = int(re.search(r"(?<=\()(.*?)(?=\))", fragment_code).group(1))  # get charge state
    formula = re.search(r"(?<=\[)(.*?)(?=\])", fragment_code)
    if formula == None:
        formula = ""
    else:
        formula = str(re.search(r"(?<=\[)(.*?)(?=\])", fragment_code).group(1))

    return start, end, ion_cap_start, ion_cap_end, charge, formula

def match_fragments(exp_mz, theo_frag, tolerance: float):

    theo_frag = [[k, v] for k, v in sorted(theo_frag.items(), key = lambda item: item[1])]

    re_term = re.compile(r"^t:|:t")

    iter_2, last_match = tee(iter(theo_frag))

    d = {}

    fragment_theoretical_code = []
    fragment_theoretical_mz = []
    fragment_theoretical_nmatch = []

    for i in exp_mz:
        d.setdefault(i, [])
        found = False
        while True:
            j = next(iter_2, (None, None))

            # print(j)
            if j[1] is None:
                break
            if abs(i - j[1]) <= tolerance:
                k = [j[0], j[1], abs(i - j[1])]

                d[i].append(k)
                if not found:
                    iter_2, last_match = tee(iter_2)
                    found = True
            else:
                if found:
                    break

        fragment_theoretical_nmatch.append(len(d[i]))
        if len(d[i]) > 0:

            closest = None
            for frag in d[i]:

                if re_term.search(frag[0]):  # Prioritize annotation of terminal ions
                    closest = frag
                    break

            if closest is None:
                closest = min(
                    d[i], key = lambda t: t[2]
                )  # add the only the annotation with the lowest mass error

            fragment_theoretical_code.append(closest[0])
            fragment_theoretical_mz.append(closest[1])
        else:
            fragment_theoretical_code.append(None)
            fragment_theoretical_mz.append(None)

        iter_2, last_match = tee(last_match)

    return (fragment_theoretical_mz, fragment_theoretical_code, fragment_theoretical_nmatch)
