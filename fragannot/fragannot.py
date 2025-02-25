# Dependencies
import argparse
from argparse import ArgumentError
from logging import warning
import json
import numpy as np
import re
from typing import List  # type hinting
from pyteomics import mass
from itertools import tee
from random import uniform
from pyteomics import parser
import logging
import logging.config
import os
import ms_deisotope

# Local
from ._version import __version__
from . import constant
from .parser import Parser as Parser

## ------------------ UNSUPPORTED ---------------------

class Fragannot:
    def __init__(self):
        # Set up logging
        logging.config.dictConfig(self.__load_log_config())
        self.logger = logging.getLogger(__name__)

    def fragment_annotation(
        self,
        ident_file,
        spectra_file,
        tolerance,
        fragment_types,
        charges,
        losses,
        file_format,
        deisotope,
        write_file=True,
    ):
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

        P = Parser()

        psms = P.read(spectra_file, ident_file, file_format=file_format)
        i = 0

        psms_json = dict()

        for psm in psms:
            psms_json_key = psm.spectrum_id
            # skip if one psm of the same spectrum has already been done
            if psms_json_key in psms_json:
                psms_json[psms_json_key]["nr_idents_with_same_rank"] += 1
                continue

            if (i + 1) % 100 == 0:
                self.logger.info(f"{i + 1} spectra annotated")

            if charges == "auto":  # if charges to consider not specified: use precursor charge as max charge
                charges_used = range(1, abs(psm.get_precursor_charge()), 1)
            else:
                charges_used = charges

            theoretical_fragment_code = self.compute_theoretical_fragments(
                sequence_length=len(psm.peptidoform.sequence),
                fragment_types=fragment_types,
                charges=charges_used,
                neutral_losses=losses,
            )

            theoretical_fragment_dict = {
                f: self.theoretical_mass_to_charge(f, psm.peptidoform) for f in theoretical_fragment_code
            }

            if deisotope:  # deisotoping TODO check desotoping method to optimize
                mzs, intensities = self.deisotope_peak_list(
                    psm.spectrum["mz"].tolist(), psm.spectrum["intensity"].tolist()
                )
            else:
                mzs = psm.spectrum["mz"]
                intensities = psm.spectrum["intensity"].tolist()

            annotation_mz, annotation_code, annotation_count = self.match_fragments(
                mzs, theoretical_fragment_dict, tolerance=tolerance
            )

            psm.spectrum["intensity"] = intensities
            psm.spectrum["mz"] = mzs
            psm.spectrum["theoretical_mz"] = annotation_mz
            psm.spectrum["theoretical_code"] = annotation_code
            psm.spectrum["matches_count"] = annotation_count

            psms_json[psms_json_key] = \
                {
                    "sequence": psm.peptidoform.sequence,
                    "proforma": psm.peptidoform.proforma,
                    "annotation": psm.spectrum,
                    "spectrum_id": psm.spectrum_id,
                    "identification_score": psm.score,
                    "rank": psm.rank,
                    "nr_idents_with_same_rank": 1,
                    # "precursor_charge": int(psm.get_precursor_charge()),
                    "precursor_intensity": 666, # what??
                }
            i += 1

            # if i == 1000:
            #     break

        if write_file:
            with open(P.output_fname, "w", encoding="utf8") as f:
                json.dump(psms_json, f)

        return psms_json

    # Function

    def deisotope_peak_list(self, mzs, intensities):
        peaks = ms_deisotope.deconvolution.utils.prepare_peaklist(zip(mzs, intensities))
        deconvoluted_peaks, targeted = ms_deisotope.deconvolute_peaks(
            peaks, averagine=ms_deisotope.peptide, scorer=ms_deisotope.MSDeconVFitter(10.0), verbose=True
        )
        mzs = [p.mz for p in deconvoluted_peaks.peaks]
        intensities = [p.intensity for p in deconvoluted_peaks.peaks]
        return mzs, intensities

    def compute_theoretical_fragments(
        self,
        sequence_length: int,
        fragment_types: List[str],
        charges: List[int] = [-1],
        neutral_losses: List[str] = [],
        internal: bool = True,
    ) -> List[str]:

        ion_directions = constant.ion_direction

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

        internal_frags_with_nl = []

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

    def theoretical_mass_to_charge(self, fragment_code, peptidoform):

        start, end, ion_cap_start, ion_cap_end, charge, formula = self.parse_fragment_code(fragment_code)

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
        fragment_mass = (P + M + SI + EI + (H * charge) - L) / np.abs(charge)

        return fragment_mass

    def parse_fragment_code(self, fragment_code: str):

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

    def matching(self, mz1, mz2, tol):
        if abs(mz1 - mz2) <= tol:
            return True
        return False

    def match_fragments(self, exp_mz, theo_frag, tolerance):

        theo_frag = [[k, v] for k, v in sorted(theo_frag.items(), key=lambda item: item[1])]

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
                if self.matching(i, j[1], tolerance):
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
                        d[i], key=lambda t: t[2]
                    )  # add the only the annotation with the lowest mass error

                fragment_theoretical_code.append(closest[0])
                fragment_theoretical_mz.append(closest[1])
            else:
                fragment_theoretical_code.append(None)
                fragment_theoretical_mz.append(None)

            iter_2, last_match = tee(last_match)

        return (fragment_theoretical_mz, fragment_theoretical_code, fragment_theoretical_nmatch)

    def print_parameters(self, args):
        params = ""
        for arg, val in args.items():
            params = params + f"{arg} : {val} \n\t"
        self.logger.info(params)

    def __load_log_config(self):
        """Load log configurations from log_conf.json

        Returns:
        -------
        config : Configuration loaded fro, log_conf.json
        """
        config = {}
        with open(os.path.join(os.path.dirname(__file__), "log_conf.json"), "r", encoding="utf-8") as fd:
            config = json.load(fd)
        return config


# function called on start
def main(parser=argparse.ArgumentParser()):

    parser.add_argument("-v", "--version", action="store_true", help="Shows the app version.")
    parser.add_argument(
        "-i",
        "--identification",
        type=str,
        required=True,
        help="Path to spectra identification file in 'mzid' format",
    )
    parser.add_argument(
        "-s", "--spectra", type=str, required=True, help="Path to spectra file in 'mzml' or 'mgf' format"
    )

    parser.add_argument(
        "-t", "--tolerance", type=float, default=0.05, required=False, help="MS2 tolerance in Da"
    )

    parser.add_argument(
        "-f",
        "--fragment_types",
        type=str,
        default=["b", "y"],
        nargs="+",
        required=False,
        help="list of fragment ion types to be considered",
    )

    parser.add_argument(
        "-c",
        "--charges",
        type=str,
        default="auto",
        nargs="+",
        required=False,
        help="list of fragment charges to be considered",
    )
    parser.add_argument(
        "-F", "--format", type=str, default="infer", required=False, help="Identification file format"
    )

    parser.add_argument(
        "-l",
        "--losses",
        type=str,
        default=[],
        nargs="+",
        required=False,
        help="list molecular formula to be considered a possible neutral loss (e.g H20 for water loss)",
    )

    parser.add_argument(
        "-d",
        "--deisotope",
        type=bool,
        default=False,
        action=argparse.BooleanOptionalAction,
        required=False,
        help="Whether or not to perform deisotoping of the spectra",
    )

    args = parser.parse_args()

    if args.version:
        return __version__
    else:
        fragannot = Fragannot()

        welcome_msg = (
            "........................................\n\t"
            + "..lo:::'.;:c:;..........................\n\t"
            + "..odc:..;d:.:d;.........................\n\t"
            + "..ol'...;xo:ox,.........................\n\t"
            + "..oc....;d;.;d;.........................\n\t"
            + "..:;....':'.':'...'''............''.....\n\t"
            + ".................:ccc;..........,::'....\n\t"
            + ".....................'::'.....;:,.......\n\t"
            + "......................:l:;,;,;cl,.......\n\t"
            + "......................:l:,;:lc;,,''''...\n\t"
            + "......................:l' .;l;. 'clc;'..\n\t"
            + "......................:lc::clc::clc'....\n\t"
            + "......................:lclccllc,;:'.....\n\t"
            + ".......'..........'..':cccccc:;;;,'.....\n\t"
            + "......;l;......';:ccc::::::::;;,........\n\t"
            + "......;l;....';cccclllcccccccc;.........\n\t"
            + "......;c;..,ccccllcclllllllll:..........\n\t"
            + "........';:ccllllllllllllllll:..........\n\t"
            + ".........',;cllcc:;'';cl:,,;c:..........\n\t"
            + "...........,cllc:;....:l,...;:..........\n\t"
        )
        fragannot.logger.info("-= Starting fragannot =-")
        fragannot.logger.info(welcome_msg)

        fragannot.print_parameters(vars(args))
        fragannot.fragment_annotation(
            ident_file=args.identification,
            spectra_file=args.spectra,
            tolerance=args.tolerance,
            fragment_types=args.fragment_types,
            charges=args.charges,
            losses=args.losses,
            file_format=args.format,
            deisotope=args.deisotope,
        )


def annotate_fragments(
    ident_file, spectra_file, tolerance, fragment_types, charges, losses, file_format, write_file=True
):
    fragannot = Fragannot()
    return fragannot.fragment_annotation(
        ident_file=ident_file,
        spectra_file=spectra_file,
        tolerance=tolerance,
        fragment_types=fragment_types,
        charges=charges,
        losses=losses,
        file_format=file_format,
        write_file=write_file,
    )
