import os
import re
import json
import pandas as pd
from functools import lru_cache
from typing import List, Dict, Union, BinaryIO

__version = "1.0.0"
__date = "2023-02-07"


class JSONConverter:
    """
    Convert JSON files to different (file) formats.
    Example usage:
        json_file = "data.json"
        converter = JSONConverter()
        converter.load(json_file)
        dataframes = converter.to_dataframes()
        converter.to_csv(prefix = "data")
        dataframes[0].head()
        dataframes[1].head()
    """
    fragment_code_pattern = re.compile(r"(?P<ion_cap_start>[^:]+):(?P<ion_cap_end>[^:]+)@(?P<start>[0-9]+):(?P<end>[0-9]+)\((?P<charge>(?:\+|\-)?[0-9]+)\)(\[(?P<formula>.*)\])?")

    def __init__(self, input: Union[str, BinaryIO] = None) -> None:
        """
        Constructor.
        If 'input' is provided, self.load(input) will be called.
        """
        self.data = None
        if input is not None:
            self.load(input)

    def load(self, input: Union[str, BinaryIO]) -> Dict:
        """
        Load JSON from file or filestream (auto-detect).
        """

        if os.path.isfile(input):
            self.data = self.load_from_file(input)
        else:
            self.data = self.load_from_stream(input)

        return self.data

    def load_from_file(self, json_file: str) -> Dict:
        """
        Load JSON from file.
        """

        with open(json_file, "r", encoding="utf-8") as f:
            self.data = json.load(f)
            f.close()

        return self.data

    def load_from_stream(self, json_filestream: BinaryIO) -> Dict:
        """
        Load JSON from filestream.
        """

        self.data = json.load(json_filestream)

        return self.data

    def to_csv(self, prefix: str, data: Dict = None) -> List[pd.DataFrame]:
        """
        Export JSON to CSV format.
        """

        if data is None:
            data = self.data

        dataframes = self.to_dataframes(data)
        dataframes[0].to_csv(prefix + "_fragment.csv", index=False)
        dataframes[1].to_csv(prefix + "_spectrum.csv", index=False)

        return dataframes

    def to_excel(self, prefix: str, data: Dict = None) -> List[pd.DataFrame]:
        """
        Export JSON to EXCEL format.
        """

        if data is None:
            data = self.data

        dataframes = self.to_dataframes(data)
        dataframes[0].to_excel(prefix + "_fragment.xlsx")
        dataframes[1].to_excel(prefix + "_spectrum.xlsx")

        return dataframes

    def to_dataframes(self, data: Dict = None) -> List[pd.DataFrame]:
        """
        Re-format JSON to dataframes.
        """

        if data is None:
            data = self.data

        # fragment-centric dataframe structure
        fragment = {"frag_code": [],
                    "frag_seq": [],
                    "frag_type1": [],
                    "frag_type2": [],
                    "position_frag_type1": [],
                    "position_frag_type2": [],
                    "frag_length": [],
                    "frag_charge": [],
                    "frag_intensity": [],
                    "frag_mz": [],
                    "perc_of_total_intensity": [],
                    "prop_intensity_to_base_peak": [],
                    "modification": [],
                    "spectrum_id": [],
                    "ambiguity": [],
                    "nr_idents_with_same_rank": []}

        # spectrum-centric dataframe structure
        spectrum = {"perc_internal": [],
                    "perc_terminal": [],
                    "perc_other": [],
                    "total_int_internal": [],
                    "total_int_terminal": [],
                    "total_int_other": [],
                    "top1_internal_ion_code": [],
                    "top1_internal_seq": [],
                    "top2_internal_ion_code": [],
                    "top2_internal_seq": [],
                    "top3_internal_ion_code": [],
                    "top3_internal_seq": [],
                    "intensity_explained_aa": [],
                    "intensity_explained_precursor": [],
                    "spectrum_id": [],
                    "score": [],
                    "peptide_seq": [],
                    "proforma": [],
                    "peptide_length": [],
                    "nr_idents_with_same_rank": []}

        for key in data.keys():
            entry = data[key]
            pep_seq = entry["sequence"]

            # Fragment-centric
            for i, code in enumerate(entry["annotation"]["theoretical_code"]):
                # skip non-annonated fragments
                if code is None:
                    fragment["frag_code"].append("n:n@0:0(+0)")
                    fragment["frag_seq"].append("")
                    fragment["frag_type1"].append("n")
                    fragment["frag_type2"].append("n")
                    fragment["position_frag_type1"].append(0)
                    fragment["position_frag_type2"].append(0)
                    fragment["frag_length"].append(0)
                    fragment["frag_charge"].append(0)
                    fragment["frag_intensity"].append(entry["annotation"]["intensity"][i])
                    fragment["frag_mz"].append(entry["annotation"]["mz"][i])
                    total_intensity = sum(entry["annotation"]["intensity"])
                    fragment["perc_of_total_intensity"].append(entry["annotation"]["intensity"][i] / total_intensity)
                    intensity_base_peak = max(entry["annotation"]["intensity"])
                    fragment["prop_intensity_to_base_peak"].append(entry["annotation"]["intensity"][i] / intensity_base_peak)
                    fragment["modification"].append("")
                    fragment["spectrum_id"].append(self._get_spectrum_id(entry))
                    fragment["ambiguity"].append(self._get_ambiguity(entry["annotation"], i))
                    fragment["nr_idents_with_same_rank"].append(entry["nr_idents_with_same_rank"])
                else:
                    start, end, ion_cap_start, ion_cap_end, charge, formula = self._parse_fragment_code(code)
                    frag_seq = pep_seq[start - 1:end]
                    fragment["frag_code"].append(code)
                    fragment["frag_seq"].append(frag_seq)
                    fragment["frag_type1"].append(ion_cap_start)
                    fragment["frag_type2"].append(ion_cap_end)
                    fragment["position_frag_type1"].append(start)
                    fragment["position_frag_type2"].append(end)
                    fragment["frag_length"].append(end - start + 1)
                    fragment["frag_charge"].append(charge)
                    fragment["frag_intensity"].append(entry["annotation"]["intensity"][i])
                    fragment["frag_mz"].append(entry["annotation"]["mz"][i])
                    total_intensity = sum(entry["annotation"]["intensity"])
                    fragment["perc_of_total_intensity"].append(entry["annotation"]["intensity"][i] / total_intensity)
                    intensity_base_peak = max(entry["annotation"]["intensity"])
                    fragment["prop_intensity_to_base_peak"].append(entry["annotation"]["intensity"][i] / intensity_base_peak)
                    fragment["modification"].append(self._parse_modfication(entry["proforma"], start, end))
                    fragment["spectrum_id"].append(self._get_spectrum_id(entry))
                    fragment["ambiguity"].append(self._get_ambiguity(entry["annotation"], i))
                    fragment["nr_idents_with_same_rank"].append(entry["nr_idents_with_same_rank"])

            # Spectrum-centric
            percentages_and_total_intensities = self._calculate_internal_terminal_non_annotated_ions(entry["annotation"]["theoretical_code"], entry["annotation"]["intensity"])
            spectrum["perc_internal"].append(percentages_and_total_intensities["internal"])
            spectrum["perc_terminal"].append(percentages_and_total_intensities["terminal"])
            spectrum["perc_other"].append(percentages_and_total_intensities["non_annotated"])
            spectrum["total_int_internal"].append(percentages_and_total_intensities["total_int_internal"])
            spectrum["total_int_terminal"].append(percentages_and_total_intensities["total_int_terminal"])
            spectrum["total_int_other"].append(percentages_and_total_intensities["total_int_non"])
            top_3 = self._find_top3_most_intense_internal_ions(entry["annotation"]["theoretical_code"], entry["annotation"]["intensity"], pep_seq)
            spectrum["top1_internal_ion_code"].append(top_3[0][0])
            spectrum["top1_internal_seq"].append(top_3[1][0])
            spectrum["top2_internal_ion_code"].append(top_3[0][1])
            spectrum["top2_internal_seq"].append(top_3[1][1])
            spectrum["top3_internal_ion_code"].append(top_3[0][2])
            spectrum["top3_internal_seq"].append(top_3[1][2])
            spectrum["intensity_explained_aa"].append(self._find_explained_by_aa(entry["annotation"]["theoretical_code"], entry["annotation"]["intensity"]))
            spectrum["intensity_explained_precursor"].append(self._find_explained_precursor(entry))
            spectrum["spectrum_id"].append(self._get_spectrum_id(entry))
            spectrum["score"].append(self._get_identification_score(entry))
            spectrum["peptide_seq"].append(pep_seq)
            spectrum["proforma"].append(entry["proforma"])
            spectrum["peptide_length"].append(len(pep_seq))
            spectrum["nr_idents_with_same_rank"].append(entry["nr_idents_with_same_rank"])

        return [pd.DataFrame(fragment), pd.DataFrame(spectrum)]

    def _get_spectrum_id(self, entry: Dict) -> Union[int, float, str]:
        """
        Get the spectrum ID if it exists.
        """

        if "spectrum_id" in entry:
            return entry["spectrum_id"]
        else:
            return ""

    def _get_identification_score(self, entry: Dict) -> Union[int, float, str]:
        """
        Get the identification score if it exists.
        """

        if "identification_score" in entry:
            return entry["identification_score"]
        else:
            return 0.0

    def _get_ambiguity(self, entry_annotation: Dict, index: int) -> Union[int, float, str]:
        """
        Get the ambiguity if information is available.
        """

        if "matches_count" in entry_annotation:
            return entry_annotation["matches_count"][index]
        else:
            return -1

    def _find_explained_by_aa(self, fragments: List[str], intensities: List[float]) -> float:
        """
        Calculate the explained intensity by single amino acids.
        """

        intensities_single_aa = 0

        for i, frag in enumerate(fragments):
            if frag is not None:
                start, end, ion_cap_start, ion_cap_end, charge, formula = self._parse_fragment_code(frag)
                if start == end:
                    intensities_single_aa += intensities[i]

        return intensities_single_aa / sum(intensities)

    def _find_explained_precursor(self, entry: Dict) -> Union[int, float]:
        """
        Calculate the explained intensity by the precursor.
        """

        precursor_intensity = 0

        if "precursor_intensity" in entry:
            precursor_intensity = float(entry["precursor_intensity"])

        if precursor_intensity > 0:
            return precursor_intensity / sum(entry["annotation"]["intensity"])

        return -1

    def _find_top3_most_intense_internal_ions(self, fragments: List[str], intensities: List[float], pep_seq: str) -> List[List[str]]:
        """
        Get the top 3 most intense internal ions.
        """

        mapping = dict()
        for i, intensity in enumerate(intensities):
            # Add only the internal ions
            if fragments[i] is not None:
                start, end, ion_cap_start, ion_cap_end, charge, formula = self._parse_fragment_code(fragments[i])
                if "t" not in (ion_cap_start, ion_cap_end):
                    mapping[intensity] = fragments[i]

        top_3_intensities = sorted(list(mapping.keys()), reverse=True)[:3]
        top_3_codes = [mapping[intensity] for intensity in top_3_intensities]

        top_3_sequences = []

        for code in top_3_codes:
            start, end, ion_cap_start, ion_cap_end, charge, formula = self._parse_fragment_code(code)
            top_3_sequences.append(pep_seq[start-1:end])

        if len(top_3_codes) < 3:
            if len(top_3_codes) == 0:
                top_3_codes = ["", "", ""]
                top_3_sequences = ["", "", ""]
            if len(top_3_codes) == 1:
                top_3_codes.append("")
                top_3_codes.append("")
                top_3_sequences.append("")
                top_3_sequences.append("")
            else:
                top_3_codes.append("")
                top_3_sequences.append("")

        return [top_3_codes, top_3_sequences]

    def _calculate_internal_terminal_non_annotated_ions(self, fragments: List[str], intensities: List[float]) -> Dict:
        """
        Calculate ion counts and frequencies.
        """

        # Number of different ion types
        internal = 0
        terminal = 0
        non_annotated = 0

        # Total intensities of ion types
        total_int_internal = 0
        total_int_terminal = 0
        total_int_non = 0

        for i, fragment in enumerate(fragments):
            if fragment is not None:
                start, end, ion_cap_start, ion_cap_end, charge, formula = self._parse_fragment_code(fragment)
                if "t" in (ion_cap_start, ion_cap_end):
                    terminal += 1
                    total_int_terminal += intensities[i]
                else:
                    internal += 1
                    total_int_internal += intensities[i]
            else:
                non_annotated += 1
                total_int_non += intensities[i]

        return {"internal": internal / len(fragments),
                "terminal": terminal / len(fragments),
                "non_annotated": non_annotated / len(fragments),
                "total_int_internal": total_int_internal,
                "total_int_terminal": total_int_terminal,
                "total_int_non": total_int_non}

    def _parse_modfication(self, proforma: str, start: int, end: int) -> str:
        """
        Get modification string.
        """

        # Finds the modifications in the proforma string
        pattern = r"\[([\+?\-?A-Za-z0-9_\.?0-9]+)\]"
        mods = re.findall(pattern, proforma)

        positions = self._parse_modification_positions(proforma)

        modifications = ""

        # Finds the positions of the modifications
        for i, position in enumerate(positions):
            if position >= start and position <= end:
                modifications += mods[i] + ";"

        return modifications.rstrip(";")

    def _parse_modification_positions(self, seq: str) -> List[int]:
        """
        Get positions of the modifications.
        """

        i = 1
        is_modification = False
        annotation = []
        for aa in seq:
            if aa == "[":
                is_modification = True
                annotation.append(0)
            elif aa == "]":
                is_modification = False
            else:
                if not is_modification:
                    annotation.append(i)
                    i += 1

        positions = []

        for i, anno in enumerate(annotation):
            if anno == 0:
                positions.append(annotation[i-1])

        return positions

    @staticmethod
    @lru_cache(maxsize=10000)
    def _parse_fragment_code(fragment_code: str):
        """
        Parse the fragment code.
        """

        # test if fragment code format is valid*

        try:
            groups = JSONConverter.fragment_code_pattern.match(fragment_code).groupdict()
        except AttributeError:
            raise RuntimeError("Incorrect fragment code format: {0}".format(fragment_code))

        return int(groups["start"]), int(groups["end"]), groups["ion_cap_start"], groups["ion_cap_end"], int(groups["charge"]), groups["formula"] or ""
