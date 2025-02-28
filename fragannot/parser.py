from psm_utils.io import read_file
from psm_utils import Peptidoform, PSM
import traceback
import logging
import logging.config
import json
import requests
from pyteomics import mzid
from urllib.parse import urlparse
from os.path import splitext, split
import os
import re
import gzip
from .spectrumfile import SpectrumFile

RANK_LIMIT = 1


# Excerpt from MS:1001143 items (PSM-level search engine specific statistic)
STANDARD_SEARCHENGINE_SCORES = [
    "Amanda:AmandaScore",
    "Andromeda:score",
    "Byonic:Score",
    "Comet:Xcorr",
    "Comet:xcorr",
    "DeBunker:score",
    "IdentityE Score",
    "KSDP score",
    "MS-GF:RawScore",
    "MSFit:Mowse score",
    "MSPathFinder:RawScore",
    "MSPepSearch:score",
    "Mascot:score",
    "MetaMorpheus:score",
    "OMSSA:evalue",
    "OpenPepXL:score",
    "PEAKS:peptideScore",
    "PeptideShaker PSM score",
    "Phenyx:Pepzscore",
    "ProLuCID:xcorr",
    "ProSight:specral C-score",
    "Profound:z value",
    "ProteinProspector:score",
    "ProteinScape:SequestMetaScore",
    "ProteomeDiscoverer:Delta Score",
    "SEQUEST:xcorr",
    "SIM-XL score ",
    "SQID:score ",
    "Sonar:Score",
    "SpectrumMill:Score",
    "TopMG:spectral E-Value",
    "X!Tandem:hyperscore",
    "ZCore:probScore:",
    "Percolator:score",
    "xi:score",
]


class Parser:
    def __init__(self, parser_pattern: str, is_streamlit: bool = False):
        # Set up logging
        self.logger = logging.getLogger(__name__)
        self.spectra = None
        self.psm_list = None
        self.is_streamlit = is_streamlit
        self.parser_pattern = parser_pattern

    def read(self, raw_file, ident_file, file_format, max_rank=RANK_LIMIT):
        """Read and process raw file and identification file.

        Parameters
        ----------
        raw_file : str
            Path or url to the raw file
        ident_file_path : str
            Path or url to the identification file
        """

        self.max_rank = max_rank
        try:
            if self.__is_url(raw_file):
                self.logger.info("Raw file is not local, try to download.")
                raw_file_name = self.__get_file_from_url(raw_file)
            else:
                if os.path.exists(raw_file):
                    raw_file_name = raw_file
                else:
                    self.logger.error(f"File doesn't exist: {raw_file}")
                    raise Exception("File doesn't exist")
            if self.__is_url(ident_file):
                self.logger.info("Ident file is not local, try to download.")
                ident_file_name = self.__get_file_from_url(ident_file)
            else:
                if os.path.exists(ident_file):
                    ident_file_name = ident_file
                else:
                    self.logger.error(f"File doesn't exist: {raw_file}")
                    raise Exception("File doesn't exist")
            self.__load(raw_file_name, ident_file_name, file_format)
        except Exception:
            self.logger.error(f"Couldn't read file. Exception:\n{traceback.format_exc()}")

        self.logger.info(f"Read {len(self.psm_list)} PSMs from identification file")
        if self.is_streamlit:
            print(f"Read {len(self.psm_list)} PSMs from identification file")

        count = 0
        for psm in self.psm_list:
            try:
                spectrum = self.spectra.get_by_id(psm["spectrum_id"])
                psm.spectrum = {"mz": spectrum["m/z array"], "intensity": spectrum["intensity array"]}
                count += 1
                if count % 500 == 0:
                    self.logger.info(f"{count} spectra processed")
                    if self.is_streamlit:
                        print(f"{count} spectra processed")

            except KeyError:
                self.logger.warning(f'SpectrumId - {psm["spectrum_id"]} not found')
        output_fpath = os.path.splitext(raw_file_name)[0] + ".json"
        self.output_fname = os.path.basename(output_fpath)
        return self.psm_list

    def __load(self, raw_file_path, ident_file_path, file_format):
        """Load raw file and identification file.

        Parameters
        ----------
        raw_file_path : str
            Path to the raw file
        ident_file_path : str
            Path to the identification file
        """
        self.spectra = self.__read_raw_file(raw_file_path)
        self.psm_list = self.__read_id_file(ident_file_path, file_format)

    def __read_raw_file(self, file_path):
        """Read raw file

        Parameters
        ----------
        file_path : str
            Path to the raw file
        """
        return SpectrumFile(file_path, self.parser_pattern)

    def __read_id_file(self, file_path, file_format):
        """Read identification file more generously then psm_utils

        Parameters
        ----------
        file_path : str
            Path to the raw file
        file_format : str
            Identification file format
        """
        extension = splitext(file_path)[1]

        if extension.lower() == ".mzid" or extension.lower() == ".mzidentml":
            result = []

            score_name = ""
            for psm in mzid.MzIdentML(file_path):

                spectrumID = psm["spectrumID"]
                for spectrum_identification in psm["SpectrumIdentificationItem"]:
                    rank = spectrum_identification["rank"]
                    charge_state = spectrum_identification["chargeState"]
                    if rank <= self.max_rank:
                        if score_name == "":
                            score_name = self.__infer_score_name(spectrum_identification.keys())
                        score = spectrum_identification[score_name]
                        sequence = spectrum_identification["PeptideEvidenceRef"][0]["PeptideSequence"]
                        modifications = spectrum_identification["PeptideEvidenceRef"][0].get(
                            "Modification", None
                        )

                        filename = split(psm["location"])[1]
                        if modifications is not None:
                            aas = [""] + [aa for aa in sequence] + [""]
                            for mod in modifications:
                                loc = mod["location"]
                                mass = mod["monoisotopicMassDelta"]
                                if "residues" in mod.keys():
                                    res = mod["residues"][0]
                                    if loc > 0 and not res == aas[loc]:
                                        raise Exception(
                                            f"Mismatch {modifications} {sequence} {res} {aas[loc]} {loc}"
                                        )
                                try:
                                    aas[loc] += f"[+{mass}]" if mass > 0 else f"[{mass}]"
                                except Exception:
                                    self.logger.error(f"Mass error at psm: {psm}")
                                    break

                            sequence = "".join(aas[1:-1])

                            if aas[0] != "":
                                sequence = f"{aas[0]}-{sequence}"

                            if aas[-1] != "":
                                sequence = f"{sequence}-{aas[-1]}"

                        # add charge state of precursor
                        sequence = sequence + "/" + str(charge_state)

                        result.append(
                            PSM(
                                peptidoform=Peptidoform(sequence),
                                run=filename,
                                spectrum_id=spectrumID,
                                score=score,
                                rank=rank,
                                charge=charge_state,
                            )
                        )

            return result

        else:
            return read_file(file_path, filetype=file_format)

    def __load_log_config(self):
        """Load log configurations from log_conf.json

        Returns:
        -------
        config : Configuration loaded fro, log_conf.json
        """
        config = {}
        with open("log_conf.json", "r", encoding="utf-8") as fd:
            config = json.load(fd)
        return config

    def __uncompress(self, file_path, block_size=65536):
        new_file_path = file_path
        if file_path.endswith(".gz"):
            new_file_path = file_path[:-3]
            with gzip.open(file_path, "rb") as s_file, open(new_file_path, "wb") as d_file:
                while True:
                    block = s_file.read(block_size)
                    if not block:
                        break
                    else:
                        d_file.write(block)
        return new_file_path

    # Check if the given filename is a url or not
    def __is_url(self, url):
        is_url = True
        url_parsed = urlparse(url)
        if url_parsed.scheme in ("file", ""):
            is_url = False
        return is_url

    def __get_file_from_url(self, url):
        """Download file from url and return the path to the file

        Parameters
        ----------
        url : str
            URL to the file to download
        """
        if not os.path.exists(r".\downloads"):
            os.makedirs(r".\downloads")

        r = requests.get(url, allow_redirects=True)
        if url.find("/"):
            file_path = ".\\downloads\\" + url.rsplit("/", 1)[1]
        else:
            file_path = ".\\downloads\\" + self.__getFilename_fromCd(r.headers.get("content-disposition"))

        open(file_path, "wb").write(r.content)
        file_path = self.__uncompress(file_path)

        return file_path

    def __getFilename_fromCd(self, cd):
        """Gets filename from content disposition

        Parameters
        ----------
        cd : str
            Content disposition
        """
        if not cd:
            return None
        fname = re.findall("filename=(.+)", cd)
        self.logger.info(fname)
        if len(fname) == 0:
            return None
        return fname[0]

    def __infer_score_name(self, keys):
        """Infer the score from the known list of PSM scores."""
        all_scores = []
        for score in STANDARD_SEARCHENGINE_SCORES:
            if score in keys:
                all_scores.append(score)
        if len(all_scores) != 0:
            if len(all_scores) > 1:
                self.logger.warning(
                    "More than one score was found in the identification"
                    f"file {all_scores}, used score: {all_scores[0]}"
                )
            return all_scores[0]
        else:
            self.logger.info(keys)
            raise Exception("No known score metric found in mzIdentML file.")
