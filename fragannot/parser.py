from psm_utils.psm_list import PSMList
from typing import BinaryIO
import logging
import logging.config
import os
from collections import Counter
from .spectrumfile import SpectrumFile

RANK_LIMIT = 1  # maximum PSM rank allowed
HIT_LIMIT = 1  # actual max number of PSMs per spectrum allowed (sometimes there are many PSMs with the same rank)


class Parser:
    def __init__(self, psm_list: PSMList, is_streamlit: bool = False):
        # Set up logging
        self.logger = logging.getLogger(__name__)
        self.spectra = None
        self.psm_list = psm_list
        self.is_streamlit = is_streamlit

    def read(self, raw_file: BinaryIO, max_rank=RANK_LIMIT, max_hits=HIT_LIMIT):
        self.__load(raw_file)
        self.logger.info(f"Read {len(self.psm_list)} PSMs from identification file")
        if self.is_streamlit:
            print(f"Read {len(self.psm_list)} PSMs from identification file")

        count = 0
        spec_counts = Counter()
        for psm in self.psm_list:
            if (psm.rank or 1) <= max_rank:
                try:
                    spectrum = self.spectra.get_by_id(psm["spectrum_id"])
                except KeyError:
                    self.logger.warning(f'SpectrumId - {psm["spectrum_id"]} not found')
                else:
                    if spec_counts.get(psm.spectrum_id, 0) < max_hits:
                        spec_counts[psm.spectrum_id] += 1
                        psm.spectrum = {"mz": spectrum["m/z array"], "intensity": spectrum["intensity array"]}
                        count += 1
                        if count % 500 == 0:
                            self.logger.info(f"{count} spectra processed")
                            if self.is_streamlit:
                                print(f"{count} spectra processed")

        output_fpath = os.path.splitext(raw_file.name)[0] + ".json"
        self.output_fname = os.path.basename(output_fpath)
        return self.psm_list

    def __load(self, raw_file: BinaryIO):
        self.spectra = self.__read_raw_file(raw_file)

    def __read_raw_file(self, file):
        return SpectrumFile(file)
