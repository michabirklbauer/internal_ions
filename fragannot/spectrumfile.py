import re
import logging
from pyteomics import mzml, mgf
from os.path import splitext
from collections import defaultdict

import sys
sys.path.append("..")
from util.spectrumio import parse_scannr


class SpectrumFile:
    def __init__(self, file_path, parser_pattern=r'index(?:\s+)?=(?:\s+)?(\d+)'):
        self.logger = logging.getLogger(__name__)
        self.indices = defaultdict(dict)
        self.spectra_source = None
        self.file_format = None
        self._build_index = None
        self.get_by_id = None
        self.parser_pattern = parser_pattern

        self._load(file_path)
        self._build_index()

    def __repr_int(self, var: str) -> bool:
        if type(var) is not str:
            return False
        try:
            int(var)
        except ValueError:
            return False
        return True

    def _load(self, file_path):
        extension = splitext(file_path)[1]

        if extension.lower() == ".mzml":
            self.logger.info(f"Inferred mzML format from {file_path}")
            self.spectra_source = mzml.MzML(file_path)
            self.file_format = 'mzml'
            self._build_index = self._index_MZML
            self.get_by_id = self._get_by_id_MZML

        elif extension.lower() == ".mgf":
            self.logger.info(f"Inferred MGF format from {file_path}")
            self.spectra_source = mgf.IndexedMGF(file_path)
            self.file_format = 'mgf'
            self._build_index = self._index_MGF
            self.get_by_id = self._get_by_id_MGF

        else:
            self.logger.info(
                f"Cannot infer format from {file_path}, only mzML and MGF formats are supported"
            )
            raise Exception("Unsupported spectra file format")

    def _get_by_id_MGF(self, id_string):
        if type(id_string) is int or self.__repr_int(id_string):
            if type(id_string) is str:
                id_string = int(id_string)
            return self.spectra_source.get_by_index(self.indices['scan'][id_string])
        elif type(id_string) is str:
            # in case of whitespaces
            id_string = id_string.strip()
            try:
                # the default way
                return self.spectra_source.get_by_id(id_string)
            except KeyError:
                # trying to recover
                scan_nr = parse_scannr(id_string, -1, self.parser_pattern)
                if scan_nr[0] == 0:
                    return self.spectra_source.get_by_index(self.indices['scan'][scan_nr[1]])

                match_index = re.match(r'index(?:\s+)?=(?:\s+)?(\d+)', id_string)

                if match_index is not None:
                    return self.spectra_source.get_by_index(int(match_index.group(1)))

                if id_string in self.indices['title'].keys():
                    return self.spectra_source.get_by_index(self.indices['title'][id_string])

                if type(id_string) is int or id_string.isdigit():
                    return self.spectra_source.get_by_index(self.indices['scan'][int(id_string)])

                raise KeyError(f'Cannot infer MGF scan from spectrum_id={id_string}')

        else:
            raise TypeError(f'Unsupported id type ({type(id_string)}): should be int or string')

    def _get_by_id_MZML(self, id_string):
        if type(id_string) is int or self.__repr_int(id_string):
            if type(id_string) is str:
                id_string = int(id_string)
            return self.spectra_source.get_by_index(self.indices['scan'][id_string])
        elif type(id_string) is str:
            # in case of whitespaces
            id_string = id_string.strip()
            try:
                # the default way
                return self.spectra_source.get_by_id(id_string)
            except KeyError:
                # trying to recover
                scan_nr = parse_scannr(id_string, -1, self.parser_pattern)
                if scan_nr[0] == 0:
                    return self.spectra_source.get_by_index(self.indices['scan'][scan_nr[1]])

                match_index = re.match(r'index(?:\s+)?=(?:\s+)?(\d+)', id_string)

                if match_index is not None:
                    return self.spectra_source.get_by_index(int(match_index.group(1)))

                # if id_string in self.indices['title'].keys():
                #     return self.spectra_source.get_by_index(self.indices['title'][id_string])

                if type(id_string) is int or id_string.isdigit():
                    return self.spectra_source.get_by_index(self.indices['scan'][int(id_string)])

                raise KeyError(f'Cannot infer MzML scan from spectrum_id={id_string}')

        else:
            raise TypeError(f'Unsupported id type ({type(id_string)}): should be int or string')

    def _index_MGF(self):
        for index in range(len(self.spectra_source)):
            params = self.spectra_source[index]['params']
            if 'title' in params.keys():
                self.indices['title'][params['title']] = index

            if 'scans' in params.keys():
                self.indices['scan'][int(params['scans'])] = index

    def _index_MZML(self):
        for spectrum in self.spectra_source:
            spectrumID = spectrum['id']
            index = spectrum['index']

            scan_match = re.match(r'.+scan(?:\s+)?=(?:\s+)?(\d+)', spectrumID)
            if scan_match is not None:
                self.indices['scan'][int(scan_match.group(1))] = index

            self.indices['id'][spectrumID] = index
