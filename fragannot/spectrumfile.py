import logging
from pyteomics import mgf
from os.path import splitext
from collections import Counter
from io import StringIO
from typing import TextIO


class SpectrumFile:
    def __init__(self, uploaded_file):
        self.logger = logging.getLogger(__name__)
        self.file_format = None
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
        use_scans = len(scans) > len(titles)
        if max(len(scans), len(titles)) < i + 1:
            self.logger.warning("Not all spectra can be indexed in %s: %d scans, %d titles", file, len(scans), len(titles))
        return use_scans

    def _load(self, uploaded_file):
        uploaded_file.seek(0)
        extension = splitext(uploaded_file.name)[1]

        if extension.lower() == ".mgf":
            self.logger.info(f"Inferred MGF format from {uploaded_file.name}")
            index_by_scans = self.get_mgf_index_by_scans(StringIO(uploaded_file.getvalue().decode('utf-8')))
            self.spectra_source = mgf.IndexedMGF(uploaded_file, index_by_scans=index_by_scans)
            self.file_format = 'mgf'
            self.index = self.spectra_source.index

        else:
            self.logger.info(
                f"Cannot infer format from {uploaded_file}, only mzML and MGF formats are supported"
            )
            raise Exception("Unsupported spectra file format")

        self.logger.info(f"Loaded a SpectrumFile with {len(self.spectra_source)} spectra.")

    def get_by_id(self, i):
        return self.spectra_source[i]
