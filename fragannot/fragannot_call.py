from .fragannot_numba import FragannotNumba

from typing import Dict
from typing import List
from typing import BinaryIO


def fragannot_call(spectrum_file: BinaryIO,
                   identifications_file: BinaryIO,
                   tolerance: float,
                   fragment_types: List[str],
                   charges: List[str],
                   losses: List[str],
                   deisotope: bool,
                   file_format: str = "infer",
                   verbose: bool = False) -> Dict:

    frag = FragannotNumba()
    fragannot_dict = frag.fragment_annotation(identifications_file,
                                              spectrum_file,
                                              tolerance,
                                              fragment_types,
                                              charges,
                                              losses,
                                              file_format,
                                              deisotope,
                                              write_file=False)

    return fragannot_dict
