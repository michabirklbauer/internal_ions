from .fragannot_numba import FragannotNumba

from typing import Dict
from typing import List
from typing import BinaryIO
from psm_utils.psm_list import PSMList


def fragannot_call(spectrum_file: BinaryIO,
                   psms: PSMList,
                   tolerance: float,
                   fragment_types: List[str],
                   charges: List[str],
                   losses: List[str],
                   deisotope: bool,
                   verbose: bool = False) -> Dict:

    frag = FragannotNumba()
    fragannot_dict = frag.fragment_annotation(psms,
                                              spectrum_file,
                                              tolerance,
                                              fragment_types,
                                              charges,
                                              losses,
                                              deisotope,
                                              write_file=False)

    return fragannot_dict
