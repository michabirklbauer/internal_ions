#!/usr/bin/env python3

from pyteomics.proforma import ProForma

from typing import List

def read_proforma(proforma_string: str) -> List[int]:

    result = []
    p = ProForma.parse(proforma_string)

    if not p.n_term is None:
        for mod in p.n_term:
            result.append(0)

    for i, (aa, mods) in enumerate(p.sequence):
        if not mods is None:
            for mod in mods:
                result.append(i + 1)

    if not p.c_term is None:
        for mod in p.c_term:
            result.append(-1)

    return result
