#!/usr/bin/env python3

import pandas as pd
from typing import List


def filter_dataframes(dataframes: List[pd.DataFrame],
                      ion_filter: List[bool]) -> List[pd.DataFrame]:
    """
    Filters dataframes based on the given values.
    """

    ion_filter_translation = ["n", "a", "b", "c", "cdot", "c-1", "c+1", "x", "y", "z", "zdot", "z+1", "z+2", "z+3"]
    ions_considered = []

    for i, val in enumerate(ion_filter):
        if val:
            ions_considered.append(ion_filter_translation[i])

    fragments_df = dataframes[0].copy()
    spectra_df = dataframes[1].copy()

    # filter by ion type
    fragments_df = fragments_df[fragments_df["frag_code"].str.contains("|".join(ions_considered))]

    return [fragments_df, spectra_df]
