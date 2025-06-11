import pandas as pd
from typing import List


def filter_dataframes(dataframes: List[pd.DataFrame],
                      N_ion: bool) -> List[pd.DataFrame]:
    """
    Filters non-annotated ions
    """

    fragments_df = dataframes[0]
    spectra_df = dataframes[1]

    if not N_ion:
        fragments_df = fragments_df.loc[(fragments_df.frag_type1 != "n") & (fragments_df.frag_type2 != "n")].copy()

    return [fragments_df, spectra_df]
