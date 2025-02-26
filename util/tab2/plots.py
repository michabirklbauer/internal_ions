import numpy as np
import pandas as pd

import logomaker
import matplotlib.pyplot as plt
import plotly.graph_objects as go

from typing import Callable

DEFAULT_MARGIN = dict(t=0, b=0, l=0, r=0, pad=0)


def common_type_gen(fragments_dataframe: pd.DataFrame) -> pd.Series:
    common_type = fragments_dataframe["frag_type1"].astype(str).str.cat(fragments_dataframe["frag_type2"], sep="-")
    common_type = common_type.replace("n", "not annotated", regex=True)
    common_type = common_type.replace("t-", "", regex=True)
    common_type = common_type.replace("-t", "", regex=True)
    fragments_dataframe["frag_types"] = common_type
    return common_type


def common_type_hist(fragments_dataframe: pd.DataFrame) -> go.Figure:
    common_type = common_type_gen(fragments_dataframe)
    fig = go.Figure([go.Histogram(x=common_type)])
    fig.update_layout(margin=DEFAULT_MARGIN)
    return fig


def common_type_pie(fragments_dataframe: pd.DataFrame) -> go.Figure:
    common_type = common_type_gen(fragments_dataframe)
    counts = common_type.value_counts()

    pie_fig = go.Figure([go.Pie(labels=counts.keys(), values=counts)])
    pie_fig.update_layout(margin=DEFAULT_MARGIN)
    return pie_fig


def fragments_hist_plotter(colname: str, xaxis_title: str, transform: Callable = lambda x: x):
    def func(fragments_dataframe: pd.DataFrame) -> go.Figure:
        types = fragments_dataframe["frag_types"].unique()
        histograms = list()
        for t in types:
            histograms.append(go.Histogram(x=transform(fragments_dataframe.loc[fragments_dataframe.frag_types == t, colname]),
                                           histnorm="probability", name=t, nbinsx=50))

        fig = go.Figure(histograms)
        fig.update_layout(barmode="group", xaxis_title=xaxis_title, yaxis_title="Probability", margin=DEFAULT_MARGIN)
        return fig
    return func


log_ion_intens_dist = fragments_hist_plotter("frag_intensity", "Log2(Intensity)", np.log)
rel_ion_intens_perc = fragments_hist_plotter("perc_of_total_intensity", "Intensity")
rel_ion_intens_prop = fragments_hist_plotter("prop_intensity_to_base_peak", "Percentage per base peak")
mz_dist_ion_type = fragments_hist_plotter("frag_mz", "m/z")


def spectrum_hist_plotter(col_prefix, xaxis_title):
    def func(spectra_dataframe: pd.DataFrame) -> go.Figure:
        types = ["internal", "terminal", "other"]

        histograms = list()
        for t in types:
            histograms.append(go.Histogram(x=spectra_dataframe[col_prefix + t],
                                           histnorm="probability", name=t, nbinsx=50))

        fig = go.Figure(histograms)
        fig.update_layout(barmode="group", xaxis_title="Percentage", yaxis_title="Probability", margin=DEFAULT_MARGIN)

        return fig
    return func


per_spec_ion_type = spectrum_hist_plotter("perc_", "Percentage")
per_spec_ion_intens = spectrum_hist_plotter("total_int_", "Intensity Percentage")


# currently not used
def log_ion_intens_ridge(fragments_dataframe: pd.DataFrame) -> go.Figure:

    types = fragments_dataframe["frag_types"].unique()

    fig = go.Figure()
    for t in types:
        fig.add_trace(go.Violin(x=np.log(fragments_dataframe[fragments_dataframe.frag_types == t].frag_intensity), name=t))

    fig.update_traces(orientation="h", side="positive", width=3, points=False)
    fig.update_layout(barmode="group", xaxis_title="Log2(Intensity)", yaxis_title="Probability", margin=DEFAULT_MARGIN)
    return fig


def rel_ion_intens_ridge(fragments_dataframe: pd.DataFrame) -> go.Figure:

    types = fragments_dataframe["frag_types"].unique()

    fig = go.Figure()
    for t in types:
        fig.add_trace(go.Violin(x=fragments_dataframe[fragments_dataframe.frag_types == t].perc_of_total_intensity, name=t))

    fig.update_traces(orientation="h", side="positive", width=3, points=False)
    fig.update_layout(barmode="group", xaxis_title="Intensity", yaxis_title="Probability", margin=DEFAULT_MARGIN)

    return fig


# currently not used
def rel_ion_intens_prop_ridge(fragments_dataframe: pd.DataFrame) -> go.Figure:

    types = fragments_dataframe["frag_types"].unique()

    fig = go.Figure()
    for t in types:
        fig.add_trace(go.Violin(x=fragments_dataframe[fragments_dataframe.frag_types == t].prop_intensity_to_base_peak, name=t))

    fig.update_traces(orientation="h", side="positive", width=3, points=False)
    fig.update_layout(barmode="group", xaxis_title="Intensity", yaxis_title="Probability", margin=DEFAULT_MARGIN)

    return fig


def logo_of_fraction(spectra_dataframe: pd.DataFrame,
                     fragments_dataframe: pd.DataFrame,
                     topn: int = 3,
                     max_length: int = 10,
                     min_length: int = 10) -> plt.figure:

    fig1, ax1 = plt.subplots(facecolor="white", figsize=(8, 4), dpi=1000)

    if topn == 0:
        df = fragments_dataframe.frag_seq
    else:
        df = spectra_dataframe.top1_internal_seq
    if topn > 1:
        df = pd.concat([df, spectra_dataframe.top2_internal_seq])
    if topn > 2:
        df = pd.concat([df, spectra_dataframe.top3_internal_seq])

    # Finding the maximum length of the strings in the column
    # Creating a new dataframe to store the results
    df_frequency = pd.DataFrame()

    # Iterating over each position up to the maximum length
    tt = df[(df.str.len() <= max_length) & (df.str.len() >= min_length)]
    for i in range(0, max_length):
        # Counting the frequency of each letter at position i
        frequency = tt.str[i].value_counts()
        # Adding the frequency to the new dataframe
        df_frequency[i] = frequency / sum(frequency)

    # Displaying the frequency dataframe
    df_frequency = df_frequency.fillna(0)
    # Transposing and plotting logos
    ww_df = df_frequency.T

    if ww_df.shape[0] == 0 or ww_df.shape[1] == 0:
        return fig1

    logomaker.Logo(ww_df, ax=ax1, color_scheme="NajafabadiEtAl2017", vpad=0.1, width=0.8, show_spines=False)
    ax1.set_xlabel("Position in subset")
    ax1.set_ylabel("Amino acid frequency")
    ax1.set_facecolor("white")

    return fig1
