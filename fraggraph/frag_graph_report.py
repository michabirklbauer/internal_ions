# Standard Library Imports
import json
import math
import itertools
from itertools import combinations, chain

# Third-party Library Imports
import numpy as np
import pandas as pd
from scipy.special import rel_entr, softmax
from numpy.linalg import norm
from sklearn.mixture import GaussianMixture
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
from sklearn.metrics.pairwise import cosine_similarity
import plotly.graph_objects as go
from tqdm import tqdm
import networkx as nx
from pyteomics import mass, parser as pyteomics_parser
import brainpy as bp
import fraggraph.constant as constant
import psm_utils
from scipy.optimize import minimize
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from plotnine import (
    ggplot,
    aes,
    geom_tile,
    geom_jitter,
    geom_hline,
    scale_fill_manual,
    position_stack,
    geom_text,
    geom_abline,
    geom_smooth,
    geom_bar,
    geom_boxplot,
    geom_violin,
    geom_point,
    geom_rect,
    scale_fill_gradientn,
    labs,
    theme_bw,
    theme,
    element_blank,
    element_line,
    element_text,
    annotate,
    scale_x_reverse,
    scale_y_continuous,
    ggtitle,
    ggsave,
    ylim,
    xlim,
)
import statsmodels.api as sm
import statsmodels.formula.api as smf
import statsmodels.tools as sm_tools
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.colors as matcolors

""" script to generate visualizations in ggplot2 from a FragGraph instance """


def draw_fragment_coverage_matrix(
    FG, x="intensity", x_min=1, x_max=10, filename=None, peptidoform_index=0
):
    """
    This function draws a fragment coverage matrix using ggplot.

    Parameters:
    FG: Fragment Graph object
    x (str): The column name to use for the intensity. Default is "intensity".
    x_min (int): The minimum limit for the color scale. Default is 1.
    x_max (int): The maximum limit for the color scale. Default is 10.
    filename (str): The filename to save the plot. If None, the plot is returned instead of being saved. Default is None.
    peptidoform_index (int): The index of the peptidoform to use. Default is 0.

    Returns:
    ggplot object if filename is None, else None.
    """
    try:
        # Get fragment table I0 from FG (assuming it returns a DataFrame)
        frag_df = FG.get_fragment_table_I0(peptidoform_index)

        # Data processing
        frag_df["intensity"] = np.log(frag_df[x]).replace(-np.inf, np.nan)
        frag_df[["start_pos", "end_pos"]] = frag_df[["start_pos", "end_pos"]].astype(
            int
        )

        # Define reverse viridis color palette
        colors = [
            "#fde725",
            "#b5de2b",
            "#6ece58",
            "#35b779",
            "#1f9e89",
            "#26828e",
            "#31688e",
            "#3e4989",
            "#482878",
            "#440154",
        ]

        # Create ggplot
        gg = (
            ggplot(frag_df, aes(x="end_pos", y="start_pos", fill=x))
            + geom_tile(colour="black", size=0.1)
            + scale_fill_gradientn(
                colors=colors, na_value="grey", limits=(x_min, x_max)
            )
            + labs(
                x="End Position (AA index)",
                y="Start Position (AA index)",
                fill=x,
            )
            + theme_bw()
            + theme(
                panel_grid=element_blank(),
                panel_border=element_blank(),
                axis_text=element_text(size=5),
                axis_title=element_text(size=8),
                legend_text=element_text(size=10),
                legend_title=element_text(size=12),
                axis_text_x=element_text(angle=45, hjust=1),
            )
            + scale_x_reverse(
                breaks=range(
                    int(min(frag_df["start_pos"])), int(max(frag_df["end_pos"])) + 1
                )
            )
            + scale_y_continuous(
                breaks=range(
                    int(min(frag_df["start_pos"])), int(max(frag_df["end_pos"])) + 1
                )
            )
        )

        # Save to file or return plot
        if filename:
            gg.save(filename=filename)
            print(f"Plot saved to {filename}")

        return gg

    except Exception as e:
        print(f"An error occurred: {e}")


def draw_fragment_coverage_matrix_binned_mz(
    FG,
    mz_bin_size=0.01,
    filename=None,
    peptidoform_index=0,
):
    """
    Generate a fragment matrix where a color is attributed for each mz bin.

    Parameters:
    - FG: FragGraph instance
    - mz_bin_size: float, optional (default=0.01)
        The size of the mz bin.
    - filename: str, optional (default=None)
        The filename to save the plot. If None, the plot will be displayed directly.
    - peptidoform_index: int, optional (default=0)
        The index of the peptidoform to use.
    """

    # colors to cycle through
    colors = [
        "red",
        "blue",
        "green",
        "purple",
        "orange",
        "yellow",
        "pink",
        "brown",
        "grey",
        "black",
    ]

    # get fragment table I0 from FG (assuming it returns a DataFrame)
    frag_df = FG.get_fragment_table_I0(peptidoform_index=peptidoform_index)

    # log the intensity column
    frag_df["intensity"] = np.log(frag_df["intensity"])
    # make sure that start_pos and end_pos are integers
    frag_df["start_pos"] = frag_df["start_pos"].astype(int)
    frag_df["end_pos"] = frag_df["end_pos"].astype(int)

    # create color column
    # iterate over mz values in ascending order
    # if the difference from one mz value to the other is bigger than the bin size, change color
    # if the difference is smaller, keep the same color

    # sort the dataframe by mz
    frag_df = frag_df.sort_values(by=["mz"])
    # reset index
    frag_df = frag_df.reset_index(drop=True)
    # iterate over mz values
    color = ["red"]
    current_color = 0
    bin = [0]
    current_bin = 0
    for i in range(0, len(frag_df) - 1):
        if abs(frag_df["mz"][i + 1] - frag_df["mz"][i]) > mz_bin_size:
            current_color += 1
            current_bin += 1
            if current_color == len(colors):
                current_color = 0
        color.append(colors[current_color])
        bin.append(current_bin)
    # add color column to dataframe
    frag_df["color"] = color
    frag_df["bin"] = bin

    # add a column with the number of fragment in each bin
    frag_df["n_frag"] = frag_df.groupby("bin")["bin"].transform("count")
    # add a column bin_str for the labels, if n_frag > 1, add the number of fragment in the bin
    frag_df["bin_str"] = frag_df["bin"].astype(str)

    # replace bin_str by "." if n_frag == 1
    frag_df.loc[frag_df["n_frag"] == 1, "bin_str"] = "."

    # plot
    gg = (
        ggplot(frag_df, aes(x="end_pos", y="start_pos", fill="mz"))
        + geom_tile(colour="black", size=0.1)
        + labs(x="End Pos", y="Start Pos", fill="Total Intensity")
        + theme_bw()
        + labs(x="End Position (AA index)", y="Start Position (AA index)")
        + theme(
            panel_grid=element_blank(),
            panel_border=element_blank(),
            axis_text=element_text(size=5),
            axis_title=element_text(size=8),
            legend_text=element_text(size=10),
            legend_title=element_text(size=12),
            axis_text_x=element_text(angle=45, hjust=1),
        )
        + scale_x_reverse(
            breaks=range(
                int(min(frag_df["start_pos"])), int(max(frag_df["end_pos"])) + 1
            )
        )
        + scale_y_continuous(
            breaks=range(
                int(min(frag_df["start_pos"])), int(max(frag_df["end_pos"])) + 1
            )
        )
        + theme(legend_position="none")
        + geom_text(aes(label="bin_str"), size=2, color="white")
        + theme(legend_position="none")
    )

    # save or display the plot
    if filename is not None:
        gg.save(filename)
    else:
        gg.draw()

    return gg

def draw_fragment_coverage_matrix_difference(
    FG1,
    FG2,
    x="intensity",
    filename=None,
    mod1_position=None,
    mod2_position=None,
    cell_width=1,
    cell_height=1,
    FG1_peptidoform_index=0,
    FG2_peptidoform_index=0,
):
    try:
        # Extracting fragment data
        frag_df_1 = FG1.get_fragment_table_I0(peptidoform_index=FG1_peptidoform_index)
        frag_df_2 = FG2.get_fragment_table_I0(peptidoform_index=FG2_peptidoform_index)

        # Data type conversion
        frag_df_1["start_pos"] = frag_df_1["start_pos"].astype(int)
        frag_df_1["end_pos"] = frag_df_1["end_pos"].astype(int)
        frag_df_2["start_pos"] = frag_df_2["start_pos"].astype(int)
        frag_df_2["end_pos"] = frag_df_2["end_pos"].astype(int)

        # Merging dataframes
        frag_df = pd.merge(
            frag_df_1,
            frag_df_2,
            on=["start_pos", "end_pos"],
            how="outer",
            suffixes=("_1", "_2"),
        )

        # Handling missing values
        frag_df["intensity_1"] = frag_df["intensity_1"].fillna(1)
        frag_df["intensity_2"] = frag_df["intensity_2"].fillna(1)

        # Computing log2 fold change
        frag_df["FC"] = np.log2(frag_df["intensity_2"] / frag_df["intensity_1"]).clip(
            -10, 10
        )

        # Special cases
        frag_df.loc[
            (frag_df["intensity_1"] > 1) & (frag_df["intensity_2"] == 1), "FC"
        ] = -10
        frag_df.loc[
            (frag_df["intensity_2"] > 1) & (frag_df["intensity_1"] == 1), "FC"
        ] = 10

        # Annotations
        frag_df["annotation"] = frag_df["FC"].apply(
            lambda x: "✖" if x == -10 else ("●" if x == 10 else "")
        )

        # Building the plot
        gg = (
            ggplot(frag_df, aes(x="end_pos", y="start_pos", fill="FC"))
            + geom_tile(colour="black", size=0.1)
            + geom_text(aes(label="annotation"), size=4, color="white")
            + scale_fill_gradientn(
                colors=["#2b58ba", "#d6d6d6", "#d93838"],
                na_value="grey",
                limits=(-10, 10),
            )
            + labs(
                x="End Position (AA index)",
                y="Start Position (AA index)",
                fill="Log2FC Intensity",
            )
            + theme_bw()
            + theme(
                panel_grid=element_blank(),
                panel_border=element_blank(),
                axis_text=element_text(size=5),
                axis_title=element_text(size=8),
                legend_text=element_text(size=10),
                legend_title=element_text(size=12),
                axis_text_x=element_text(angle=45, hjust=1),
            )
            + scale_x_reverse(
                breaks=range(
                    int(min(frag_df["start_pos"])), int(max(frag_df["end_pos"])) + 1
                )
            )
            + scale_y_continuous(
                breaks=range(
                    int(min(frag_df["start_pos"])), int(max(frag_df["end_pos"])) + 1
                )
            )
        )

        # Add annotations if positions are specified
        if mod1_position is not None and mod2_position is not None:
            # Calculate positions for rectangles
            first_mod = min(mod1_position, mod2_position)
            second_mod = max(mod1_position, mod2_position) - 1

            C_term_rect, N_term_rect = calculate_rect_positions(
                first_mod, second_mod, cell_width, cell_height, frag_df
            )

            (
                gg
                + annotate(
                    "rect", **C_term_rect, fill="#FF000000", color="purple", size=0.8
                )
                + annotate(
                    "rect", **N_term_rect, fill="#FF000000", color="purple", size=0.81
                )
            )

        # Saving or returning the plot
        if filename:
            gg.save(filename)
            print(f"Plot saved to {filename}")
            return gg
        else:
            return gg

    except Exception as e:
        print(f"An error occurred: {e}")


def calculate_rect_positions(first_mod, second_mod, cell_width, cell_height, frag_df):
    C_term_rect = {
        "xmin": second_mod + cell_height / 2,
        "xmax": first_mod - cell_height / 2,
        "ymin": min(frag_df["start_pos"]) - cell_width / 2,
        "ymax": first_mod - cell_width / 2,
    }

    N_term_rect = {
        "xmin": max(frag_df["end_pos"]) + cell_width / 2,
        "xmax": second_mod + cell_height / 2,
        "ymin": first_mod - cell_height / 2 + cell_height,
        "ymax": second_mod + cell_width / 2 + cell_height,
    }

    return C_term_rect, N_term_rect


def barplot_intensity_compare(FG1, FG2, ions_type="terminal", filename=None):
    try:
        # Extracting fragment data
        frag_df_1 = FG1.get_fragment_table_I0()
        frag_df_2 = FG2.get_fragment_table_I0()
        max_end_pos = max(max(frag_df_1["end_pos"]), max(frag_df_2["end_pos"]))

        # Selecting ions type
        if ions_type == "terminal":
            frag_df_1 = frag_df_1[
                (frag_df_1["start_pos"] == 1) | (frag_df_1["end_pos"] == max_end_pos)
            ]
            frag_df_2 = frag_df_2[
                (frag_df_2["start_pos"] == 1) | (frag_df_2["end_pos"] == max_end_pos)
            ]
        elif ions_type == "internal":
            frag_df_1 = frag_df_1[
                (frag_df_1["start_pos"] != 1) & (frag_df_1["end_pos"] != max_end_pos)
            ]
            frag_df_2 = frag_df_2[
                (frag_df_2["start_pos"] != 1) & (frag_df_2["end_pos"] != max_end_pos)
            ]

        # Merging dataframes
        frag_df = pd.merge(
            frag_df_1,
            frag_df_2,
            on=["start_pos", "end_pos"],
            how="outer",
            suffixes=("_1", "_2"),
        )

        # Convert to long format
        frag_df = pd.melt(
            frag_df,
            id_vars=["start_pos", "end_pos"],
            value_vars=["intensity_1", "intensity_2"],
            var_name="group",
            value_name="intensity",
        )

        # Preparing the plot
        gg = (
            ggplot(frag_df, aes(x="group", y="intensity", fill="group"))
            + geom_boxplot()
            + labs(
                x="Group",
                y="Intensity",
                title=f"Intensity Comparison ({ions_type} ions)",
            )
            + theme_bw()
            + theme(
                axis_text_x=element_text(angle=45, hjust=1),
                axis_title=element_text(size=10),
                legend_title=element_text(size=12),
            )
        )

        # Saving or returning the plot
        output_filename = f"{filename}_{ions_type}.pdf" if filename else None
        if output_filename:
            gg.save(output_filename)
            print(f"Plot saved to {output_filename}")
        else:
            return gg

    except Exception as e:
        print(f"An error occurred: {e}")


def barplot_intensity_compare_site_determining(
    FG1,
    FG2,
    mod1_position,
    mod2_position,
    ions_type="terminal",
    filename=None,
    FG1_peptidoform_index=0,
    FG2_peptidoform_index=0,
):
    try:
        frag_df_1 = FG1.get_fragment_table_I0(peptidoform_index=FG1_peptidoform_index)
        frag_df_2 = FG2.get_fragment_table_I0(peptidoform_index=FG2_peptidoform_index)

        max_end_pos = max(max(frag_df_1["end_pos"]), max(frag_df_2["end_pos"]))

        # Merging and filtering dataframes
        frag_df = pd.merge(
            frag_df_1,
            frag_df_2,
            on=["start_pos", "end_pos"],
            how="outer",
            suffixes=("_1", "_2"),
        )

        frag_df = pd.melt(
            frag_df,
            id_vars=["start_pos", "end_pos"],
            value_vars=["intensity_1", "intensity_2"],
            var_name="group",
            value_name="intensity",
        )

        # Selecting ions based on type
        if ions_type == "terminal":
            frag_df = frag_df[
                (frag_df["start_pos"] == 1) | (frag_df["end_pos"] == max_end_pos)
            ]
        elif ions_type == "internal":
            frag_df = frag_df[
                (frag_df["start_pos"] != 1) & (frag_df["end_pos"] != max_end_pos)
            ]

        # Filtering site determining fragments
        frag_df_n = frag_df[
            (frag_df["start_pos"] < mod1_position)
            & (frag_df["end_pos"] >= mod1_position)
            & (frag_df["end_pos"] < mod2_position)
        ]
        frag_df_c = frag_df[
            (frag_df["end_pos"] > mod2_position)
            & (frag_df["start_pos"] <= mod2_position)
            & (frag_df["start_pos"] > mod1_position)
        ]
        frag_df = pd.concat([frag_df_n, frag_df_c])

        # If not enough values
        if frag_df.empty:
            return "Na", "Na", "Na"

        # Statistical analysis
        group1 = frag_df.pivot(
            index=["start_pos", "end_pos"], columns="group", values="intensity"
        ).reset_index()
        group1 = group1.fillna(0)

        # Wilcoxon test
        t_test_res_wilcoxon = stats.wilcoxon(
            x=group1["intensity_1"],
            y=group1["intensity_2"],
            zero_method="zsplit",
            alternative="two-sided",
        )

        # T-test for difference
        difference = group1["intensity_2"] - group1["intensity_1"]
        t_test_res_1samp = stats.ttest_1samp(difference, 0, alternative="less")

        # Plotting
        group_long = pd.melt(
            group1,
            id_vars=["start_pos", "end_pos"],
            value_vars=["intensity_1", "intensity_2"],
            var_name="group",
            value_name="intensity",
        )
        group_long["intensity"] = np.log10(group_long["intensity"] + 1)

        gg = (
            ggplot(group_long, aes(x="group", y="intensity", fill="group"))
            + geom_boxplot()
            + theme_bw()
            + labs(x="Group", y="Log10 Intensity")
            + ggtitle(
                f"Wilcoxon p-value: {t_test_res_wilcoxon[1]:.2e}, 1-sample t-test p-value: {t_test_res_1samp[1]:.2e}, Mean difference: {np.mean(difference):.3f}"
            )
        )

        # difference in a dataframe
        difference_df = pd.DataFrame({"intensity": difference})

        # boxplot
        gg_boxplot = (
            ggplot(difference_df, aes(x=0, y="intensity"))
            + geom_boxplot()
            + theme_bw()
            + labs(x="Difference in Intensity")
            + ggtitle(
                f"Wilcoxon p-value: {t_test_res_wilcoxon[1]:.2e}, 1-sample t-test p-value: {t_test_res_1samp[1]:.2e}, Mean difference: {np.mean(difference):.3f}"
            )
        )

        # Saving or returning the plot
        output_filename = f"{filename}_{ions_type}.pdf" if filename else None
        if output_filename:
            gg.save(output_filename)
            print(f"Plot saved to {output_filename}")
            return gg, gg_boxplot
        else:
            return gg, gg_boxplot

    except Exception as e:
        print(f"An error occurred: {e}")
        return None, None, None, None


def percentage_intensity_distribution(
    FG,
    mod1_position=0,
    mod2_position=0,
    filename=None,
):
    # if FG is a dataframe
    if isinstance(FG, pd.DataFrame):
        frag_df = FG
    else:
        frag_df = FG.get_fragment_table_I0()

    # get sum of intensities in spectrum
    total_intensity = sum(FG.its)

    # get the intensity sum for internal and terminal fragment
    frag_df_internal = frag_df[
        (frag_df["start_pos"] != 1)
        & (frag_df["end_pos"] != len(FG.peptidoforms[0].sequence))
    ]
    frag_df_terminal = frag_df[
        (frag_df["start_pos"] == 1)
        | (frag_df["end_pos"] == len(FG.peptidoforms[0].sequence))
    ]

    # site determining internal and terminal fragment
    frag_df_n = frag_df_internal[
        (frag_df_internal["start_pos"] < mod1_position)
        & (frag_df_internal["end_pos"] >= mod1_position)
        & (frag_df_internal["end_pos"] < mod2_position)
    ]
    frag_df_c = frag_df_internal[
        (frag_df_internal["end_pos"] > mod2_position)
        & (frag_df_internal["start_pos"] <= mod2_position)
        & (frag_df_internal["start_pos"] > mod1_position)
    ]
    frag_df_internal_sdi = pd.concat([frag_df_n, frag_df_c])

    frag_df_n = frag_df_terminal[
        (frag_df_terminal["start_pos"] < mod1_position)
        & (frag_df_terminal["end_pos"] >= mod1_position)
        & (frag_df_terminal["end_pos"] < mod2_position)
    ]
    frag_df_c = frag_df_terminal[
        (frag_df_terminal["end_pos"] > mod2_position)
        & (frag_df_terminal["start_pos"] <= mod2_position)
        & (frag_df_terminal["start_pos"] > mod1_position)
    ]
    frag_df_terminal_sdi = pd.concat([frag_df_n, frag_df_c])

    # get the sum of intensities
    total_internal_intensity_non_sdi = sum(frag_df_internal["intensity"]) - sum(
        frag_df_internal_sdi["intensity"]
    )
    total_terminal_intensity_non_sdi = sum(frag_df_terminal["intensity"]) - sum(
        frag_df_terminal_sdi["intensity"]
    )
    total_internal_intensity_site_determining = sum(frag_df_internal_sdi["intensity"])
    total_terminal_intensity_site_determining = sum(frag_df_terminal_sdi["intensity"])

    # print
    print("total intensity: ", total_intensity)

    # dataframe as five rows
    df = pd.DataFrame(
        {
            "group": [
                "internal_non_sdi",
                "internal_sdi",
                "terminal_non_sdi",
                "terminal_sdi",
            ],
            "intensity": [
                total_internal_intensity_non_sdi,
                total_internal_intensity_site_determining,
                total_terminal_intensity_non_sdi,
                total_terminal_intensity_site_determining,
            ],
        }
    )
    # add percentage column
    df["percentage"] = df["intensity"] / total_intensity * 100

    # add color column (blue drakblue red darkred)
    df["color"] = ["#5fa1e8", "#1561b3", "#eb7575", "#a82222"]

    # plot stacked bar plot
    gg = (
        ggplot(df, aes(x=0, y="percentage", fill="group"))
        + geom_bar(position="stack", stat="identity")
        + theme_bw()
        + labs(x="Group", y="Percentage Intensity")
        + scale_fill_manual(values=df["color"])
        + ylim(0, 100)
        + geom_text(aes(label="percentage"), position=position_stack(vjust=0.5), size=5)
    )

    if filename is not None:
        ggsave(gg, filename=filename)

    return gg


# plot the experimental spectrum stored in the FG object
def plot_spectrum(FG, category=["unassigned", "internal", "terminal", "both"]):
    df = FG.get_peak_table()

    # Add new column with color
    df["category"] = "unassigned"
    df.loc[df["internal"] > 0, "category"] = "internal"
    df.loc[df["terminal"] > 0, "category"] = "terminal"
    df.loc[(df["internal"] > 0) & (df["terminal"] > 0), "category"] = "both"

    # Filter dataframe based on category parameter
    df = df[df["category"].isin(category)]

    # Plot the spectrum
    plt.figure(figsize=(14, 5))
    markerline, stemlines, baseline = plt.stem(
        df["mz"],
        df["intensity"],
        linefmt="black",
        markerfmt=" ",
        basefmt=" ",
    )
    plt.setp(stemlines, "linewidth", 0.4)
    plt.xlabel("m/z")
    plt.ylabel("Intensity")

    return plt

def plot_intensity_distribution(FG, filename=None):
    # return a stacked bar plot of the percentage of intensity for internal and terminal fragments

    frag_df = FG.get_fragment_table_I0()

    # get sum of intensities in spectrum
    total_intensity = sum(FG.its)

    # get the intensity sum for internal and terminal fragments
    frag_df_internal = frag_df[
        (frag_df["start_pos"] != 1)
        & (frag_df["end_pos"] != len(FG.peptidoforms[0].sequence))
    ]
    frag_df_terminal = frag_df[
        (frag_df["start_pos"] == 1)
        | (frag_df["end_pos"] == len(FG.peptidoforms[0].sequence))
    ]

    # get the sum of intensities
    total_internal_intensity = sum(frag_df_internal["intensity"])
    total_terminal_intensity = sum(frag_df_terminal["intensity"])

    # calculate unassigned intensity
    total_unassigned_intensity = (
        total_intensity - total_internal_intensity - total_terminal_intensity
    )

    # print
    print("total intensity: ", total_intensity)
    print("total internal intensity: ", total_internal_intensity)
    print("total terminal intensity: ", total_terminal_intensity)
    print("total unassigned intensity: ", total_unassigned_intensity)

    # dataframe as three rows
    df = pd.DataFrame(
        {
            "group": ["internal", "terminal", "unassigned"],
            "intensity": [
                total_internal_intensity,
                total_terminal_intensity,
                total_unassigned_intensity,
            ],
        }
    )
    # add percentage column
    df["percentage"] = df["intensity"] / total_intensity * 100
    df["percentage"] = df["percentage"].apply(
        lambda x: round(x, 2)
    )  # round up the percentage number

    # add color column
    df["color"] = [
        "#154360",
        "#c83636",
        "#ffffff",
    ]  # white color for unassigned intensities

    # plot stacked bar plot
    gg = (
        ggplot(df, aes(x=0, y="percentage", fill="group"))
        + geom_bar(position="stack", stat="identity", color="black", size=0.5)
        + theme_bw()
        + labs(x="Group", y="Percentage Intensity")
        + scale_fill_manual(values=df["color"])
        + ylim(0, 100)
        + geom_text(
            aes(label="percentage"), position=position_stack(vjust=0.5), size=10
        )  # increase font size for percentage values
        + theme(
            panel_grid_major_x=element_blank(), panel_grid_minor_x=element_blank()
        )  # remove vertical grid lines
        + theme(
            panel_grid_major_y=element_line(color="gray", linetype="dashed", size=0.5)
        )  # add dashed horizontal grid lines
        + theme(axis_text_x=element_blank())  # remove x-axis labels
        + theme(axis_ticks_major_x=element_blank())  # remove x-axis ticks
        + theme(
            axis_text_y=element_text(size=12)
        )  # increase font size for y-axis labels
    )

    # save or display the plot
    if filename is not None:  # as svg
        gg.save(filename, width=5, height=10)

    return gg


def plot_fragment_count_distribution(FG, filename=None):
    # return a stacked bar plot of the percentage of fragment count for internal, terminal, and unassigned

    frag_df = FG.get_peak_table()
    print(frag_df)

    # get total number of peaks
    total_peaks = len(frag_df)

    # get the number of peaks assigned to internal and terminal fragments
    frag_df_internal = frag_df[(frag_df["internal"] > 0) & (frag_df["terminal"] == 0)]
    frag_df_terminal = frag_df[(frag_df["terminal"] > 0) & (frag_df["internal"] == 0)]
    frag_df_both = frag_df[(frag_df["internal"] > 0) & (frag_df["terminal"] > 0)]
    # get the number of peaks assigned to internal, terminal, and both fragments
    total_internal_peaks = len(frag_df_internal)
    total_terminal_peaks = len(frag_df_terminal)
    total_both_peaks = len(frag_df_both)

    # calculate unassigned peaks
    total_unassigned_peaks = (
        total_peaks - total_internal_peaks - total_terminal_peaks - total_both_peaks
    )

    # print
    print("total peaks: ", total_peaks)
    print("total internal peaks: ", total_internal_peaks)
    print("total terminal peaks: ", total_terminal_peaks)
    print("total both peaks: ", total_both_peaks)
    print("total unassigned peaks: ", total_unassigned_peaks)

    # dataframe as four rows
    df = pd.DataFrame(
        {
            "group": ["internal", "terminal", "both", "unassigned"],
            "peaks": [
                total_internal_peaks,
                total_terminal_peaks,
                total_both_peaks,
                total_unassigned_peaks,
            ],
        }
    )

    # add percentage column
    df["percentage"] = df["peaks"] / total_peaks * 100
    df["percentage"] = df["percentage"].apply(
        lambda x: round(x, 2)
    )  # round up the percentage number

    # add color column
    df["color"] = [
        "#004D40",
        "#c83636",
        "#ffcc00",
        "#ffffff",
    ]  # yellow color for both fragments, white color for unassigned intensities

    print(df)

    # plot stacked bar plot
    gg = (
        ggplot(df, aes(x=0, y="percentage", fill="group"))
        + geom_bar(position="stack", stat="identity", color="black", size=0.5)
        + theme_bw()
        + labs(x="Group", y="Percentage Peaks")
        + scale_fill_manual(values=df["color"])
        + ylim(0, 100)
        + geom_text(
            aes(label="percentage"), position=position_stack(vjust=0.5), size=10
        )  # increase font size for percentage values
        + theme(
            panel_grid_major_x=element_blank(), panel_grid_minor_x=element_blank()
        )  # remove vertical grid lines
        + theme(
            panel_grid_major_y=element_line(color="gray", linetype="dashed", size=0.5)
        )  # add dashed horizontal grid lines
        + theme(axis_text_x=element_blank())  # remove x-axis labels
        + theme(axis_ticks_major_x=element_blank())  # remove x-axis ticks
        + theme(
            axis_text_y=element_text(size=12)
        )  # increase font size for y-axis labels
    )

    # save or display the plot
    if filename is not None:
        gg.save(filename, width=5, height=10)

    return gg


def plot_mass_spectrum(FG, color_map=None, filename=None):
    """
    Plot a mass spectrum from a given fraggraph.

    :param data: A fraggraph object
    :param color_map: Optional dictionary mapping category1 values to colors
    :param show_arrows: Whether to show arrows for peaks in category2
    """

    data = FG.get_peak_table()

    # Add category1 column
    data["category1"] = "unassigned"
    data.loc[data["internal"] > 0, "category1"] = "internal"
    data.loc[data["terminal"] > 0, "category1"] = "terminal"
    data.loc[(data["internal"] > 0) & (data["terminal"] > 0), "category1"] = "both"

    # Add category2 column
    data["category2"] = True

    # Create a figure and axis for the plot
    fig, ax = plt.subplots(figsize=(16, 6))  # Increase the width of the plotting window

    if color_map is None:
        color_map = {
            "internal": "#FFC107",
            "terminal": "#D81B60",
            "both": "#004D40",
            "unassigned": "#a9aaab",
        }

    # Iterate through the DataFrame and plot each peak
    for index, row in data.iterrows():
        color = color_map[row["category1"]]
        ax.plot(
            [row["mz"], row["mz"]], [0, row["intensity"]], color=color, linewidth=0.5
        )  # Reduce the bar width

        # Optionally add an arrow for category2 peaks
        # if row["category2"]:
        #     ax.annotate(
        #         "",
        #         xy=(row["mz"], row["intensity"]),
        #         xytext=(row["mz"], row["intensity"] + 2),  # Adjust text position
        #         arrowprops=dict(
        #             facecolor="green",
        #             shrink=0.05,
        #             headwidth=1,  # Adjust the width of the arrow head
        #             headlength=1,  # Adjust the length of the arrow head
        #             width=0.5,  # Adjust the width of the arrow shaft
        #         ),
        #     )

    # Set labels and title
    ax.set_xlabel("m/z")
    ax.set_ylabel("Intensity")
    ax.set_title("Mass Spectrum")

    # Show the plot or save the plot
    if filename is not None:  # save as svg
        plt.savefig(filename, format="svg")

    return plt


def plot_mass_spectrum_site_determining(
    FG, exclusion_list=None, color_map=None, filename=None
):
    """
    Plot a mass spectrum from a given fraggraph, highlighting site determining peaks.

    :param FG: A fraggraph object
    :param position1: The first position in the sequence
    :param position2: The second position in the sequence
    :param color_map: Optional dictionary mapping category1 values to colors
    """

    data = FG.get_peak_table()

    # Add category1 column
    data["category1"] = "unassigned"
    data.loc[data["internal"] > 0, "category1"] = "internal"
    data.loc[data["terminal"] > 0, "category1"] = "terminal"
    data.loc[(data["internal"] > 0) & (data["terminal"] > 0), "category1"] = "both"

    # Exclude peaks based on category 1
    if exclusion_list is not None:
        data = data[~data["category1"].isin(exclusion_list)]

    # define the category2 from the proteoform_indicescolulmn
    # unpack all indices, if it contains only the same index, then add this index to the category2 column, if contains different indices, specify "multiple" in category2 column
    category2 = []
    for index, row in data.iterrows():
        # unpack list of list in proteoform_indices column
        proteoform_indices = row["proteoform_indices"]
        proteoform_indices = [
            item for sublist in proteoform_indices for item in sublist
        ]
        # check if it contains only the same index
        if len(set(proteoform_indices)) == 1:
            category2.append(proteoform_indices[0])
        else:
            category2.append("multiple")

    data["category2"] = category2

    # Create a figure and axis for the plot
    fig, ax = plt.subplots(figsize=(16, 6))  # Increase the width of the plotting window

    print("category2: ", category2)

    if color_map is None:
        color_map = {
            "internal": "#FFC107",
            "terminal": "#D81B60",
            "both": "#004D40",
            "unassigned": "#a9aaab",
        }
    # Iterate through the DataFrame and plot each peak
    for index, row in data.iterrows():
        color = color_map[row["category1"]]
        ax.plot(
            [row["mz"], row["mz"]], [0, row["intensity"]], color=color, linewidth=0.5
        )  # Reduce the bar width
        # add symbols according to category2

        if row["category2"] == "multiple":
            pass
        elif row["category2"] != "multiple":
            if int(row["category2"]) == 0:
                marker = "o"
            elif int(row["category2"]) == 1:
                marker = "v"
            elif int(row["category2"]) == 2:
                marker = "s"

            ax.scatter(
                row["mz"],
                row["intensity"],
                marker=marker,
                s=15,
                color="darkgrey",
                linewidths=0.5,
                edgecolors="black",
            )

    # Set labels and title
    ax.set_xlabel("m/z")
    ax.set_ylabel("Intensity")
    ax.set_title("Mass Spectrum (Site Determining)")

    # Show the plot or save the plot
    if filename is not None:  # save as svg
        plt.savefig(filename, format="svg")

    return plt


import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker


def plot_fragment_coverage_sequence(FG, filename=None):

    plots_st = list()

    # get the fragment table
    frag_df = FG.get_fragment_table_I0()

    # get the sequence
    sequence = FG.peptidoforms[0].sequence

    # represent the sequence as a list of letters
    sequence_list = list(sequence)

    # calculate the total spectrum intensity
    total_intensity = frag_df["intensity"].sum()

    # represent percentage of total spectrum intensity per position as a list of zeros
    intensity_list_terminal = [0] * len(sequence)
    intensity_list_internal = [0] * len(sequence)

    # for each terminal fragment, add the percentage of total spectrum intensity to the position covered by that fragment
    for index, row in frag_df.iterrows():
        if row["start_pos"] == 1 or row["end_pos"] == len(sequence):
            for i in range(row["start_pos"] - 1, row["end_pos"]):
                intensity_list_terminal[i] += row["intensity"] / total_intensity * 100

    # for each internal fragment
    for index, row in frag_df.iterrows():
        if row["start_pos"] != 1 and row["end_pos"] != len(sequence):
            for i in range(row["start_pos"] - 1, row["end_pos"]):
                intensity_list_internal[i] += row["intensity"] / total_intensity * 100

    print("intensity_list_terminal: ", intensity_list_terminal)
    print("intensity_list_internal: ", intensity_list_internal)

    # plot the result as heatmap representing the sequence (use letters in the cell), the sequence should be displayed across several rows
    # plot max 15 amino acids per row, if the sequence is longer, add more rows, if the sequence is shorter, add empty cells
    # plot the intensity as color, use a color map from white to blue, the more intense the fragment, the darker the color

    n = 10

    # make matrix of sequence, all row must be the same length, add 0 to the end of the list if the length is not enough
    sequence_matrix = []
    for i in range(0, len(sequence), n):
        sequence_matrix.append(sequence_list[i : i + n])
    sequence_matrix[-1] += [""] * (n - len(sequence_matrix[-1]))

    # make matrix of intensity, all row must be the same length, add 0 to the end of the list if the length is not enough
    intensity_matrix_terminal = []
    for i in range(0, len(sequence), n):
        intensity_matrix_terminal.append(intensity_list_terminal[i : i + n])
    intensity_matrix_terminal[-1] += [0] * (n - len(intensity_matrix_terminal[-1]))

    intensity_matrix_internal = []
    for i in range(0, len(sequence), n):
        intensity_matrix_internal.append(intensity_list_internal[i : i + n])
    intensity_matrix_internal[-1] += [0] * (n - len(intensity_matrix_internal[-1]))

    print("sequence_matrix: ", sequence_matrix)
    print("intensity_matrix: ", intensity_matrix_terminal)

    sequence_matrix = np.array(sequence_matrix).T
    intensity_matrix_terminal = np.array(intensity_matrix_terminal).T
    intensity_matrix_internal = np.array(intensity_matrix_internal).T

    z = intensity_matrix_terminal
    # ffab03
    nx, ny = z.shape
    indx, indy = np.arange(nx), np.arange(ny)
    x, y = np.meshgrid(indx, indy)

    fig, ax = plt.subplots()

    cmap = matcolors.LinearSegmentedColormap.from_list("", ["white", "#D81B60"])

    im = ax.imshow(z.T, interpolation="nearest", cmap=cmap)

    for xval, yval in zip(x.flatten(), y.flatten()):
        zval = z[xval, yval]
        t = sequence_matrix[xval][yval]
        c = "black"
        ax.text(xval, yval, t, color=c, va="center", ha="center")

    xlabels = "abcdefghij"
    ylabels = "0123456789"
    ax.set_xticks(indx + 0.5)
    ax.set_yticks(indy + 0.5)
    ax.grid(ls="-", lw=2)

    for a, ind, labels in zip((ax.xaxis, ax.yaxis), (indx, indy), (xlabels, ylabels)):
        a.set_major_formatter(ticker.NullFormatter())
        a.set_minor_locator(ticker.FixedLocator(ind))
        a.set_minor_formatter(ticker.FixedFormatter(labels))

    ax.xaxis.tick_top()

    # Add colorbar
    cbar = fig.colorbar(im)
    cbar.set_label("Percentage of Total Spectrum Intensity")

    # output if filemane is specified
    if filename is not None:
        plt.savefig(str(filename) + "_terminal.svg", format="svg")

    plots_st.append(fig)

    # plot internal
    z = intensity_matrix_internal

    nx, ny = z.shape
    indx, indy = np.arange(nx), np.arange(ny)
    x, y = np.meshgrid(indx, indy)

    fig, ax = plt.subplots()
    cmap = matcolors.LinearSegmentedColormap.from_list("", ["white", "#FFC107"])
    im = ax.imshow(z.T, interpolation="nearest", cmap=cmap)

    for xval, yval in zip(x.flatten(), y.flatten()):
        zval = z[xval, yval]
        t = sequence_matrix[xval][yval]
        c = "black"
        ax.text(xval, yval, t, color=c, va="center", ha="center")

    xlabels = "abcdefghij"
    ylabels = "0123456789"
    ax.set_xticks(indx + 0.5)
    ax.set_yticks(indy + 0.5)
    ax.grid(ls="-", lw=2)

    for a, ind, labels in zip((ax.xaxis, ax.yaxis), (indx, indy), (xlabels, ylabels)):
        a.set_major_formatter(ticker.NullFormatter())
        a.set_minor_locator(ticker.FixedLocator(ind))
        a.set_minor_formatter(ticker.FixedFormatter(labels))

    ax.xaxis.tick_top()

    # Add colorbar
    cbar = fig.colorbar(im)
    cbar.set_label("Percentage of Total Spectrum Intensity")

    # output if filemane is specified
    if filename is not None:
        plt.savefig(str(filename) + "_internal.svg", format="svg")

    plots_st.append(fig)

    return plots_st
