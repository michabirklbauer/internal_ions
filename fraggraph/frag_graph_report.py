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
        else:
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
    filename="percentage_intensity_distribution.pdf",
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
    plt.show()

    return plt
