#!/usr/bin/env python3

import pandas as pd
import plotly.graph_objects as go
import numpy as np

def plot_spectrum(intensity_values, mz_values) -> go.Figure:
    fig = go.Figure()

    for mz, i in zip(mz_values, intensity_values):
        fig.add_trace(go.Scatter(x=[mz, mz], y=[0, i], mode='lines', line=dict(color='black')))

    fig.update_layout(xaxis_title='m/z', yaxis_title='intensity', yaxis=dict(range=[0, max(intensity_values)]))

    for trace in fig.data:
        trace.showlegend = False
    
    return fig


import plotly.graph_objects as go

def plot_consensus_spectrum(intensity_values, mz_values, coverage) -> go.Figure:
    fig = go.Figure()
    # Define a color gradient from black (cov = 1) to yellow (cov = 0).
    get_color = lambda v: (int(255 * (1 - v)), int(255 * (1 - v)), 0) if v < 1 else (0, 0, 0)

    for mz, i, cov in zip(mz_values, intensity_values, coverage):
        fig.add_trace(go.Scatter(x=[mz, mz], y=[0, i], mode='lines', line=dict(color='rgb' + str(get_color(cov)))))

    fig.update_layout(xaxis_title='m/z', yaxis_title='intensity', yaxis=dict(range=[0, max(intensity_values)]))

    for trace in fig.data:
        trace.showlegend = False

    return fig




def plot_spectra_chromatogram(spectra: dict) -> pd.DataFrame:
    """
    Returns a DataFrame with the following columns:
    "scan_nr", "precursor", "charge", "rt", "max_intensity", "peaks"
    """
    
    fig = go.Figure()
    
    #get the min and max intensity
    min_intensity = min([spectrum["precursor"][1] for spectrum in spectra.values()])
    max_intensity = max([spectrum["precursor"][1] for spectrum in spectra.values()])
    print("min/max intensity: ", min_intensity, max_intensity)
    colorscale = [[0, 'blue'], [1, 'red']]  # Define your desired colorscale
    
    # check whether min and max intensity are the same
    if min_intensity == max_intensity:
        min_intensity = 0
        max_intensity = 1
        
    
    for scan_nr, spectrum in spectra.items():
        intensity = spectrum["precursor"][1]
        color = (intensity - min_intensity) / (max_intensity - min_intensity)  # Calculate the color value based on intensity
        
        fig.add_trace(go.Scatter(
            x=[spectrum["rt"]],
            y=[spectrum["precursor"][0]],
            mode='markers',
            marker=dict(size=10, color=color, colorscale=colorscale),
            text=scan_nr,
            name=scan_nr
        ))
    
    fig.update_layout(
        title='Spectra in retention time / precursor mass',
        xaxis_title='retention time',
        yaxis_title='precursor mass',
        showlegend=False
    )
    
    return fig


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


#fragment coverage matrix plotly

def draw_fragment_coverage_matrix_plotly(
    FG, x="intensity", x_min=1, x_max=10, filename=None, peptidoform_index=0
):
    """
    This function draws a fragment coverage matrix using plotly.

    Parameters:
    FG: Fragment Graph object
    x (str): The column name to use for the intensity. Default is "intensity".
    x_min (int): The minimum limit for the color scale. Default is 1.
    x_max (int): The maximum limit for the color scale. Default is 10.
    filename (str): The filename to save the plot. If None, the plot is returned instead of being saved. Default is None.
    peptidoform_index (int): The index of the peptidoform to use. Default is 0.

    Returns:
    plotly.graph_objects.Figure object if filename is None, else None.
    """
    
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

    # Create plotly figure
    fig = go.Figure(
        data=go.Heatmap(
            z=frag_df[x],
            x=frag_df["end_pos"],
            y=frag_df["start_pos"],
            colorscale=colors,
            zmin=x_min,
            zmax=x_max,
            colorbar=dict(title=x),
        )
    )

    # Update layout
    fig.update_layout(
        xaxis_title = "End Position (AA index)",
        yaxis_title = "Start Position (AA index)",
        height = 500,
        width = 500,
        title = "Fragment Coverage Matrix",
    )

    #return plot:
    fig['layout']['xaxis']['autorange'] = "reversed"
    
    return fig


#fragment coverage matrix difference plotly

def draw_fragment_coverage_matrix_difference_plotly(
    FG,
    peptidoform_index_1=0,
    peptidoform_index_2=1,
):
    """
    This function draws a fragment coverage matrix difference using plotly.

    Parameters:
    FG1: Fragment Graph object 1
    FG2: Fragment Graph object 2
    x (str): The column name to use for the intensity. Default is "intensity".
    FG1_peptidoform_index (int): The index of the peptidoform to use for FG1. Default is 0.
    FG2_peptidoform_index (int): The index of the peptidoform to use for FG2. Default is 0.

    Returns:
    plotly.graph_objects.Figure object if filename is None, else None.
    """
    
    # Extracting fragment data
    frag_df_1 = FG.get_fragment_table_I0(peptidoform_index=peptidoform_index_1)
    frag_df_2 = FG.get_fragment_table_I0(peptidoform_index=peptidoform_index_2)

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
    
    # # Annotations
    # frag_df["annotation"] = frag_df["FC"].apply(
    #     lambda x: "✖" if x == -10 else ("●" if x == 10 else "")
    # )
    
    # Create plotly figure
    
    #print head of frag_df
    

    
    fig = go.Figure(
        data=go.Heatmap(
            z=frag_df["FC"],
            x=frag_df["end_pos"],
            y=frag_df["start_pos"],
            colorscale="RdBu",
            zmin=-10,
            zmax=10,
            colorbar=dict(title="Log2FC Intensity"),
        )
    )
    
    # Update layout
    
    fig.update_layout(
        xaxis_title="End Position (AA index)",
        yaxis_title="Start Position (AA index)",
        height=600,
        width=600,
        title="Fragment Coverage Matrix Difference",
    )
    
    fig['layout']['xaxis']['autorange'] = "reversed"
    
    
    
    return fig


# define the function draw_barplor_intensity_SDI, this function take a fragment graph object and a list of position ranges,
#if the position range is not provided, the function will plot the intensity of site determining ions automatically 
# for the peptidoform of peptidoform_index.
#the site determine ions are ions of the sae position but whose mass are different across the peptidoforms

def draw_barplot_intensity_SDI(FG, position_range = None):
    """
    This function draws a barplot of the intensity of site determining ions (SDI) for a given Fragment Graph object.

    Parameters:
    FG: Fragment Graph object
    position_ranges (list): A list of tuples with the start and end positions of the SDI. If None, the SDI are automatically determined.

    Returns:
    plotly.graph_objects.Figure object
    """
    
    
    
    
    # Get fragment table I0 from FG (assuming it returns a DataFrame)
    frag_df = FG.get_all_fragment_table_I0()
    
    if position_range is None:
        frag_df = frag_df.drop_duplicates(subset=["start_pos", "end_pos", "mz"], keep=False)
    else: 
        #use the position ranges to filter the fragment table
        #TODO check that this is correct
        frag_df = frag_df[(frag_df["start_pos"] >= position_range[0]) & (frag_df["end_pos"] <= position_range[1])]
        
        
    
    #boxplot intensiyt in fucntion of peptidoform index
    
    fig = go.Figure(
        data=[
            go.Box(
                y=frag_df["intensity"],
                x=frag_df["peptidoform_index"],
                name="Intensity",
                marker_color="blue",
            )
        ]
    )
    
    # Update layout
    
    fig.update_layout(
        xaxis_title="Peptidoform Index",
        yaxis_title="Intensity",
        title="Barplot of Intensity SDI",
    )
    
    return fig


# create a jitter plot of the intensity of the site determining ions for each peptidoform and connect the points with a line if it is the same position

def draw_jitterplot_intensity_SDI(FG, position_range = None):
    
    """
    This function draws a jitter plot of the intensity of site determining ions (SDI) for a given Fragment Graph object.

    Parameters:
    FG: Fragment Graph object
    position_ranges (list): A list of tuples with the start and end positions of the SDI. If None, the SDI are automatically determined.

    Returns:
    plotly.graph_objects.Figure object
    """
    
    # Get fragment table I0 from FG (assuming it returns a DataFrame)
    frag_df = FG.get_all_fragment_table_I0()
    
    if position_range is None:
        frag_df = frag_df.drop_duplicates(subset=["start_pos", "end_pos", "mz"], keep=False)
    else: 
        #use the position ranges to filter the fragment table
        #TODO check that this is correct
        frag_df = frag_df[(frag_df["start_pos"] >= position_range[0]) & (frag_df["end_pos"] <= position_range[1])]
    
    #create a jitter plot of the intensity of the site determining ions for each peptidoform and connect the points with a line if it is the same position
    
    fig = go.Figure()
    
    for i in range(len(frag_df)):
        fig.add_trace(go.Scatter
                      
                        (x = frag_df["peptidoform_index"],
                         y = frag_df["intensity"],
                         mode = "lines+markers",
                         name = "Intensity",
                         marker_color = "blue"))
        
    # Update layout
    
    
    fig.update_layout(
        xaxis_title="Peptidoform Index",
        yaxis_title="Intensity",
        title="Jitterplot of Intensity SDI",
    )
    
    return fig
    
    
    
    

    
    
    
    
        

        
        

    