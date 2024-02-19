#!/usr/bin/env python3

import pandas as pd
import plotly.graph_objects as go

def plot_spectrum(intensity_values, mz_values) -> go.Figure:
    fig = go.Figure()

    for mz, i in zip(mz_values, intensity_values):
        fig.add_trace(go.Scatter(x=[mz, mz], y=[0, i], mode='lines', line=dict(color='black')))

    fig.update_layout(xaxis_title='m/z', yaxis_title='intensity', yaxis=dict(range=[0, max(intensity_values)]))

    for trace in fig.data:
        trace.showlegend = False
    
    return fig
