#!/usr/bin/env python3

import pandas as pd
import plotly.graph_objects as go

def plot_spectrum(intensity_values, mz_values) -> go.Figure:
    fig = go.Figure()

    for mz, i in zip(mz_values, intensity_values):
        fig.add_trace(go.Scatter(x=[mz, mz], y=[0, i], mode='lines', line=dict(color='black')))

    fig.update_layout(xaxis_title='m/z', 
                      yaxis_title='Intensity', 
                      yaxis=dict(range=[0, max(intensity_values)]),
                      xaxis=dict(gridcolor='transparent'), 
        margin=dict(
            t=0,  
            b=0, 
            l=0,  
            r=0,  
            pad=0  
        ))

    for trace in fig.data:
        trace.showlegend = False
    
    return fig
