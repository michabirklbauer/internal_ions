#!/usr/bin/env python3

import pandas as pd
import plotly.graph_objects as go

def example_plot(dataframe: pd.DataFrame) -> go.Figure:

    return go.Figure([go.Histogram(x = dataframe["values"])])
