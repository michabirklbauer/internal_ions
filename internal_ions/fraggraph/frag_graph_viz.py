# import
from pyvis import network as net
import plotly.graph_objects as go

def draw_graph3(
    networkx_graph,
    output_filename="graph.html",
    show_buttons=True,
    customizable=False
):
    """
    This function accepts a networkx graph object,
    converts it to a pyvis network object preserving its node and edge attributes,
    and both returns and saves a dynamic network visualization.
    Valid node attributes include:
        "size", "value", "title", "x", "y", "label", "color".
        (For more info: https://pyvis.readthedocs.io/en/latest/documentation.html#pyvis.network.Network.add_node)
    Valid edge attributes include:
        "arrowStrikethrough", "hidden", "physics", "title", "value", "size"
        (For more info: https://pyvis.readthedocs.io/en/latest/documentation.html#pyvis.network.Network.add_edge)
    Args:
        networkx_graph: The graph to convert and display
        notebook: Display in Jupyter?
        output_filename: Where to save the converted network
        show_buttons: Show buttons in saved version of network?
        only_physics_buttons: Show only buttons controlling physics of network?
        height: height in px or %, e.g, "750px" or "100%
        size: size in px or %, e.g, "750px" or "100%
        bgcolor: background color, e.g., "black" or "#222222"
        font_color: font color,  e.g., "black" or "#222222"
        pyvis_options: provide pyvis-specific options (https://pyvis.readthedocs.io/en/latest/documentation.html#pyvis.options.Options.set)
    """

    # make a pyvis network
    network_class_parameters = {
        "notebook": False,
        "height": "900px",
        "width": "900px",
        "bgcolor": "white",
        "font_color": "black",
        "cdn_resources": "remote"
    }
    pyvis_graph = net.Network(
        **{
            parameter_name: parameter_value
            for parameter_name, parameter_value in network_class_parameters.items()
            if parameter_value
        },
    )

    # for each node and its attributes in the networkx graph
    for node, node_attrs in networkx_graph.nodes(data=True):
        pyvis_graph.add_node(node, **node_attrs)

    # for each edge and its attributes in the networkx graph
    for source, target, edge_attrs in networkx_graph.edges(data=True):
        # if value/size not specified directly, and weight is specified, set 'value' to 'weight'
        if not "value" in edge_attrs and not "size" in edge_attrs and "weight" in edge_attrs:
            # place at key 'value' the weight of the edge
            edge_attrs["value"] = edge_attrs["weight"]
        # add the edge
        pyvis_graph.add_edge(source, target, **edge_attrs)

    # turn buttons on
    # if show_buttons:
    #     pyvis_graph.show_buttons()

    # pyvis-specific options
    options = """
        var options = {
        "configure": {
                "enabled": false
        },
        "nodes": {
            "font": {
                "size": 3
            },
            "size": 1
        },
        "edges": {
            "color": {
            "inherit": true
            },
            "smooth": false
        },
        "physics": {
            "barnesHut": {
            "gravitationalConstant": -12050
            },
            "minVelocity": 5,
            "timestep": 0.2
        }
        }
        """

    if customizable:
        options = """
            var options = {
            "configure": {
                    "enabled": true
            },
            "nodes": {
                "font": {
                    "size": 3
                },
                "size": 1
            },
            "edges": {
                "color": {
                "inherit": true
                },
                "smooth": false
            },
            "physics": {
                "barnesHut": {
                "gravitationalConstant": -12050
                },
                "minVelocity": 5,
                "timestep": 0.2
            }
            }
            """

    pyvis_graph.set_options(options)

    # return and also save
    return pyvis_graph.write_html(output_filename)


def draw_annotated_spectrum(graph):

    # create a empty plotly figure
    fig = go.Figure()

    # plot the peaks of the spectrum as vertical lines
    # x axis is mz and height of the line is intensity
    for mz, intensity in zip(graph.mzs, graph.its):
        peak_name = str(mz) + "_0"
        # check whether peak node is annotated (more than one connected node)
        if len(list(graph.neighbors(peak_name))) == 1:
            fig.add_trace(
                go.Scatter(x=[mz, mz], y=[0, intensity], mode="lines", line=dict(color="black", width=1))
            )
        elif len(list(graph.neighbors(peak_name))) > 1:
            # get types of connected nodes:
            types = [graph.nodes[neighbor]["frag_dir"] for neighbor in graph.neighbors(peak_name)]

            # if terminal node, color it green
            if ("C" in types or "N" in types) and "internal" not in types:
                fig.add_trace(
                    go.Scatter(
                        x=[mz, mz], y=[0, intensity], mode="lines", line=dict(color="#0a9396", width=1)
                    )
                )
            elif "C" in types or "N" in types and "internal" in types:
                fig.add_trace(
                    go.Scatter(x=[mz, mz], y=[0, intensity], mode="lines", line=dict(color="green", width=1))
                )
            elif "internal" in types:
                fig.add_trace(
                    go.Scatter(
                        x=[mz, mz], y=[0, intensity], mode="lines", line=dict(color="#ee9b00", width=1)
                    )
                )

    fig.update_layout(plot_bgcolor="white")

    return fig
