
import networkx as nx
import matplotlib.axes
import matplotlib.pyplot as plt

def draw_DAG(G: nx.DiGraph, 
             ax: matplotlib.axes.Axes | None=None
             ) -> None:
    """
    Input:
        G: A networkx DiGraph where internal nodes are tuples of leaves
        ax: 
    Output:
        None
    
    Draws the network with layers aligned using matplotlib and networkx multipartite layout.
    """

    for layer, nodes in enumerate(nx.topological_generations(G)):
        for node in nodes:
            G.nodes[node]["layer"] = layer

    # TODO: This is bugged and sometimes draws the network completley wrong. THe idea is to put all leaves on the final layer.
    # leaf_list = []
    # for node in G.nodes:
    #     if type(node) != tuple and node != "rho":   # When constructing the canonical dag/network for R, we relabeled all leaves to not be tuples anymore.
    #         leaf_list.append(node)
    # leaf_layer = max(G.nodes[node]["layer"] for node in leaf_list)

    # for node in leaf_list:
    #     G.nodes[node]["layer"] = leaf_layer

    pos = nx.multipartite_layout(G, subset_key="layer", align="horizontal")

    nx.draw(G, pos=pos, ax=ax)
    nx.draw_networkx_labels(G, pos=pos, ax=ax)
    plt.show()


def plot_results(realizable: bool, results_list: list) -> None:
    if realizable:
        plot_realizable_results(results_list)
    else:
        plot_non_realizable_results(results_list)


def plot_realizable_results(results_list) -> None:
    pass


def plot_non_realizable_results(results_list) -> None:
    pass
