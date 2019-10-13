import os
from collections import ChainMap
import glob

import scipy.stats as stats
import pandas as pd

import hetio.readwrite
import hetmech.matrix

def load_hetnets(hetnet_file, permuted_directory, subset_genes, metaedge_abbrev='GpXCELL'):
    """
    Load in real and permuted hetnets and store in a dictionary. A hetnet is a
    "heterogeneous network" described in https://neo4j.het.io/browser/
    Arguments:
    hetnet_file - the file path of the real data hetnet
    permuted_directory - the directory where permuted hetnets are stored
    subset_genes - the gene identifiers to use in the adjacency matrices
    metaedge_abbrev - the abbreviation to use for loading metaedge graph
                      (default: 'GpXCELL' - cell type genesets)
    Output:
    A dictionary of real and permuted hetnets
    """

    paths = os.listdir(permuted_directory)
    idx = 0
    hetnet_dict = {}
    for path in paths:
        if path != 'stats.tsv':

            graph = hetio.readwrite.read_graph(
                os.path.join(permuted_directory, path)
            )

            graph = hetmech.matrix.metaedge_to_adjacency_matrix(
                graph, metaedge=metaedge_abbrev
                )

            graph = pd.DataFrame(graph[2], index=graph[0], columns=graph[1])
            graph.index = graph.index.map(str)
            graph = graph.reindex(subset_genes).fillna(0) * 1
            hetnet_dict[idx] = graph
            idx += 1

    graph = hetio.readwrite.read_graph(hetnet_file)
    graph = hetmech.matrix.metaedge_to_adjacency_matrix(
        graph, metaedge=metaedge_abbrev
        )
    graph = pd.DataFrame(graph[2], index=graph[0], columns=graph[1])
    graph.index = graph.index.map(str)
    graph = graph.reindex(subset_genes).fillna(0) * 1
    hetnet_dict['real'] = graph

    return hetnet_dict