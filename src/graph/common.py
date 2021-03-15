#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import re

# External imports
import numpy as np
import networkx as nx

# Internal imports
from .. import distance

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def constructFromDistance(kmerId, kmerCount, **kwargs):
    ## Drop columns that are the same in every sequence and 
    ## Evaluate distances between Kmer frequency distributions
    kmerDist  = distance.calculate(kmerId, kmerCount, **kwargs)

    ## Construct a directed graph
    kmerGraph = nx.from_pandas_edgelist(kmerDist, 'id_x', 'id_y',
        edge_attr='distance', create_using=nx.DiGraph)

    ## Add node attributes
    nodeAttrs = kmerId.set_index('id').to_dict('index')
    nx.set_node_attributes(kmerGraph, nodeAttrs)
    return kmerGraph

def filterSelfEdges(kmerGraph):
    newKmerGraph = kmerGraph.copy()
    edgeIds      = list(nx.selfloop_edges(newKmerGraph))
    newKmerGraph.remove_edges_from(edgeIds)
    return newKmerGraph

def filterByDistance(kmerGraph, d):
    newKmerGraph = kmerGraph.copy()
    f            = lambda x: x[2]['distance'] < d
    edges        = filter(f, newKmerGraph.edges(data=True))
    edgeIds      = [e[:2] for e in edges]
    newKmerGraph = newKmerGraph.edge_subgraph(edgeIds)
    return newKmerGraph

def filterOutliers(kmerGraph):
    dists = [a['distance'] for s, t, a in kmerGraph.edges(data=True)]
    dists = np.array(dists)

    ## Calculate the quartiles of a boxplot for distances
    uq  = np.percentile(dists, 75)
    lq  = np.percentile(dists, 25)
    iqr = uq - lq

    ## Calculate the end whiskers of a boxplot for distances
    uw  = float(dists[dists <= uq + 1.5 * iqr].max())
    lw  = float(dists[dists >= lq - 1.5 * iqr].min())

    newKmerGraph = kmerGraph.copy()
    f            = lambda x: ((x[2]['distance'] > lw) & (x[2]['distance'] < uw))
    edges        = filter(f, newKmerGraph.edges(data=True))
    edgeIds      = [e[:2] for e in edges]
    newKmerGraph = newKmerGraph.edge_subgraph(edgeIds)
    return newKmerGraph

def findById(kmerGraph, sIds, tIds=None):
    p = '|'.join(sIds)
    f = lambda x: re.search(p, x[0])

    nodes   = list(filter(f, kmerGraph.nodes(data=True)))
    nodeIds = [n[0] for n in nodes]
    df      = nx.to_pandas_edgelist(kmerGraph, nodelist=nodeIds)
    if (tIds is not None):
        p    = '|'.join(tIds)
        cond = df['target'].str.contains(p)
        return df[cond]
    return df

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
