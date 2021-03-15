#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import json
import os
from pathlib import Path

# External imports
import networkx as nx
import pandas as pd

# Internal imports
from .common import createDirIfNone

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def write(filepath, kmerGraph, prefix='', formatAs='tsv'):
    if (formatAs == 'tsv'):
        _writeTsv(filepath, kmerGraph, prefix)

    elif (formatAs == 'json'):
        _writeJson(filepath, kmerGraph)

#------------------- Private Classes & Functions ------------#

def _writeTsv(filepath, kmerGraph, prefix):
    createDirIfNone(filepath)

    ## Write edges
    edgeFilename = '{}_edges.txt'.format(prefix)
    edgeFilepath = Path(filepath, edgeFilename)
    edgeDf = nx.to_pandas_edgelist(kmerGraph)
    edgeDf.to_csv(edgeFilepath, sep='\t', index=False)

    ## Write nodes
    nodeFilename = '{}_nodes.txt'.format(prefix)
    nodeFilepath = Path(filepath,nodeFilename)
    nodeDf = dict(kmerGraph.nodes(data=True))
    nodeDf = pd.DataFrame.from_dict(nodeDf, orient='index').reset_index()
    nodeDf = nodeDf.rename(columns={'index':'id'})
    nodeDf.to_csv(nodeFilepath, sep='\t', index=False)

def _writeJson(filepath, kmerGraph):
    with open(filepath, 'w') as f:
        s = nx.json_graph.cytoscape_data(kmerGraph)
        json.dumps(s, filepath, indent=4)

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()
