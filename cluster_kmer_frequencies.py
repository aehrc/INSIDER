#!/bin/python

#------------------- Description & Notes --------------------#

'''
For larger data sets, we seem to run out of memory before running
out of time. So we should be thinking of a different/better approach
if we want to run this on a large-scale. Another thing would be to
distinguish between read, contig, and genome sequences as the underlying
basis of each sequence differs (e.g., non-overlapping Vs overlapping, small
number Vs high number of sequences)
'''

#------------------- Dependencies ---------------------------#

# Standard library imports
import sys

# External imports
import pandas as pd

# Internal imports
from src.kmer import analyse as kmeranalyse
from src.util import params
from src.util import spark
from src import io
from src import kmer

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

#------------------- Private Classes & Functions ------------#

def distanceKmerFrequencies(obsDirs, expDirs, oFile, algo, pFile):
    print("Distance clustering")
    with spark.getSparkSession() as ss:
        with ss.sparkContext as sc:
            ## Instead of reading the Kmer frequencies back into memory
            ## we could extend the pipeline here so that all the computation
            ## and analysis is done in parallel. This would essentially be
            ## a merging of the two Dev scripts (along with a few
            ## other changes)

            ## Read Kmer frequencies
            print("Reading Kmers")
            oKmerDf = io.kmer.read(*obsDirs)
            eKmerDf = io.kmer.read(*expDirs)

            ## Assign each sequence to a cluster
            labels = kmeranalyse.cluster.distance.getLabels(oKmerDf, algo, eKmerDf, **kwargs)

            ## Write cluster labels to file
            labels.to_csv(oFile, sep='\t', index=False)

def consensusKmerFrequencies(obsDirs, oFile, pFile):
    ## Read algorithm params
    params = io.readJson(pFile) if pFile is not None else pd.DataFrame()

    print("Consensus clustering")
    with spark.getSparkSession() as ss:
        with ss.sparkContext as sc:
            ## Instead of reading the Kmer frequencies back into memory
            ## we could extend the pipeline here so that all the computation
            ## and analysis is done in parallel. This would essentially be
            ## a merging of the two Dev scripts (along with a few
            ## other changes)

            ## Read Kmer frequencies
            print("Reading Kmers")
            kmerDf = io.kmer.read(*obsDirs)

            ## Assign each sequence to a cluster
            print("Finding clusters")
            kmerDf    = kmeranalyse.setup(kmerDf)
            clusterDf = kmeranalyse.cluster.consensus.getLabels(kmerDf, params)

            print("Cleaning up output")
            clusterDf    = clusterDf.repartition(sc.defaultParallelism)
            clusterDf    = clusterDf.toPandas()

            ## Write cluster labels to file
            print("Writing output")
            clusterDf.to_csv(oFile, sep='\t', index=False)

def main():
    argParser = params.ClusterKmerArgParser()
    argParser.parse()

    ## There are two ways we can do clustering:
    ## 1. Compute a distance matrix between Kmer frequencies and 
    ##    perform clustering on the distance matrix
    ## 2. Reduce Kmer frequencies into a smaller number of features and
    ##    perform clustering on the reduced features
    if (argParser.cmd == 'distance'):
        argParser.parseDistanceArgs()
        argParser.printArgs()

        distanceKmerFrequencies(argParser.obsDirs, argParser.expDirs,
            argParser.oFile, argParser.algo, argParser.pFile)

    elif (argParser.cmd == 'consensus'):
        argParser.printArgs()

        consensusKmerFrequencies(argParser.obsDirs,
            argParser.oFile, argParser.pFile)

    print("DONE")

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
