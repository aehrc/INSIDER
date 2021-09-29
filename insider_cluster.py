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
import argparse
import sys
from pathlib import Path

# External imports
import pandas as pd
import json

# Internal imports
THIS_DIR = Path(__file__).resolve().parent
sys.path.append(str(THIS_DIR.parent))       ## Allow us to import from SRC
from src.kmer import analyse
from src.kmer import transform
from src.util import spark
from src import io

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def consensusKmerFrequencies(obsDirs, oFile, pFile):
    def _setup(kmerDf):
        ## **********
        ## *** Sequences that are relatively short are likely to be significantly
        ## *** different to other sequences because their counts are proportionally
        ## *** more important. Conversely, sequences that are relatively long are
        ## *** likely to be similar to other sequences because their counts
        ## *** ultimately dictate what the estimated/true signature will look like.
        ## ***
        ## *** Because of this, there could be biases due to sequence length
        ## *** differences. So to improve confidence, we may have to include
        ## *** a size selection step, whereby we group sequences by length
        ## *** and analyze each group individually. But for now, we'll just
        ## *** remove the relatively short sequences
        ## **********
        kmerDf = transform.filter.removeShortSequences(kmerDf)

        ## **********
        ## *** We assume that the Kmer frequencies were calculated in only one
        ## *** orientation. However, we don't actually know the correct orientation of
        ## *** each sequence. So to account for this, we will calculate the counts
        ## *** in both directions, sum up the counts of complementary Kmers and
        ## *** collapse them into a single feature.
        ## **********
        kmerDf = transform.agg.clusterRevCompCounts(kmerDf)
        kmerDf.show()
        return kmerDf

    print("Consensus clustering")
    params = get_spark_params()
    with spark.getSparkSession(params) as ss:
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
            kmerDf    = _setup(kmerDf)
            params    = readParams(pFile)
            clusterDf = analyse.cluster.getLabels(kmerDf, params)

            print("Cleaning up output")
            clusterDf = clusterDf.repartition(sc.defaultParallelism)
            clusterDf = clusterDf.toPandas()

            ## Write cluster labels to file
            print("Writing output")
            clusterDf.to_csv(oFile, sep='\t', index=False)

#------------------- Private Classes & Functions ------------#

def readParams(pFile):
    ## Default parameters
    drKwargs   = {}
    dcKwargs   = {'eps':0.05}
    dconKwargs = {'eps':0.5}
    params = {'rkwargs':drKwargs, 'ckwargs':dcKwargs, 'conkwargs':dconKwargs}

    if (pFile is not None):
        with open(pFile) as f:
            data = json.load(f)
            params['rkwargs']   = params['rkwargs'] if 'rkwargs' not in data else data['rkwargs']
            params['ckwargs']   = params['ckwargs'] if 'ckwargs' not in data else data['ckwargs']
            params['conkwargs'] = params['conkwargs'] if 'conkwargs' not in data else data['conkwargs']
    return params

def get_spark_params():
    params = [
        ## Driver
        ('spark.driver.memory', '128G'),
        ('spark.driver.maxResultSize', '0'),

        ## Executor
        # ('spark.executor.cores', '5'),         ## Didn't seem to have any effect...
        # ('spark.executor.instances', '21'),    ## Didn't seem to have any effect...
        ('spark.executor.memory', '64G'),
        ('spark.executor.heartbeatInterval', '60s'),

        ## SQL
        ('spark.sql.broadcastTimeout', '600s'),
        ('spark.sql.execution.arrow.pyspark.enabled', 'true'),
        ('spark.sql.shuffle.partitions', '200'),

        ## Misc
        ('spark.local.dir', THIS_DIR),
        ('spark.network.timeout', '600s'),
        ('spark.default.parallelism', '8')
    ]
    return params

def make_parser():
    def initArgs(p):
        p.add_argument("--freqDir", help="Observed Kmer frequencies directory",
            nargs='*', type=str, required=True)
        p.add_argument("-o", help="Output file (.tsv)",
            type=str, required=True)
        p.add_argument("--params", help="Clustering algorithm \
            parameters (.json)", type=str)

    parser    = argparse.ArgumentParser(description='Compute \
        oligonucleotide frequencies of sequence/s')
    subparser = parser.add_subparsers(dest='command')
    distance  = subparser.add_parser('distance')
    distance.add_argument("--exp", help="Expected Kmer frequencies",
        nargs='*', type=str, required=True)
    consensus = subparser.add_parser('consensus')
    initArgs(distance)
    initArgs(consensus)
    return parser

def main(parser):
    args = parser.parse_args()
    if (args.command is None):
        parser.print_help()
        sys.exit(1)

    else:
        ## There are two ways we can do clustering:
        ## 1. Compute a distance matrix between sequences and
        ##    perform clustering on the distance matrix
        ## 2. Randomly cluster sequences and identify sequences
        ##    that consistently group together
        if (args.command == 'distance'):
            print(args)
            # distanceKmerFrequencies(args.obs, args.exp,
            #     args.o, args.params)
            pass

        elif (args.command == 'consensus'):
            print(args)
            consensusKmerFrequencies(args.freqDir, args.o, args.params)

    print("DONE")

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    parser = make_parser()
    main(parser)

#------------------------------------------------------------------------------
