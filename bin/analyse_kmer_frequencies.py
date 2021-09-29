#!/bin/python

#------------------- Description & Notes --------------------#

'''
Description:
    Given a directory containing a list of files, output (to file)
    the results of analyses. Each file in the directory is compressed
    in Parquet format and contains the frequencies for each oligonucleotide
    sequence of length K in a sequence/s.
    Analyses that can be performed are:
        * PCA
        * t-SNE
        * Network
        * Phylogeny (Experimental)
        * Length

Args:
    freqDir (dir):
        Directory containing a list of files. Each file is compressed
        in Parquet format and contains the frequencies for each
        oligonucleotide sequence of length K in a sequence/s.

Returns:
    oDir (dir):
        Directory containing a list of files. Each file is TSV format and
        contains the results of an analysis.
'''

#------------------- Dependencies ---------------------------#

# Standard library imports
import argparse
import sys
from pathlib import Path

# External imports
import pandas as pd
import pyspark.sql.functions as F

# Internal imports
THIS_DIR = Path(__file__).resolve().parent
sys.path.append(str(THIS_DIR.parent))       ## Allow us to import from SRC
from src.util import spark
from src.kmer import analyse
from src.kmer import transform
from src.kmer.pairwise import distance as kmerdistance
from src import io
from src import kmer
from src import ml

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def runPca(iDirs, oDir, stranded):
    def _getTable(kmerDf):
        ## Convert the counts into Spark vectors and convert the dataframe into
        ## Spark RDD of (K, V) pairs in preparation for our analysis
        kmerDf = transform.agg.countsToSortedList(kmerDf)
        kmerPdfRdd = transform.toPdfRdd(kmerDf)

        ## Run PCA
        kmerId    = kmerPdfRdd.keys()
        kmerCount = kmerPdfRdd.values()
        kmerPca   = ml.feature.sparkIncrementalPcaReduce(kmerCount, n_components=3)

        ## Join tables
        kmerPca = kmerId.zip(kmerPca).map(lambda x: pd.concat(x, axis=1).fillna(''))
        cols    = list(kmerPca.first().columns)
        kmerPca = kmerPca.flatMap(lambda x: x.to_numpy().tolist()).toDF(cols)
        return kmerPca

    print("Running PCA")
    params = get_spark_params()
    with spark.getSparkSession(params) as ss:
        with ss.sparkContext as sc:
            ## Read Kmer frequencies, optimise parallelisation and QC
            kmerDf = io.kmer.read(*iDirs)
            kmerDf = qc(kmerDf, stranded)

            ## Get a table containing the PCA coordinates and
            kmerPca = _getTable(kmerDf)

            ## Write the table to disk
            print("Writing output")
            kmerPca.write.csv(oDir, mode='overwrite', sep='\t', header=True)
            ## Instead of outputing, we could extend the pipeline here
            ## so that it goes directly to analysing the PCA coordinates

def runTsne(iDirs, oDir, stranded):
    def _getTable(kmerDf):
        ## Convert the dataframe into Spark RDD of (K, V) pairs in preparation
        ## for our analysis
        kmerPdfRdd = transform.toPdfRdd(kmerDf)
        kmerPdfRdd = (
            kmerPdfRdd
            .map(transform.rotatePdf)
            .map(transform.splitPdf)
        )

        ## Load all the data into memory since we don't
        ## have a parallelised or incremental version
        kmerId    = kmerPdfRdd.keys()
        kmerCount = kmerPdfRdd.values()
        f = lambda x, y: pd.concat([x, y])
        kmerId    = kmerId.reduce(f).reset_index(drop=True)
        kmerCount = kmerCount.reduce(f).fillna(0).reset_index(drop=True)

        ## Run PCA and join tables
        kmerTsne  = ml.feature.sklearnReduce(kmerCount, 'TSNE', n_components=3)
        kmerTsne  = ml.feature.scale(kmerTsne)
        kmerTsne  = pd.concat([kmerId, kmerTsne], axis=1)
        return kmerTsne

    print("Running TSNE")
    params = get_spark_params()
    with spark.getSparkSession(params) as ss:
        with ss.sparkContext as sc:
            ## Read Kmer frequencies, optimise parallelisation and QC
            kmerDf = io.kmer.read(*iDirs)
            kmerDf = qc(kmerDf, stranded)

            ## Get a table containing the TSNE coordinates
            kmerTsne = _getTable(kmerDf)

            ## Write the table to disk
            print("Writing output")
            f = '{}.tsv'.format(oDir)
            kmerTsne.to_csv(f, sep='\t', index=False)
            ## Instead of outputing, we could extend the pipeline here
            ## so that it goes directly to analysing the PCA coordinates

def runNetwork(iDirs, oDir, stranded):
    def _getTable(kmerDf):
        ## Calculate the distance between each sequence and
        ## remove non-self pairs
        kmerDist = kmerdistance.traditional.getDistances(kmerDf)
        kmerDist = kmerDist.filter(kmerDist.seqId_1 != kmerDist.seqId_2)
        kmerDist = kmerdistance.scale(kmerDist)
        kmerDist = kmerDist.repartition(kmerDf.rdd.getNumPartitions())

        kmerId   = kmerDf.select(kmer.KMERDF_SEQINFO_COL_NAMES).distinct()
        kmerId   = kmerId.repartition(kmerDf.rdd.getNumPartitions())
        return (kmerDist, kmerId)

    print("Constructing network")
    params = get_spark_params()
    with spark.getSparkSession(params) as ss:
        with ss.sparkContext as sc:
            ## Read Kmer frequencies, optimise parallelisation and QC
            kmerDf = io.kmer.read(*iDirs)
            kmerDf = qc(kmerDf, stranded)

            ## Get tables of edges and nodes
            (kmerEdges, kmerNodes) = _getTable(kmerDf)

            ## Write the table to disk
            print("Writing output")
            e = '{}/edges'.format(oDir)
            n = '{}/nodes'.format(oDir)
            kmerEdges.write.csv(e, mode='overwrite', sep='\t', header=True)
            kmerNodes.write.csv(n, mode='overwrite', sep='\t', header=True)
            ## Instead of outputing, we could extend the pipeline here
            ## so that it goes directly to analysing the PCA coordinates

def runPhlogeny(iDirs, oDir, stranded):
    def _getDendrogram(kmerDf):
        ## There are 2 ways we can construct phylogenetic trees:
        ## * Based on PCA coordinates and distances
        ## * Directly based on distances
        ## My feeling is that they are basically the same. In fact, working directly
        ## on the distances should theoretically be better since we don't lose
        ## information. But regardless, we have to load all the data into memory
        ## since we don't have a parallelised way of doing this

        ## Approach 1.
        # kmerPca = runPCA(kmerDf)
        # kmerPca = kmerPca.toPandas()
        # kmerPca.index = kmerPca[kmer.SEQID_COL_NAME]
        # kmerPca = kmerPca[['Component1', 'Component2', 'Component3']]
        # ml.cluster.visualiseHierarchicalClusters(kmerPca, 0,
        #     linkage='average')
        # from matplotlib import pyplot as plt
        # plt.savefig("test.svg")

        ## Approach 2.
        ## Calculate the distance between each sequence and
        ## remove non-self pairs
        kmerDist = kmerdistance.traditional.getDistances(kmerDf)
        kmerDist = kmerDist.filter(kmerDist.seqId_1 != kmerDist.seqId_2)
        kmerDist = kmerdistance.scale(kmerDist)
        kmerDist = kmerDist.toPandas()
        kmerDist = kmerdistance.tableToSymMatrix(kmerDist)

        ## Draw plot
        ml.cluster.visualiseHierarchicalClusters(kmerDist, 0,
            linkage='average', affinity='precomputed')
        from matplotlib import pyplot as plt
        plt.savefig("test.svg")

    print("Drawing dendrogram")
    params = get_spark_params()
    with spark.getSparkSession(params) as ss:
        with ss.sparkContext as sc:
            ## Read Kmer frequencies, optimise parallelisation and QC
            kmerDf = io.kmer.read(*iDirs)
            kmerDf = qc(kmerDf, stranded)

            ## Export dendrogram of sequences
            ## Experimental...
            _getDendrogram(kmerDf)

def runLength(iDirs, oFile):
    print("Running K-mer length analysis")
    params = get_spark_params()
    with spark.getSparkSession(params) as ss:
        with ss.sparkContext as sc:
            ## Read Kmer frequencies, optimise parallelisation and QC
            kmerDf = io.kmer.read(*iDirs)

            ## Calculate the average K-mer length based on the length
            ## of each sequence
            avgLength = analyse.length.getLowerLimit(kmerDf)
            print("Average K-mer length:\t{}".format(avgLength))

            acf = analyse.length.getAcf(kmerDf)
            print("Average number of common K-mers:\t{}".format(acf))

            fck = analyse.length.getFck(kmerDf)
            print("Fraction of common K-mers:\t{}".format(fck))

            fuk = analyse.length.getFuk(kmerDf)
            print("Fraction of unique K-mers:\t{}".format(fuk))

def runPlot(iDirs, oFile, stranded, seqIds, pType):
    def _getPlot(kmerDf):
        ## Convert the dataframe into Spark RDD of (K, V) pairs in preparation
        ## for our analysis
        kmerPdfRdd = transform.toPdfRdd(kmerDf)
        kmerPdfRdd = (
            kmerPdfRdd
            .map(transform.rotatePdf)
            .map(transform.splitPdf)
        )

        ## Load all the data into memory since we don't
        ## have a parallelised or incremental version
        kmerId    = kmerPdfRdd.keys()
        kmerCount = kmerPdfRdd.values()
        f = lambda x, y: pd.concat([x, y])
        kmerId    = kmerId.reduce(f).reset_index(drop=True)
        kmerCount = kmerCount.reduce(f).fillna(0).reset_index(drop=True)
        kmerCount.index = kmerId[kmer.SEQID_COL_NAME]

        ## Draw plot
        from matplotlib import pyplot as plt
        if (pType):
            analyse.viz.histogram(kmerCount)
            plt.savefig("test.svg")

        else:
            analyse.viz.radar(kmerCount)
            plt.savefig("test.svg")

    print("Drawing K-mer frequencies plot")
    params = get_spark_params()
    with spark.getSparkSession(params) as ss:
        with ss.sparkContext as sc:
            ## Read Kmer frequencies, optimise parallelisation and QC
            kmerDf = io.kmer.read(*iDirs)
            kmerDf = qc(kmerDf, stranded)
            if (seqIds is not None):
                seqIds = '|'.join(seqIds)
                cond   = (F.col(kmer.SEQID_COL_NAME).rlike(seqIds))
                kmerDf = kmerDf.filter(cond)

            ## Export plot of counts
            _getPlot(kmerDf)

#------------------- Private Classes & Functions ------------#

def qc(kmerDf, stranded):
    ## Repartition when we have a relatively small number of sequences
    # nParts = 2

    ## **********
    ## *** Re-organise Kmer frequencies depending on what sequences we're
    ## *** analysing. This might involve summing counts of complementary
    ## *** Kmers (if sequences are not in the same orientation)
    ## **********
    if (stranded):
        kmerDf = transform.agg.clusterRevCompCounts(kmerDf)

    ## Normalise frequencies
    kmerDf = transform.convert.countsToNormalised(kmerDf)
    # kmerDf = transform.convert.countsToProbabilities(kmerDf)
    return kmerDf

def get_spark_params():
    # import os                                  ## For HPC
    params = [
        ## Driver
        # ('spark.driver.cores', '5'),           ## Same as executor
        ('spark.driver.memory', '27G'),        ## Same as executor
        ('spark.driver.maxResultSize', '0'),

        ## Executor
        # ('spark.executor.cores', '5'),         ## Didn't seem to have any effect...
        # ('spark.executor.instances', '99'),    ## Didn't seem to have any effect...
        ('spark.executor.memory', '27G'),
        ('spark.executor.heartbeatInterval', '60s'),

        ## SQL
        ('spark.sql.broadcastTimeout', '600s'),
        ('spark.sql.execution.arrow.pyspark.enabled', 'true'),
        ('spark.sql.shuffle.partitions', '200'),

        ## Misc
        ('spark.local.dir', THIS_DIR),
        # ('spark.local.dir', os.environ['MEMDIR']), ## For HPC
        ('spark.network.timeout', '600s'),
        # ('spark.default.parallelism', '8')
    ]
    return params

def make_parser():
    def _addPlotArgs(p):
        p.add_argument("-i", help="IDs", nargs='+', type=str)
        p.add_argument("-t", help="Plot histogram", action='store_true')
        _initArgs(p)

    def _initArgs(p):
        p.add_argument("--freqDir", help="Observed Kmer frequencies directory",
            nargs='*', type=str, required=True)
        p.add_argument("-o", help="Output file/directory",
            type=str, required=True)
        p.add_argument("-d", help="Analyse sequences as double-stranded \
            instead of single-stranded", action='store_true')

    parser    = argparse.ArgumentParser(description='Analyse \
        oligonucleotide frequencies')
    subparser = parser.add_subparsers(dest='command')
    pca     = subparser.add_parser('pca')
    tsne    = subparser.add_parser('tsne')
    network = subparser.add_parser('network')
    phylo   = subparser.add_parser('phylo')
    length  = subparser.add_parser('length')
    plot    = subparser.add_parser('plot')
    _initArgs(pca)
    _initArgs(tsne)
    _initArgs(network)
    _initArgs(phylo)
    _initArgs(length)
    _addPlotArgs(plot)
    return parser

def main(parser):
    args = parser.parse_args()
    if (args.command is None):
        parser.print_help()
        sys.exit(1)

    else:
        if (args.command == 'pca'):
            print(args)
            runPca(args.freqDir, args.o, args.d)

        elif (args.command == 'tsne'):
            print(args)
            runTsne(args.freqDir, args.o, args.d)

        elif (args.command == 'network'):
            print(args)
            runNetwork(args.freqDir, args.o, args.d)

        elif (args.command == 'phylo'):
            print(args)
            runPhylogeny(args.freqDir, args.o, args.d)

        elif (args.command == 'length'):
            print(args)
            runLength(args.freqDir, args.d)

        elif (args.command == 'plot'):
            print(args)
            runPlot(args.freqDir, args.o, args.d, args.i, args.t)

    print("DONE")

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    parser = make_parser()
    main(parser)

#------------------------------------------------------------------------------
