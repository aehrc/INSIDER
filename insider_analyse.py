#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import argparse
import sys
from pathlib import Path 

# External imports
import pandas as pd
import pyspark.sql.functions as sparkF
from pyspark.sql import SparkSession

# Internal imports
THIS_DIR = Path(__file__).resolve().parent
sys.path.append(str(THIS_DIR.parent))       ## Allow us to import from SRC
from src.util import spark
from src.kmer import analyse
from src.kmer import transform
from src.kmer.pairwise import stats as kmerstats
from src.kmer.pairwise import distance as kmerdistance
from src import io
from src import kmer
from src import ml

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def mainKmerFrequencies(iDirs, oFile, cIdFile):
    def _setup(kmerDf, cIdFile):
        ## Calculate the genome signature by summing the counts of ALL sequences
        kmerDfY = transform.agg.sumCounts(kmerDf)
        kmerDfY = transform.agg.clusterRevCompCounts(kmerDfY)

        ## Calculate the cluster/contig signatures of relatively longer sequences
        kmerDfX = kmerDf if cIdFile is None else _getClusterDf(kmerDf, cIdFile)
        kmerDfX = transform.filter.removeShortSequences(kmerDfX)
        kmerDfX = transform.agg.clusterRevCompCounts(kmerDfX)
        return (kmerDfX, kmerDfY)

    def _getClusterDf(kmerDf, cIdFile):
        cIdDf = readClusterLabels(cIdFile)

        ## Cluster Kmer frequencies
        cKmerDf = transform.agg.clusterCountsByColumn(kmerDf, cIdDf)
        cKmerDf = cKmerDf.repartition(kmerDf.rdd.getNumPartitions(),
            kmer.SEQID_COL_NAME)

        print("Calculating cluster info")
        csDf  = analyse.cluster.getClusterSizes(cKmerDf, cIdDf)
        cmcDf = analyse.cluster.getClusterMeanCounts(cKmerDf, cIdDf)
        return cKmerDf

    def _analyse(kmerDfX, kmerDfY):
        def _runStats(kmerDfX, kmerDfY):
            ## Calculate Pvalue
            print("Running Chi-squared test")
            kmerDfY   = getExpectedCounts(kmerDfX, kmerDfY)
            kmerStats = kmerstats.chi.getPvalues(kmerDfX, kmerDfY)
            kmerStats = kmerStats.filter(kmerStats.seqId_1 == kmerStats.seqId_2)
            kmerStats = kmerstats.getZScore(kmerStats, 'cramers_v')
            kmerStats = kmerStats.drop('seqId_2') \
                .withColumnRenamed('seqId_1', kmer.SEQID_COL_NAME) \
                .withColumnRenamed('p_value', 'chi_p_value')
            return kmerStats

        def _runPCA(kmerDfX, kmerDfY):
            ## Get PCA coordinates
            print("Running PCA")
            nParts  = kmerDfX.rdd.getNumPartitions()
            kmerDfX = kmerDfX.unionByName(kmerDfY)
            kmerDfX = transform.convert.countsToNormalised(kmerDfX)
            kmerDfX = transform.agg.countsToSortedList(kmerDfX)
            kmerDfX = kmerDfX.coalesce(2)

            ## Run PCA
            kmerPdfRdd = transform.toPdfRdd(kmerDfX)
            kmerId     = kmerPdfRdd.keys()
            kmerCount  = kmerPdfRdd.values()
            kmerPca    = ml.feature.sparkIncrementalPcaReduce(kmerCount, n_components=3)

            ## Join tables and repartition the dataframe
            kmerPca = kmerId.zip(kmerPca).map(lambda x: pd.concat(x, axis=1).fillna(''))
            cols    = list(kmerPca.first().columns)
            kmerPca = kmerPca.flatMap(lambda x: x.to_numpy().tolist()).toDF(cols)
            return kmerPca

        ## Run our analysis
        kmerStats   = _runStats(kmerDfX, kmerDfY)
        kmerPca     = _runPCA(kmerDfX, kmerDfY)
        kmerOutlier = analyse.outlier.getLabels(kmerDfX)

        ## Join the above tables. We shouldn't have any
        ## problems if each table is relatively small.
        kmerDf = kmerPca.join(kmerOutlier, kmer.SEQID_COL_NAME, 'left') \
            .join(kmerStats, kmer.SEQID_COL_NAME, 'left')
        kmerDf = kmerDf.fillna(0)
        return kmerDf

    def _cleanup(kmerDf, cIdFile):
        if (cIdFile is not None):
            cIdDf = readClusterLabels(cIdFile)
            kmerDf = kmerDf.withColumnRenamed(kmer.SEQID_COL_NAME, 'clusterId')
            kmerDf = kmerDf.join(cIdDf, 'clusterId', 'left')

        kmerDf = kmerDf.orderBy('chi_p_value')
        return kmerDf

    print("Non-reference frequencies analysis")
    params = get_spark_params()
    with spark.getSparkSession(params) as ss:
        with ss.sparkContext as sc:
            ## Read Kmer frequencies
            print("Reading Kmers")
            kmerDf = io.kmer.read(*iDirs)

            ## Setup our DataFrames:
            ## * kmerDfY = Frequencies across all sequences
            ## * kmerDfX = Frequencies of sequences/clusters
            ##             that meet our filtering criteria
            (kmerDfX, kmerDfY) = _setup(kmerDf, cIdFile)

            ## Run our analysis
            kmerDf = _analyse(kmerDfX, kmerDfY)

            ## Cleanup our DataFrame for writing
            kmerDf = _cleanup(kmerDf, cIdFile)

            ## Write the table to disk
            print("Writing output")
            kmerDf  = kmerDf.coalesce(sc.defaultParallelism)
            kmerDf.write.csv(oFile, mode='overwrite', sep='\t', header=True)
            ## Instead of outputing, we could extend the pipeline here
            ## so that it goes directly to analysing the PCA coordinates

def newKmerFrequencies(iDirsX, iDirsY, oFile, cIdFile):
    def _setup(kmerDfX, kmerDfY):
        ## Calculate the genome signature by summing the counts of ALL sequences
        nParts  = kmerDfY.rdd.getNumPartitions()
        kmerDfY = transform.agg.sumCounts(kmerDfY)
        kmerDfY = transform.agg.clusterRevCompCounts(kmerDfY)
        kmerDfY = transform.convert.countsToNormalised(kmerDfY)
        kmerDfY = kmerDfY.repartition(nParts, kmer.SEQID_COL_NAME)

        ## Calculate the contig signatures of relatively longer sequences
        nParts  = kmerDfX.rdd.getNumPartitions()
        kmerDfX = transform.filter.removeShortSequences(kmerDfX)
        kmerDfX = transform.agg.clusterRevCompCounts(kmerDfX)
        kmerDfX = transform.convert.countsToNormalised(kmerDfX)
        kmerDfX = kmerDfX.repartition(nParts, kmer.SEQID_COL_NAME)
        return (kmerDfX, kmerDfY)

    def _analyse(kmerDfX, kmerDfY):
        ## Calculate the distance between each contig signature
        ## and the genome signature
        print("Calculating distance")
        ## We don't have to use a distance. In fact, we could use a chi2 again.
        kmerDist = kmerdistance.traditional.getDistances(kmerDfX, kmerDfY)
        kmerDist = kmerdistance.scale(kmerDist)
        kmerDist = kmerstats.getZScore(kmerDist, 'distance')
        kmerDist = kmerDist.drop('seqId_2') \
            .withColumnRenamed('seqId_1', kmer.SEQID_COL_NAME)

        ## Join the above tables
        kmerDf  = kmerDfX.select(kmer.KMERDF_SEQINFO_COL_NAMES).distinct() \
            .join(kmerDist, kmer.SEQID_COL_NAME, 'left')
        return kmerDf

    print("Reference frequencies analysis")
    params = get_spark_params()
    with spark.getSparkSession(params) as ss:
        with ss.sparkContext as sc:
            ## Read Kmer frequencies
            print("Reading Kmers")
            kmerDfX = io.kmer.read(*iDirsX)
            kmerDfY = io.kmer.read(*iDirsY)

            ## Set up our DataFrames:
            ## * kmerDfY = Frequencies across all sequences
            ## * kmerDfX = Frequencies of sequences/clusters
            ##             that meet our filtering criteria
            (kmerDfX, kmerDfY) = _setup(kmerDfX, kmerDfY)

            ## Run our analysis
            kmerDf = _analyse(kmerDfX, kmerDfY)

            ## Write the table to disk
            print("Writing output")
            kmerDf = kmerDf.coalesce(sc.defaultParallelism)
            kmerDf.write.csv(oFile, mode='overwrite', sep='\t', header=True)

#------------------- Private Classes & Functions ------------#

def readClusterLabels(cIdFile):
    ss  = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Read cluster labels
    print("Reading cluster labels")
    cIdDf   = ss.read.csv(cIdFile, sep='\t', header=True)
    cIdDf   = cIdDf.select(kmer.SEQID_COL_NAME, 'clusterId')
    return cIdDf

def getExpectedCounts(kmerDfX, kmerDfY):
    def _f(iterator):
        for pdf in iterator:
            pdf[kmer.COUNT_COL_NAME] = pdf[kmer.COUNT_COL_NAME] * pdf['t']
            pdf = pdf.drop(columns=['t'])
            yield pdf

    nParts  = kmerDfX.rdd.getNumPartitions()
    kmerDfX = kmerDfX.groupby(kmer.KMERDF_SEQINFO_COL_NAMES) \
        .agg(sparkF.sum(kmer.COUNT_COL_NAME).alias('t'))

    ## Calculate the probabilities for each Kmer
    kmerDf = transform.convert.countsToProbabilities(kmerDfY)
    kmerDf = kmerDf.alias('r').crossJoin(kmerDfX.alias('l'))
    kmerDf = kmerDf.select(sparkF.col('l.seqId').alias(kmer.SEQID_COL_NAME),
        sparkF.col('l.seqLen').alias(kmer.SEQLEN_COL_NAME),
        sparkF.col('l.filename').alias(kmer.FILE_COL_NAME),
        sparkF.col('r.kmer').alias(kmer.KMER_COL_NAME),
        sparkF.col('r.count').alias(kmer.COUNT_COL_NAME),
        sparkF.col('t').alias('t'))

    kmerDfY = kmerDf.mapInPandas(_f, schema=kmerDfY.schema)
    kmerDfY = kmerDfY.repartition(kmerDfY.rdd.getNumPartitions())
    return kmerDfY

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
    def _addNewArgs(p):
        p.add_argument("--refDir", help="Observed Kmer frequencies directory",
            nargs='*', type=str, required=True)
        _initArgs(p)

    def _initArgs(p):
        p.add_argument("--freqDir", help="Observed Kmer frequencies directory",
            nargs='*', type=str, required=True)
        p.add_argument("-o", help="Output file (.tsv)",
            type=str, required=True)
        p.add_argument("--cIdFile", help="Cluster labels",
            type=str)

    parser    = argparse.ArgumentParser(description='Run INSIDER analysis')
    subparser = parser.add_subparsers(dest='command')
    main      = subparser.add_parser('main')
    new       = subparser.add_parser('new')
    _initArgs(main)
    _addNewArgs(new)
    return parser

def main(parser):
    args = parser.parse_args()
    if (args.command is None):
        parser.print_help()
        sys.exit(1)

    else:
        if (args.command == 'main'):
            print(args)
            mainKmerFrequencies(args.freqDir, args.o, args.cIdFile)

        if (args.command == 'new'):
            print(args)
            newKmerFrequencies(args.freqDir, args.refDir, args.o, args.cIdFile)

    print("DONE")

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    parser = make_parser()
    main(parser)

#------------------------------------------------------------------------------
