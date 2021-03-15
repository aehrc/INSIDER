#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import math
import hashlib
from operator import add

# External imports
import pyspark.sql.functions as sparkF
from pyspark.sql import SparkSession

# Internal imports
from .constants import *
from .common import *
from .. import transform

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def setup(rdd, kmerLength):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Partition sequences. This should be proportional to
    ## the number of sequence records we're querying.
    sc     = ss.sparkContext
    eParts = math.floor(rdd.count() / CHUNK_SIZE)
    nParts = sc.defaultParallelism + (eParts * N_PARTS_MULTI)
    nParts = nParts * kmerLength
    return rdd.repartition(nParts)

def getCounts(rdd, kmerLength, ignoreNs, countExp):
    ss  = SparkSession.getActiveSession()
    rdd = rdd.map(lambda x: (x.id, str(x.seq).upper()))

    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    if (countExp == 'ZOM'):
        raise NotImplementedError('ZOM counts not implemented!')

    elif (countExp == 'MCM'):
        ## P(W) = [P(W[0:K-1]) * P(W[1:K])] / P(W[1:K-1])
        if (kmerLength < 3):
            raise ValueError('Cannot calculate counts for kmerLength < 3')

        kImerLength  = kmerLength - 1
        kIImerLength = kmerLength - 2

        kImerDf   = _getObservedCounts(rdd, kImerLength, ignoreNs)
        totalDf   = _getTotalCounts(kImerDf.rdd)
        kImerDf   = _createKmerDf(kImerDf, totalDf, asCounts=False)

        kIImerDf  = _getObservedCounts(rdd, kIImerLength, ignoreNs)
        totalDf   = _getTotalCounts(kIImerDf.rdd)
        kIImerDf  = _createKmerDf(kIImerDf, totalDf, asCounts=False)

        kmerDf    = kImerDf.groupby(kImerDf.id) \
                           .cogroup(kIImerDf.groupby(kIImerDf.id)) \
                           .applyInPandas(getMCMCounts, schema=kImerDf.schema)

    else:
        ## **********
        ## *** Not sure if this is the "best" approach for split counts.
        ## *** But seems to run the fastest!
        ## ***
        ## *** Using RDD ReduceByKey once to minimise shuffling can encounter
        ## *** a few errors (not sure why).
        ## ***
        ## *** I've tried:
        ## ***  * Using fewer lambda functions to reduce double serialisation.
        ## ***  * Repartitioning before reducing. However, this doesn't seem to
        ## ***    result in any major improvements since we're shuffling data twice.
        ## **********
        kmerDf   = _getObservedCounts(rdd, kmerLength, ignoreNs)
        totalDf  = _getTotalCounts(kmerDf.rdd)
        kmerDf   = _createKmerDf(kmerDf, totalDf)

    ## Group (K, V) pairs. This makes reading the files a bit easier.
    kmerDf = transform.agg.countsToDict(kmerDf)
    return kmerDf

def cleanup(kmerDf, kmerLength):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Calculate the expected number of partitions/files based on the
    ## data. This should be proportional
    ## to the number of sequence records and the kmer length.
    nSeqs  = kmerDf.select(SEQID_COL_NAME).distinct().count()
    eParts = math.ceil(nSeqs / CHUNK_SIZE)
    eParts = eParts + math.floor((kmerLength / 5))
    eParts = min(MAX_N_FILES, eParts)

    ## Compare the expected and observed number of partitions/files.
    ## Repartition the data accordingly.
    oParts = kmerDf.rdd.getNumPartitions()
    kmerDf = kmerDf.coalesce(eParts) if oParts > eParts else kmerDf.repartition(eParts)
    return kmerDf

#------------------- Private Classes & Functions ------------#

def _getObservedCounts(rdd, kmerLength, ignoreNs):
    kmerRdd  = _getKmerRdd(rdd, kmerLength)
    kmerDf   = kmerRdd.toDF([SEQID_COL_NAME, KMER_COL_NAME, COUNT_COL_NAME])
    if (ignoreNs):
        kmerDf   = getValidRows(kmerDf)

    return kmerDf

def _getKmerRdd(rdd, kmerLength):
    ## Get the Kmer counts for each unique record
    ## Counts are (K, V) pairs, where K = sequence and V = frequency
    ## (Hash, Seq) => (Hash, kmerSeq)
    ##             => ((Hash, kmerSeq), 1)
    ##             => ((Hash, kmerSeq), count)
    f = lambda x: getObservedSequences(x, kmerLength)
    kmerRdd = rdd.flatMapValues(f) \
                 .map(lambda x: (x, 1)) \
                 .reduceByKey(add) \
                 .map(lambda x: (*x[0], x[1]))
    ## We only count Kmers that occur at least once. Kmers that
    ## do not occur (zero counts) are ignored, but need to be
    ## added when we're analysing results. This will save us quite
    ## a bit of space.
    kmerRdd.cache()
    return kmerRdd

def _getTotalCounts(rdd):
    ## Get the total number of Kmer for each records
    ## (Hash, kmerSeq, count) => (Hash, count)
    ##                        => (Hash, total)
    f = lambda x: (x[0], x[2])
    totalRdd = rdd.map(f).reduceByKey(add)
    totalDf  = totalRdd.toDF([SEQID_COL_NAME, TOTAL_COL_NAME])
    return totalDf

def _createKmerDf(kmerDf, totalDf, asCounts=True):
    if (asCounts):
        kmerDf   = kmerDf.join(sparkF.broadcast(totalDf), SEQID_COL_NAME, 'left') \
                         .select(kmerDf.id, TOTAL_COL_NAME,
                                 KMER_COL_NAME, COUNT_COL_NAME)

    else:
        f      = transform.convert.countsToProbabilities(COUNT_COL_NAME, TOTAL_COL_NAME)
        kmerDf = kmerDf.join(sparkF.broadcast(totalDf), SEQID_COL_NAME, 'left') \
                       .select(kmerDf.id, TOTAL_COL_NAME,
                               KMER_COL_NAME, f.alias(COUNT_COL_NAME))

    return kmerDf

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
