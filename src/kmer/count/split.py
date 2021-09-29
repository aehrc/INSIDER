#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import math
from operator import add

# External imports
from pyspark.sql import SparkSession

# Internal imports
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
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    rdd = rdd.map(lambda x: (x.id, str(x.seq).upper()))
    if (countExp):
        print("Expected counts")
        kImerLength  = kmerLength - 1
        kIImerLength = kmerLength - 2

        kImerDf     = _getObservedCounts(rdd, kImerLength, ignoreNs)
        seqLengthDf = _getSeqLengths(rdd)
        kImerDf     = _createKmerDf(kImerDf, seqLengthDf)
        kImerDf     = kImerDf.withColumnRenamed(KMER_COL_NAME, 'kImer') \
            .withColumnRenamed('count', 'kImerCount')

        kIImerDf    = _getObservedCounts(rdd, kIImerLength, ignoreNs)
        seqLengthDf = _getSeqLengths(rdd)
        kIImerDf    = _createKmerDf(kIImerDf, seqLengthDf)
        kIImerDf    = kIImerDf.withColumnRenamed(KMER_COL_NAME, 'kIImer') \
            .withColumnRenamed('count', 'kIImerCount')

        kmerDf = getExpectedCounts(kImerDf, kIImerDf, kImerLength, kIImerLength)

    else:
        print("Observed counts")
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
        kmerDf      = _getObservedCounts(rdd, kmerLength, ignoreNs)
        seqLengthDf = _getSeqLengths(rdd)
        kmerDf      = _createKmerDf(kmerDf, seqLengthDf)

    ## Group (K, V) pairs. This makes reading the files a bit easier.
    kmerDf  = transform.agg.countsToDict(kmerDf)
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
    ## (ID, Seq) => (ID, kmerSeq)
    ##             => ((ID, kmerSeq), 1)
    ##             => ((ID, kmerSeq), count)
    f = lambda x: getObservedSequences(x, kmerLength)
    kmerRdd = rdd.flatMapValues(f).map(lambda x: (x, 1)) \
        .reduceByKey(add).map(lambda x: (*x[0], x[1]))
    ## We only count Kmers that occur at least once. Kmers that
    ## do not occur (zero counts) are ignored, but need to be
    ## added when we're analysing results. This will save us quite
    ## a bit of space.
    kmerRdd.cache()
    return kmerRdd

def _getSeqLengths(rdd):
    ## (ID, Seq) => (ID, len(Seq))
    seqLengthRdd = rdd.mapValues(len)
    seqLengthDf  = seqLengthRdd.toDF([SEQID_COL_NAME, SEQLEN_COL_NAME])
    return seqLengthDf

def _createKmerDf(kmerDf, seqLengthDf):
    colNames = [SEQID_COL_NAME, SEQLEN_COL_NAME, KMER_COL_NAME, COUNT_COL_NAME]
    kmerDf   = kmerDf.join(seqLengthDf, SEQID_COL_NAME, 'left').select(colNames)
    return kmerDf

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
