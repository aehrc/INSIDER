#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import math
import random
from operator import add

# External imports
import pyspark.sql.functions as F
import pyspark.sql.types as T
from pyspark.sql import SparkSession

# Internal imports
from .common import *
from .. import transform

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def setup(rdd):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Partition sequences. This should be proportional to
    ## the number of sequence records we're querying.
    sc     = ss.sparkContext
    eParts = math.floor(rdd.count() / CHUNK_SIZE)
    nParts = sc.defaultParallelism + (eParts * N_PARTS_MULTI)
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

        kImerDf   = _getObservedCounts(rdd, kImerLength, ignoreNs)
        seqLength = _getSeqLength(rdd)
        kImerDf   = _createKmerDf(kImerDf, seqLength)
        kImerDf   = kImerDf.withColumnRenamed(KMER_COL_NAME, 'kImer') \
            .withColumnRenamed('count', 'kImerCount')

        kIImerDf  = _getObservedCounts(rdd, kIImerLength, ignoreNs)
        seqLength = _getSeqLength(rdd)
        kIImerDf  = _createKmerDf(kIImerDf, seqLength)
        kIImerDf  = kIImerDf.withColumnRenamed(KMER_COL_NAME, 'kIImer') \
            .withColumnRenamed('count', 'kIImerCount')

        kmerDf = getExpectedCounts(kImerDf, kIImerDf, kImerLength, kIImerLength)

    else:
        print("Observed counts")
        ## **********
        ## *** This is (probably) the "best" approach for combined counts.
        ## *** RDD ReduceByKey runs better than DataFrame GroupBy.
        ## *** So DO NOT change!!!
        ## ***
        ## *** Using RDD ReduceByKey once to minimise shuffling can encounter
        ## *** a few errors (not sure why).
        ## ***
        ## *** My understanding is that DataFrames GroupBy is faster, but
        ## *** I don't understand why this isn't the case. I've tried
        ## ***  * Explicitly defining the schema
        ## ***  * Repartitioning the DataFrame before grouping
        ## **********
        kmerDf    = _getObservedCounts(rdd, kmerLength, ignoreNs)
        seqLength = _getSeqLength(rdd)
        kmerDf    = _createKmerDf(kmerDf, seqLength)

    ## Group (K, V) pairs. This makes reading the files a bit easier.
    kmerDf = _pivot(kmerDf, kmerLength)
    return kmerDf

def cleanup(kmerDf, kmerLength):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Calculate the expected number of partitions/files based on the
    ## data. This should be proportional to the number of sequence records and
    ## the kmer length
    eParts = 1
    if (kmerLength > 10):
        ## Start writing to more files
        ## once we start reaching > 10-mers
        r = kmerLength % 10
        eParts = 2 ** math.ceil((r / 2))
    eParts = min(MAX_N_FILES, eParts)

    ## Compare the expected and observed number of partitions/files.
    ## Repartition the data accordingly.
    oParts = kmerDf.rdd.getNumPartitions()
    kmerDf = kmerDf.coalesce(eParts) if oParts > eParts else kmerDf.repartition(eParts)
    return kmerDf

#------------------- Private Classes & Functions ------------#

def _getObservedCounts(rdd, kmerLength, ignoreNs):
    kmerRdd = _getKmerRdd(rdd, kmerLength)
    kmerDf  = kmerRdd.toDF([SEQID_COL_NAME, KMER_COL_NAME, COUNT_COL_NAME])
    if (ignoreNs):
        kmerDf = getValidRows(kmerDf)

    return kmerDf

def _getKmerRdd(rdd, kmerLength):
    ## Get the Kmer counts across all records
    ## Counts are (K, V) pairs, where K = sequence and V = frequency
    ## (ID, Seq) => (ID, kmerSeq)
    ##           => (kmerSeq, 1)
    ##           => (kmerSeq, count)
    f = lambda x: getObservedSequences(x, kmerLength)
    kmerRdd = rdd.flatMapValues(f).map(lambda x: (x[1], 1)) \
        .reduceByKey(add).map(lambda x: (0, *x))
    ## We only count Kmers that occur at least once. Kmers that
    ## do not occur (zero counts) are ignored, but need to be
    ## added when we're analysing results. This will save us quite
    ## a bit of space.
    kmerRdd.cache()
    return kmerRdd

def _getSeqLength(rdd):
    ## (ID, Seq) => (ID, len(Seq))
    seqLength = rdd.map(lambda x: len(x[1])).reduce(add)
    return seqLength

def _createKmerDf(kmerDf, seqLength):
    colNames = [SEQID_COL_NAME, F.lit(seqLength).alias(SEQLEN_COL_NAME),
                KMER_COL_NAME, COUNT_COL_NAME]
    kmerDf   = kmerDf.select(colNames)
    return kmerDf

def _pivot(kmerDf, kmerLength):
    ## Writing large Kmers (i.e., > 10) to file encounters a few warnings
    ## because of the large number of (K, V) pairs. To make it a bit easier
    ## to read and write, spread the columns across multiple rows.
    ## Just make sure that we combine the rows correctly
    ## when we're analysing results.

    ## Keep track of columns we want
    cols = kmerDf.schema.names
    cols.remove(SEQID_COL_NAME)

    ## Replace the ID column to something
    ## that can be used for grouping
    ## (ID, Kmer, Percent, *) => (ID, [Kmer -> Percent], *)
    f      = F.round(F.rand()*(kmerLength-1))
    kmerDf = kmerDf.select(f.alias(SEQID_COL_NAME), *cols)
    kmerDf = transform.agg.countsToDict(kmerDf)

    ## Replace the ID column to something
    ## more suitable for Identification purposes
    cols.remove(COUNT_COL_NAME)
    f      = F.lit(random.randint(0, 100000))
    kmerDf = kmerDf.select(f.alias(SEQID_COL_NAME).cast(T.StringType()), *cols)
    return kmerDf

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
