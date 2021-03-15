#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import pandas as pd
from pyspark.sql import SparkSession
import pyspark.sql.functions as sparkF

# Internal imports
from ..constants import *
from .. import transform
from ... import ml

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def setup(kmerDf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

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
    kmerDf = kmerDf.groupby(kmerDf.id) \
        .applyInPandas(transform.agg.clusterRevCompCounts, schema=kmerDf.schema)

    ## **********
    ## *** Perform other filtering / processing steps including:
    ## ***   * Repetitive (i.e., low complexity) Kmers
    ## **********
    kmerDf = transform.filter.removeRepetitiveKmers(kmerDf)

    ## **********
    ## *** Filtering Kmers will ultimately change the counts, so we need to
    ## *** adjust the total counts accordingly. However, when aggregating Kmers
    ## *** the total counts are adjusted automatically
    ## **********
    f = lambda x: x.assign(total=x[COUNT_COL_NAME].sum())
    kmerDf = kmerDf.groupby(SEQID_COL_NAME) \
        .applyInPandas(f, schema=kmerDf.schema)

    kmerDf.persist()
    kmerDf.show()
    return kmerDf

def getPcaDf(kmerDf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## To estimate what the true signature might be, sum the counts
    ## of each K-mer across all non-overlapping sequences
    eKmerDf    = transform.agg.sumCounts(kmerDf)
    kmerDf     = kmerDf.unionByName(eKmerDf)
    kmerDf     = kmerDf.repartition(kmerDf.rdd.getNumPartitions(),
        SEQID_COL_NAME)

    ## Convert to Pandas
    kmerPdfRdd = transform.toPdfRdd(kmerDf)
    kmerPdfRdd = kmerPdfRdd.map(transform.rotatePdf).map(transform.splitPdf)
    kmerId     = kmerPdfRdd.keys()
    kmerCount  = kmerPdfRdd.values()

    ## This will load all the data into memory
    kmerId    = kmerId.reduce(lambda x, y: pd.concat([x, y])) \
        .reset_index(drop=True)
    kmerCount = kmerCount.reduce(lambda x, y: pd.concat([x, y])) \
        .fillna(0).reset_index(drop=True)

    kmerCount = kmerCount.divide(kmerCount.sum(axis=1), axis=0)
    kmerPca   = ml.feature.sklearnReduce(kmerCount, 'PCA', n_components=0.75)
    ## Ensure that we have at least 2 components
    if (kmerPca.shape[0] == 1):
        kmerPca   = ml.feature.sklearnReduce(kmerCount, 'PCA', n_components=0.75)
    kmerPca   = ml.feature.scale(kmerPca)

    kmerPca   = pd.concat([kmerId.reset_index(drop=True), kmerPca], axis=1)
    return kmerPca

#------------------- Protected Classes & Functions ------------#

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
