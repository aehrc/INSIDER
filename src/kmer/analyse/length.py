#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import math

# External imports
import numpy as np
import pyspark.sql.functions as F
import pyspark.sql.types as T
from pyspark.sql import SparkSession

# Internal imports
from ..common import *
from ..pairwise import getPairs

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getLowerLimit(kmerDf):
    ## Based on formulas reported in:
    ## * Luczak et al. (2018)
    ## * Sims et al. (2009)
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Calculate the K-mer length based on the length of each sequence
    seqLenDf  = kmerDf.select(SEQID_COL_NAME, SEQLEN_COL_NAME).distinct()
    lengthRdd = seqLenDf.rdd.mapValues(lambda x: math.log(x, 4))

    ## Calculate the average K-mer length
    lengths   = lengthRdd.values().collect()
    avgLength = np.mean(lengths)
    return avgLength

def getAcf(kmerDf):
    ## Calculate the average number of common Kmers between sequences
    ## Based on the formulas reported in:
    ## * Zhang et al. (2017)
    ## * Pornputtapong et al. (2020)
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Find all non-self pairs
    numSeqs = kmerDf.select(SEQID_COL_NAME).distinct().count()
    kmerDf  = getPairs(kmerDf)
    kmerDf  = kmerDf.filter(F.col('l.rId') < F.col('r.rId'))

    ## Calculate the average number of common Kmers between each non-self pair
    f = F.size(F.array_intersect('l.kmer', 'r.kmer')) / (F.lit(numSeqs - 1))
    kmerDf = kmerDf.select(f.alias('commonKmers'))
    acf    = kmerDf.groupby().sum().collect()[0][0]
    return acf

def getFck(kmerDf):
    ## Calculate the fraction of Kmers present in all sequences
    ## Based on the approach reported in:
    ## * Gardner et al. (2013)
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    numSeqs = kmerDf.select(SEQID_COL_NAME).distinct().count()
    kmerDf  = kmerDf.groupby(KMER_COL_NAME).count()

    ## Calculate the number of Kmers present in all sequences
    cKmers = kmerDf.filter(F.col(COUNT_COL_NAME) == numSeqs).count()
    uKmers = kmerDf.filter(F.col(COUNT_COL_NAME) != numSeqs).count()
    fck    = cKmers / (cKmers + uKmers)
    return fck

def getFuk(kmerDf):
    ## Calculate the average fraction of unique Kmers of each sequences
    ## Based on the approach reported in:
    ## * Gardner et al. (2013)
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Count the total number of Kmers for each sequence
    f = F.sum(COUNT_COL_NAME)
    seqTotals = kmerDf.groupby(SEQID_COL_NAME).agg(f.alias('t'))

    ## Count the number of unique Kmers in each sequence
    uniqKmers = (
        kmerDf
        .filter(F.col(COUNT_COL_NAME) == 1)
        .groupby(SEQID_COL_NAME).agg(f.alias(COUNT_COL_NAME))
    )

    ## Calculate the fraction of unique Kmers in each sequence
    f = F.col(COUNT_COL_NAME) / F.col('t')
    kmerDf = (
        seqTotals
        .join(uniqKmers, on=SEQID_COL_NAME, how='left')
        .select(f.alias('uniqueKmers'))
    )

    ## Calculate the average number of unique Kmers across all sequences
    fuk = kmerDf.groupby().mean().collect()[0][0]
    return fuk

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
