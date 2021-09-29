#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import zlib

# External imports
import pandas as pd
import pyspark.sql.functions as sparkF
import pyspark.sql.types as sparkT
from pyspark.sql import SparkSession

# Internal imports
from .common import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def countsToProbabilities(kmerSdf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    kmerSdf = kmerSdf.withColumn(COUNT_COL_NAME,
        sparkF.col(COUNT_COL_NAME).cast(sparkT.FloatType()))
    kmerSdf = kmerSdf.repartition(kmerSdf.rdd.getNumPartitions(), SEQID_COL_NAME)
    kmerSdf = kmerSdf.groupby(SEQID_COL_NAME) \
        .applyInPandas(_countsToProbabilities, schema=kmerSdf.schema)
    return kmerSdf

def countsToNormalised(kmerSdf):
    ## Normalise counts according to Wang et al. 2005
    def _f(key, df, t):
        df   = _countsToProbabilities(df)
        df[COUNT_COL_NAME] = df[COUNT_COL_NAME] * t
        return df

    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Calculate the total number of possible Kmers
    kmer = kmerSdf.select(KMER_COL_NAME).first()[0]
    t    = getExpectedTotal(kmer)

    kmerSdf = kmerSdf.withColumn(COUNT_COL_NAME,
        sparkF.col(COUNT_COL_NAME).cast(sparkT.FloatType()))
    kmerSdf = kmerSdf.repartition(kmerSdf.rdd.getNumPartitions(), SEQID_COL_NAME)
    kmerSdf = kmerSdf.groupby(SEQID_COL_NAME) \
        .applyInPandas(lambda x, y: _f(x, y, t), schema=kmerSdf.schema)
    return kmerSdf

@sparkF.pandas_udf(returnType=sparkT.FloatType())
def countsToPercentagesUdf(counts: pd.Series, total: pd.Series) -> pd.Series:
    pcts = countsToProbabilitiesUdf(counts, total)
    pcts = probabilitiesToPercentagesUdf(pcts)
    return pcts

@sparkF.pandas_udf(returnType=sparkT.FloatType())
def countsToProbabilitiesUdf(counts: pd.Series, total: pd.Series) -> pd.Series:
    return ((counts / total))

@sparkF.pandas_udf(returnType=sparkT.FloatType())
def probabilitiesToPercentagesUdf(counts: pd.Series) -> pd.Series:
    return ((counts * 100))

@sparkF.pandas_udf(returnType=sparkT.IntegerType())
def percentagesToCountsUdf(pcts: pd.Series, total: pd.Series) -> pd.Series:
    counts = percentagesToProbabilitiesUdf(pcts)
    counts = probabilitiesToCountsUdf(counts, total)
    return counts

@sparkF.pandas_udf(returnType=sparkT.FloatType())
def percentagesToProbabilitiesUdf(counts: pd.Series) -> pd.Series:
    return ((counts / 100))

@sparkF.pandas_udf(returnType=sparkT.IntegerType())
def probabilitiesToCountsUdf(counts: pd.Series, total: pd.Series) -> pd.Series:
    return ((counts * total))

@sparkF.pandas_udf(returnType=sparkT.FloatType())
def kmersToComplexityUdf(kmer: pd.Series) -> pd.Series:
    f = lambda x: (len(zlib.compress(x.encode())) - len(x.encode()))
    kmer = kmer.apply(f)
    return kmer

#------------------- Protected Classes & Functions ------------#

#------------------- Private Classes & Functions ------------#

def _countsToProbabilities(kmerPdf):
    total = kmerPdf[COUNT_COL_NAME].sum()
    kmerPdf[COUNT_COL_NAME] = kmerPdf[COUNT_COL_NAME] / total
    return kmerPdf

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
