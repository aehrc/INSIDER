#!/bin/python

#------------------- Description & Notes --------------------#

'''
We can reduce the number of K-mers by:
* Removing correlated columns
* Removing constant columns

We can reduce the number of samples by:
* Removing duplicate rows (For relatively short K-mers, i.e., K < 5)
  This is actually good for highly related sequences. However, for distantly
  related sequences, this won't actually do much because the counts become
  more unique as the kmer length increases.
'''

#------------------- Dependencies ---------------------------#

# Standard library imports
import hashlib
from operator import add
import zlib

# External imports
import numpy as np
import pandas as pd
import pyspark.sql.functions as sparkF
from pyspark.sql import SparkSession
from pyspark.ml.feature import VectorAssembler
from pyspark.ml.stat import Correlation

# Internal imports
from .convert import kmersToComplexity
from ..constants import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def removeConstantCounts(df):
    ## If our counts are formatted as a Pandas DataFrame
    ## Then we remove constant Kmer columns
    if (isinstance(df, pd.DataFrame)):
        ## Find columns containing the same value across all rows
        ## (i.e., their standard deviation == 0)
        kmerCount    = df
        constantCols = kmerCount.std()[kmerCount.std() == 0]
        kmerCount    = kmerCount.drop(constantCols.index, axis=1)
        return kmerCount

    ## If our counts are formatted as a Spark DataFrame
    ## Then we remove constant Kmer rows
    else:
        ss = SparkSession.getActiveSession()
        if (ss is None):
            raise EnvironmentError("Must have an active Spark session")

        ## Find Kmer containing the same value across all rows
        ## (i.e., their standard deviation == 0)
        kmerSdf      = df
        constantCols = kmerSdf.groupby(KMER_COL_NAME) \
                              .agg(sparkF.stddev(COUNT_COL_NAME)) \
                              .filter(sparkF.col('stddev_samp(count)') == 0) \
                              .select(KMER_COL_NAME)

        kmerSdf = kmerSdf.join(constantCols,
            on=KMER_COL_NAME, how="left_anti")
        return kmerSdf

def removeCorrelatedCounts(df, corr=0.90):
    ## If our counts are formatted as a Pandas DataFrame
    ## Then we remove correlated Kmer columns
    if (isinstance(df, pd.DataFrame)):
        ## Find column pairs that are highly correlated
        kmerCount = df
        kmerCorr  = kmerCount.corr(method='pearson').abs()
        kmerCorr  = kmerCorr.where(np.triu(np.ones(kmerCorr.shape), k=1).astype(np.bool))

        ## We don't really care which one is removed
        corrCols  = [c for c in kmerCorr.columns if any(kmerCorr[c] > corr)]
        kmerCount = kmerCount.drop(corrCols, axis=1)
        return kmerCount

    ## If our counts are formatted as a Spark DataFrame
    ## Then we remove correlated Kmer rows
    else:
        ss = SparkSession.getActiveSession()
        if (ss is None):
            raise EnvironmentError("Must have an active Spark session")

        raise NotImplementedError("Not Implemented properly.")

        # df = ss.createDataFrame(kmerCount)

        # ## Convert kmerCount into Spark Vectors and
        # ## calculate the correlation matrix
        # ## Spark can't seem to handle tables containing > 65535 columns
        # ## Therefore, the maximum Kmer length we can handle (natively) is 8
        # ## (i.e., 8-mers)
        # va = VectorAssembler(inputCols=df.columns, outputCol=KMER_COL_NAME)
        # v  = va.transform(df).select(KMER_COL_NAME)
        # m  = Correlation.corr(v, KMER_COL_NAME, method='pearson')

        # ## Find idxs > corr
        # ## Spark will probably consider this a task of 'very large size'
        # ## because theres a only a single matrix. So we may encounter memory
        # ## issues at longer Kmers unless the work can be distributed
        # m  = m.rdd.map(lambda x: x[0].toArray()) \
        #           .map(np.absolute).map(lambda x: np.argwhere((x > corr))) \
        #           .flatMap(lambda x: x.tolist())

        # ## Ensure that we only get the upper triangle
        # m  = m.repartition(ss.sparkContext.defaultParallelism) \
        #       .filter(lambda x: (x[0] < x[1])).map(lambda x: x[1]) \
        #       .distinct()

        # ## Find columns that are highly correlated
        # colIdxs  = m.collect()
        # colIdxs  = sorted(list(colIdxs))
        # colNames = kmerCount.columns[colIdxs]

        # ## Remove correlated columns
        # ## We don't really care which column/s is removed, but we're
        # ## removing the first ones we encounter
        # kmerCount = kmerCount.drop(colNames, axis=1)
        # return kmerCount

def removeRepetitiveKmers(df):
    ## If our counts are formatted as a Pandas DataFrame
    ## Then we remove repetitive Kmer columns
    if (isinstance(df, pd.DataFrame)):
        ## Find repetitive Kmers; these are considered 'low complexity'.
        ## Complexity is based on the size (in bytes) of the Kmer before
        ## and after compression.
        kmerCount    = df

        f = lambda x: (len(zlib.compress(x.encode())) - len(x.encode()))
        complexity  = {c: f(c) for c in kmerCount.columns}

        threshold   = min(set(complexity.values()))
        complexCols = [k for k, v in complexity.items() if v == threshold]
        kmerCount   = kmerCount.drop(columns=complexCols)
        return kmerCount

    ## If our counts are formatted as a Spark DataFrame
    ## Then we remove repetitive Kmer rows
    else:
        ss = SparkSession.getActiveSession()
        if (ss is None):
            raise EnvironmentError("Must have an active Spark session")

        ## Find repetitive Kmers; these are considered 'low complexity'.
        ## Complexity is based on the size (in bytes) of the Kmer before
        ## and after compression. See Sims et al. (2009).
        kmerSdf     = df

        f = kmersToComplexity(KMER_COL_NAME)
        complexCols = kmerSdf.select(KMER_COL_NAME).distinct() \
                             .withColumn('complexity', f)

        threshold   = complexCols.select(sparkF.min('complexity')).collect()[0][0]
        complexCols = complexCols.filter(sparkF.col('complexity') == threshold) \
                                 .select(KMER_COL_NAME)

        kmerSdf = kmerSdf.join(complexCols, on=KMER_COL_NAME, how="left_anti")
        return kmerSdf

def removeShortSequences(kmerSdf, n=2):
    ## Find relatively short sequences (i.e., len(Seq) < threshold)
    kmerLength = len(kmerSdf.select(KMER_COL_NAME).first()[0])
    minLen     = (4 ** kmerLength) * n
    kmerSdf    = kmerSdf.filter(kmerSdf.total > minLen)
    return kmerSdf

#------------------- Protected Classes & Functions ------------#

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
