#!/bin/python

#------------------- Description & Notes --------------------#


#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import pandas as pd
from pyspark.sql import SparkSession

# Internal imports
from ..common import *
from .. import transform
from ... import ml

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getLabels(kmerDf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Assume that we can fit everything into a single partition
    ## and convert the Spark DataFrame into a Spark RDD of (K, V) pairs
    kmerDf  = kmerDf.coalesce(1)
    kmerRdd = transform.toPdfRdd(kmerDf)
    kmerRdd = kmerRdd.map(transform.rotatePdf).map(transform.splitPdf)

    ## Dimensionality reduction and outlier detection
    f = lambda x: _findOutliers(x, n=100)
    g = lambda x: _getConsensusOutliers(x)
    outlierRdd = kmerRdd.map(f).map(g)

    f = lambda x: x.to_numpy().tolist()
    outlierDf = outlierRdd.flatMap(f).toDF([SEQID_COL_NAME, 'OutlierPct'])

    ## Persist the results to avoid inconsistent results due to multiple runs
    outlierDf.persist()
    return outlierDf

#------------------- Private Classes & Functions ------------#

def _getConsensusOutliers(mOutliers):
    def _CC(v):
        ## Calculate the proportion of values
        v = v.value_counts()
        v = v / sum(v)

        ## Check that we have counts for the outlier flag
        ## If we don't, insert it into the series with 0 counts
        if (-1 not in v):
            v.loc[-1] = 0

        v = v[[-1]]
        v.index = ['OP']
        return v

    cOutliers = mOutliers.apply(_CC, axis=1)
    cOutliers = cOutliers.reset_index()
    return cOutliers

def _findOutliers(kmerPdf, n=1):
    kmerId    = kmerPdf[0]
    kmerCount = kmerPdf[1]

    ## To reduce sequence length biases, convert the counts to probabilities
    kmerCount = kmerCount.divide(kmerCount.sum(axis=1), axis=0)

    ## Assign outlier labels to each sequence
    mOutliers = [_reduce_n_outlier(kmerCount) for i in range(0, n)]
    mOutliers = [m.iloc[:, -1] for m in mOutliers]
    mOutliers = pd.concat(mOutliers, axis=1)
    mOutliers.index = kmerId[SEQID_COL_NAME]
    return mOutliers

def _reduce_n_outlier(kmerCount):
    kmerComps  = ml.feature.sklearnReduce(kmerCount, 'PCA',
        n_components=0.75)

    labels     = ml.outlier.detect(kmerComps, 'IF', n_jobs=-1)
    compLabels = pd.concat([kmerComps, labels], axis=1)
    return compLabels

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
