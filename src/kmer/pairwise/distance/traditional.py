#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import numpy as np
import pyspark.sql.functions as F
import pyspark.sql.types as T
from pyspark.sql import SparkSession
from sklearn.metrics import pairwise_distances

# Internal imports
from .common import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getDistances(kmerDfX, kmerDfY=None, **kwargs):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    metric = kwargs.pop('metric', 'euclidean')
    metric = _getMetric(metric)

    ## Find all pairs depending on the DataFrames available
    kmerDf = getPairs(kmerDfX, kmerDfY)

    ## Calculate the distance for each pair
    f = lambda x, y: float(pairwise_distances([x], [y], metric=metric))
    f = F.udf(f, T.FloatType())
    kmerDf     = kmerDf.withColumn('distance', f('l.count', 'r.count'))
    kmerDistDf = kmerDf.select(F.col('l.SeqId').alias('seqId_1'),
        F.col('r.SeqId').alias('seqId_2'), 'distance')
    return kmerDistDf

#------------------- Private Classes & Functions ------------#

def _getMetric(metric):
    if (metric == 'd2'):
        metric = _d2

    return metric

def _d2(srcDf, tgtDf):
    sumNum   = np.sum(srcDf * tgtDf)
    sq_1_sum = np.sum(srcDf ** 2)
    sq_2_sum = np.sum(tgtDf ** 2)

    denom = np.sqrt(sq_1_sum) * np.sqrt(sq_2_sum)
    d2    = 0.5 * (1.0 - (sumNum / denom))
    return d2

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
