#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import pandas as pd
from pyspark.sql import SparkSession
from scipy import stats

# Internal imports
from ...constants import *
from ... import transform

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getSignificantPairs(pvDf, sig=0.01):
    ## Find rows that are below the significance threshold
    cond = (pvDf['p_value'] < sig)
    pvDf['isSignificant'] = cond
    return pvDf

## Get the Zscore of a column
def getZScore(pvDf, colName):
    z = stats.zscore(pvDf[[colName]])
    z = pd.DataFrame(z, columns=['ZScore'], index=pvDf.index)
    zColName = 'ZScore_{}'.format(colName)
    pvDf[zColName] = z
    return pvDf

#------------------- Protected Classes & Functions ------------#

def getCountsPdf(kmerSdf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    kmerSdf    = kmerSdf.coalesce(1)
    kmerRdd    = transform.toPdfRdd(kmerSdf)
    kmerRdd    = kmerRdd.map(transform.rotatePdf).map(transform.splitPdf)
    kmerPdf    = kmerRdd.collect()[0]
    kmerCount  = kmerPdf[1]
    kmerCount.index = kmerPdf[0][SEQID_COL_NAME]
    return kmerCount

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
