#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import math

# External imports
import pyspark.sql.functions as F
from pyspark.sql import SparkSession

# Internal imports
from ..common import *
from .. import transform

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getPairs(kmerSdfX, kmerSdfY=None):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    if (kmerSdfY is None):
        ## Convert the counts into sorted lists
        kmerSdfX = transform.agg.countsToSortedList(kmerSdfX) \
            .drop(SEQLEN_COL_NAME, FILE_COL_NAME)

        nParts   = kmerSdfX.rdd.getNumPartitions()
        kmerSdfX = kmerSdfX.coalesce(int(math.sqrt(nParts)))
        kmerSdf  = getUpperTriangle(kmerSdfX)

    else:
        ## Convert the counts into sorted lists
        kmerSdfX = transform.agg.countsToSortedList(kmerSdfX) \
            .drop(SEQLEN_COL_NAME, FILE_COL_NAME)
        kmerSdfY = transform.agg.countsToSortedList(kmerSdfY) \
            .drop(SEQLEN_COL_NAME, FILE_COL_NAME)

        nParts   = kmerSdfX.rdd.getNumPartitions()
        kmerSdfX = kmerSdfX.coalesce(int(math.sqrt(nParts)))
        kmerSdf  = kmerSdfX.alias('l').crossJoin(kmerSdfY.alias('r'))
    return kmerSdf

def getUpperTriangle(kmerSdf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Insert a column containing the row number so that we only need
    ## to process the upper triangle and crossjoin the table to find all
    ## the pairs we want
    kmerSdf = kmerSdf.withColumn('rId', F.monotonically_increasing_id())
    cond = F.col('l.rId') <= F.col('r.rId')
    kmerSdf = kmerSdf.alias('l').join(kmerSdf.alias('r'), on=cond, how='outer')
    return kmerSdf

#------------------- Protected Classes & Functions ------------#

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
