#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import pyspark.sql.functions as sparkF
import pyspark.sql.types as sparkT

# Internal imports
from ..common import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getZScore(kmerStatsDf, colName):
    ## Calculate the mean and standard deviation
    df = kmerStatsDf.select(sparkF.mean(colName), 
        sparkF.stddev(colName)).collect()[0]
    (mean, std) = df

    ## Z-Score calculation
    f = lambda x: (x - mean) / std
    f = sparkF.pandas_udf(f, sparkT.FloatType())
    kmerStatsDf = kmerStatsDf.withColumn('ZScore', f(sparkF.col(colName)))
    return kmerStatsDf

#------------------- Protected Classes & Functions ------------#

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
