#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
from pyspark.sql import SparkSession
import pyspark.sql.functions as sparkF

# Internal imports

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getClusterSizes(kmerDf, cIdDf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Count the number of sequences in each cluster
    csDf = kmerDf.join(sparkF.broadcast(cIdDf), 'id', 'left') \
                   .select('id', 'ClusterId').distinct() \
                   .groupby('ClusterId').count().toPandas()
    csDf = csDf.rename(columns={'count':'NumSeqs'})
    return csDf

def getClusterMeanCounts(kmerDf, cIdDf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Count the number of Kmers in each cluster
    cmcDf = kmerDf.join(sparkF.broadcast(cIdDf), 'id', 'left') \
                     .select('id', 'ClusterId', 'total').distinct() \
                     .groupby('ClusterId').mean('total').toPandas()
    cmcDf = cmcDf.rename(columns={'avg(total)':'AvgNumKmers'})
    return cmcDf

#------------------- Protected Classes & Functions ------------#

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
