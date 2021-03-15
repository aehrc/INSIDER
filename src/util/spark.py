#!/bin/python

#------------------- Description & Notes --------------------#

'''
If we encounter ModuleNotFound, check that the file
paths for the zipped files are relative to the main scripts
'''

#------------------- Dependencies ---------------------------#

# Standard library imports
import pickle
from pathlib import Path

# External imports
import pyspark
from pyspark.sql import SparkSession
from pyspark.sql.dataframe import DataFrame

# Internal imports

#------------------- Constants ------------------------------#

THIS_FILE = Path(__file__)
THIS_DIR  = Path(THIS_FILE.resolve().parent)
SRC_DIR   = Path(THIS_DIR, '..')
SRC_ZIP   = Path(SRC_DIR, 'modules.zip')

#------------------- Public Classes & Functions -------------#

def getSparkSession():

    """
    Description:
        Generate a SparkSession

    Returns:
        sc (SparkContext)
            SparkSession environment
    """

    conf = getSparkConfiguration()
    ss   = SparkSession.builder.config(conf=conf).getOrCreate()
    ss.sparkContext.addPyFile(str(SRC_ZIP))
    return ss

#------------------- Private Classes & Functions ------------#

def getSparkConfiguration():
    conf = pyspark.SparkConf()
    conf.set('spark.driver.memory', '128G')
    conf.set('spark.driver.maxResultSize', '0')
    # conf.set('spark.executor.cores', '5')         ## Didn't seem to have any effect...
    # conf.set('spark.executor.instances', '21')    ## Didn't seem to have any effect...
    conf.set('spark.executor.memory', '64G')
    conf.set('spark.executor.heartbeatInterval', '60s')
    conf.set('spark.local.dir', THIS_DIR)
    conf.set('spark.network.timeout', '600s')
    conf.set('spark.sql.broadcastTimeout', '600s')
    conf.set('spark.sql.execution.arrow.pyspark.enabled', 'true')
    conf.set('spark.sql.shuffle.partitions', '200')
    conf.set('spark.default.parallelism', '8')
    return conf

## For chaining dataframe transformations
## https://adatis.co.uk/pyspark-dataframe-transformation-chaining/
def transform(self, f):
    return f(self)

## For improving pickling process
## https://stackoverflow.com/questions/53371112/creating-parquet-petastorm-dataset-through-spark-fails-with-overflow-error-larg
def broadcast_dump(self, value, f):
    pickle.dump(value, f, 4)  # was 2, 4 is first protocol supporting >4GB
    f.close()

    return f.name

pyspark.broadcast.Broadcast.dump = broadcast_dump

DataFrame.transform = transform

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
