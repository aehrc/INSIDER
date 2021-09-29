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

# Internal imports

#------------------- Constants ------------------------------#

THIS_FILE = Path(__file__)
THIS_DIR  = Path(THIS_FILE.resolve().parent)
SRC_DIR   = Path(THIS_DIR, '..')
SRC_ZIP   = Path(SRC_DIR, 'modules.zip')

#------------------- Public Classes & Functions -------------#

def getSparkSession(params):
    ## Spark configuration
    ## https://aws.amazon.com/blogs/big-data/best-practices-for-successfully-managing-memory-for-apache-spark-applications-on-amazon-emr/
    conf = pyspark.SparkConf()
    [conf.set(k, v) for k,v in params]

    ## Create SparkSession
    ss   = SparkSession.builder.config(conf=conf).getOrCreate()
    ss.sparkContext.addPyFile(str(SRC_ZIP))
    return ss

#------------------- Private Classes & Functions ------------#

## For improving pickling process
## https://stackoverflow.com/questions/53371112/creating-parquet-petastorm-dataset-through-spark-fails-with-overflow-error-larg
def broadcast_dump(self, value, f):
    pickle.dump(value, f, 4)  # was 2, 4 is first protocol supporting >4GB
    f.close()
    return f.name

pyspark.broadcast.Broadcast.dump = broadcast_dump

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
