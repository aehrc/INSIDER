#!/bin/python

#------------------- Description & Notes --------------------#

'''
Spark doesn't really have a problem in reading lots of large files
since they're essentially read in parts. The problem is moving it
from Spark into memory for native python processing.
'''

#------------------- Dependencies ---------------------------#

# Standard library imports
import functools
import shutil
from pathlib import Path

# External imports
import pyspark.sql.functions as sparkF
import pyspark.sql.types as sparkT
from pyspark.sql import SparkSession

# Internal imports
from .. import kmer
from .common import createDirIfNone, getTempDir

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def read(*dirs):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    filepaths = [str(f) for d in dirs for f in Path(d).glob("*.parquet")]
    kmerDfs   = [_readFile(ss, f) for f in filepaths]

    ## Construct a single Spark DataFrame
    kmerDf    = functools.reduce(lambda x,y: x.union(y), kmerDfs)
    return kmerDf

def write(filepath, kmerDf):
    ## Create DIR if it does not exist
    createDirIfNone(filepath)

    try:
        tmpDir = getTempDir(filepath)

        ## Write the table to file in PARQUET format
        ## We can change the format in the future if it doesn't work out...
        kmerDf.write.parquet(tmpDir, mode='overwrite', compression="snappy")

        ## Cannot write MapType columns in CSV format
        ## So this currently raises an exception
        # kmerDf.write.csv(tmpDir, mode='overwrite')

        ## Remove the hidden files, we don't need them
        [p.unlink() for p in Path(tmpDir).glob("*.crc")]

        ## Move the result files to the filepath
        for source in Path(tmpDir).glob("*.parquet"):
            dest = Path(filepath, source.name)
            source.replace(dest)

    except:
        raise NotImplementedError("Not Implemented!")
        sys.exit(2)

    finally:
        ## Cleanup the temp directory
        shutil.rmtree(tmpDir, ignore_errors=True)

#------------------- Private Classes & Functions ------------#

def _readFile(ss, filepath):
    kmerDf = ss.read.parquet(filepath)
    if (not _isKmerFrequencyFile(kmerDf)):
        raise NotImplementedError("Unknown kmer frequency format")

    ## Must be done before repartitioning
    kmerDf = kmerDf.withColumn(kmer.FILE_COL_NAME,
        _getDirname(sparkF.input_file_name()))

    ## **********
    ## *** Seems like Spark doesn't automatically use up all of the partitions
    ## *** available upon reading. Therefore, we'll repartition the dataframe
    ## *** to maximise the number of partitions available. This should
    ## *** (hopefully) speed up things a bit.
    ## ***
    ## *** Ideally, we should be adjusting the number of partitions based on
    ## *** the total number of sequences. However, if have to repartition
    ## *** from the start, then it's probably easier (and faster) to adjust
    ## *** Spark's defaultParallelism instead of calculating on the fly
    ## **********
    kmerDf = kmerDf.repartition(ss.sparkContext.defaultParallelism,
        kmer.SEQID_COL_NAME)

    ## **********
    ## *** We cannot directly convert the Spark dataframe into a
    ## *** Pandas dataframes because:
    ## ***  * The Kmer column is nested (which Pandas cannot handle)
    ## ***  * Theres a limit on the amount that we can load (i.e., convert
    ## ***    SparkDF into PandasDF with PyArrow) into memory.
    ## *** Therefore, we'll need to do a bit of processing first.
    ## **********
    kmerDf = kmer.transform.agg.dictToCounts(kmerDf)
    return kmerDf

def _isKmerFrequencyFile(kmerDf):
    eCols = {kmer.SEQID_COL_NAME, kmer.TOTAL_COL_NAME, kmer.KMER_COL_NAME}
    oCols = set(kmerDf.schema.names)
    if (not eCols.issubset(oCols)):
        return False

    return True

@sparkF.udf(returnType=sparkT.StringType())
def _getDirname(filepath):
    return Path(filepath).parent.name

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
