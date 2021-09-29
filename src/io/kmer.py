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
from pathlib import Path

# External imports
import pyspark.sql.functions as F
import pyspark.sql.types as T
from pyspark.sql import SparkSession

# Internal imports
from .. import kmer

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
    _printPartitioningInfo(kmerDf)
    return kmerDf

def write(filepath, kmerDf):
    ## Write the table to file in PARQUET format
    ## We can change the format in the future if it doesn't work out...
    kmerDf.write.parquet(filepath, mode='overwrite', compression="snappy")

    # ## Cannot write MapType columns in CSV format
    # ## So this currently raises an exception
    # # kmerDf.write.csv(tmpDir, mode='overwrite')

    ## Spark also saves some hidden files which we don't need.
    ## These were originally removed (for completeness), but there's
    ## no harm in keeping them around.

#------------------- Private Classes & Functions ------------#

def _readFile(ss, filepath):
    kmerDf = ss.read.parquet(filepath)
    if (not _isKmerFrequencyFile(kmerDf)):
        raise NotImplementedError("Unknown kmer frequency format")

    ## Must be done before repartitioning
    f = lambda x: Path(x).parent.name
    f = F.udf(f, T.StringType())
    kmerDf = kmerDf.withColumn(kmer.FILE_COL_NAME, f(F.input_file_name()))

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
    kmerDf = kmerDf.repartitionByRange(ss.sparkContext.defaultParallelism,
        kmer.SEQID_COL_NAME)

    ## **********
    ## *** We cannot directly convert the Spark dataframe into a
    ## *** Pandas dataframes because:
    ## ***  * The Kmer column is nested (which Pandas cannot handle)
    ## ***  * Theres a limit on the amount that we can load (i.e., convert
    ## ***    SparkDF into PandasDF with PyArrow) into memory.
    ## *** Therefore, we'll need to do a bit of processing first.
    ## **********
    kmerDf = _dictToCounts(kmerDf)
    return kmerDf

def _isKmerFrequencyFile(kmerDf):
    eCols = {kmer.SEQID_COL_NAME, kmer.SEQLEN_COL_NAME, kmer.KMER_COL_NAME}
    oCols = set(kmerDf.schema.names)
    if (not eCols.issubset(oCols)):
        return False

    return True

def _dictToCounts(kmerSdf):
    ## (ID, {Kmer -> Percent}, *) => (ID, Kmer, Percent, *)
    gCols = kmerSdf.schema.names
    gCols.remove(kmer.KMER_COL_NAME)

    f        = F.explode(kmer.KMER_COL_NAME)
    colNames = f.alias(kmer.KMER_COL_NAME, kmer.COUNT_COL_NAME)
    kmerSdf  = kmerSdf.select(*gCols, colNames)
    return kmerSdf

def _printPartitioningInfo(kmerDf):
    def _getSeqsPerParts(df):
        ## Avoids shuffling for a more accurate approximation 
        x = df.rdd.glom().first()
        x = map(list, x)
        x = map(lambda x: x[0], x)
        return set(x)

    ## **********
    ## *** N(Partitions) ~= sc.defaultParallelism * N(Files)
    ## ***               ~= N(Partitions | kmerDf)
    ## **********
    nParts = kmerDf.rdd.getNumPartitions()
    seqsPerParts = _getSeqsPerParts(kmerDf)

    ## Print some partitioning info for optimising scalability
    print("Total Num. Partitions\t{}".format(nParts))
    print("Approx Num. Seqs Per Partition\t{}".format(len(seqsPerParts)))

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
