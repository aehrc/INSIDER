#!/bin/python

#------------------- Description & Notes --------------------#



#------------------- Dependencies ---------------------------#

# Standard library imports
import random

# External imports
import pandas as pd
import pyspark.sql.functions as sparkF
import pyspark.sql.types as sparkT
from pyspark.sql import SparkSession
from pyspark.ml.linalg import Vectors, VectorUDT, _convert_to_vector
import scipy.sparse

# Internal imports
from ..constants import *
from ..common import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def toPdfRdd(kmerSdf, randomise=False):
    cColType = kmerSdf.select(COUNT_COL_NAME).dtypes[0]
    if (cColType[1] == 'vector'):
        idRdd    = kmerSdf.select(*KMERDF_COL_NAMES).rdd
        countRdd = kmerSdf.select(COUNT_COL_NAME).rdd
        kmerRdd  = kmerSdf.select(KMER_COL_NAME).distinct().rdd
        kmers    = kmerRdd.flatMap(lambda x: x).collect()[0]

        f = lambda x: [pd.DataFrame(list(x), columns=KMERDF_COL_NAMES)]
        g = lambda x: len(x) != 0
        idRdd = idRdd.mapPartitions(f).filter(g)

        f = lambda x: [pd.DataFrame(list(x), columns=kmers)]
        countRdd   = countRdd.map(lambda x: x[0].toArray().tolist()) \
            .mapPartitions(f).filter(g)

        ## Convert the Spark Dataframe into Spark RDD of (K, V) pairs
        kmerPdfRdd = idRdd.zip(countRdd)

    else:
        if (randomise):
            cols    = kmerSdf.columns
            nParts  = kmerSdf.rdd.getNumPartitions()

            ## Attach a column of random numbers
            f       = sparkF.lit(random.randint(0, nParts))
            ids     = kmerSdf.select(SEQID_COL_NAME).distinct() \
                             .select(SEQID_COL_NAME, f.alias('rand'))
            kmerSdf = kmerSdf.join(ids, SEQID_COL_NAME, 'left')

            ## Repartition the data to randomise the data
            kmerSdf = kmerSdf.repartition(nParts, SEQID_COL_NAME, 'rand')
            kmerSdf = kmerSdf.select(cols)

        cols = kmerSdf.schema.names
        f    = lambda x: [pd.DataFrame(list(x), columns=cols)]
        g    = lambda x: len(x) != 0
        kmerPdfRdd = kmerSdf.rdd.mapPartitions(f).filter(g)

    return kmerPdfRdd

def rotatePdf(kmerPdf):
    pivotColNames = set([KMER_COL_NAME, COUNT_COL_NAME])
    idxColNames   = set(kmerPdf.columns).difference(pivotColNames)

    ## ID columns for each group must contain the same values
    ## for the pivot to work properly
    kmerPdf = pd.pivot_table(kmerPdf, index=idxColNames,
        columns=KMER_COL_NAME, values=COUNT_COL_NAME, fill_value=0)
    kmerPdf = kmerPdf.reset_index()
    kmerPdf.columns.name = None
    return kmerPdf

def splitPdf(kmerPdf, colName=FILE_COL_NAME):
    ## Split Pdf by column (placed on the left)
    splitIdx  = kmerPdf.columns.get_loc(colName) + 1
    idCols    = kmerPdf.iloc[:, :splitIdx]
    countCols = kmerPdf.iloc[:, splitIdx:]
    return (idCols, countCols)

def vectoriseSdf(kmerSdf, toSparse=False):
    def _addZeroCounts_n_countsToVector(key, kmerPdf, oKmers):
        kmerPdf = addZeroCounts(key, kmerPdf, oKmers)

        ## Convert the counts into a sorted list based on the Kmers
        kmerPdf = kmerPdf.sort_values(by=[KMER_COL_NAME])
        counts  = kmerPdf.groupby([*KMERDF_COL_NAMES])[COUNT_COL_NAME] \
            .apply(list).reset_index()
        kmerPdf = kmerPdf.groupby([*KMERDF_COL_NAMES])[KMER_COL_NAME] \
            .apply(list).reset_index()
        kmerPdf[COUNT_COL_NAME] = counts[COUNT_COL_NAME]
        return kmerPdf

    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Find Kmers that do not exist in each sequence (i.e., zero count Kmers)
    sc     = ss.sparkContext
    oKmers = getOuterKmers(kmerSdf)
    bc     = sc.broadcast(oKmers)

    ## Create a new schema
    colNames = [*KMERDF_COL_NAMES, KMER_COL_NAME, COUNT_COL_NAME]
    colTypes = [sparkT.StringType(),
        sparkT.LongType(), sparkT.StringType(),
        sparkT.ArrayType(sparkT.StringType()),
        sparkT.ArrayType(sparkT.FloatType())]
    cols     = [sparkT.StructField(c, t) for c, t in zip(colNames, colTypes)]
    schema   = sparkT.StructType(cols)

    ## Normalise Kmer columns
    f       = lambda x, y: _addZeroCounts_n_countsToVector(x, y, oKmers)
    kmerSdf = kmerSdf.groupby(SEQID_COL_NAME) \
        .applyInPandas(f, schema=schema)

    ## Convert the counts into Spark DenseVectors
    f = lambda x: Vectors.dense(x)
    f = sparkF.udf(f, VectorUDT())
    kmerSdf = kmerSdf.withColumn(COUNT_COL_NAME,
        f(COUNT_COL_NAME).alias(COUNT_COL_NAME))

    ## If applicable, convert the counts into Spark SparseVectors
    if (toSparse):
        f = lambda x: _convert_to_vector(scipy.sparse.csc_matrix(x.toArray()).T)
        f = sparkF.udf(f, VectorUDT())
        kmerSdf = kmerSdf.withColumn(COUNT_COL_NAME,
            f(COUNT_COL_NAME).alias(COUNT_COL_NAME))

    return kmerSdf

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
