#!/bin/python

#------------------- Description & Notes --------------------#

'''
We can reduce the number of K-mers by:
* Clustering reverse complement counts (For relatively short K-mers, i.e. K < 5)

We can reduce the number of samples by:
* Clustering counts of related sequences
'''

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import pandas as pd
import pyspark.sql.functions as F
import pyspark.sql.types as T
from pyspark.sql import SparkSession

# Internal imports
from .convert import countsToProbabilitiesUdf
from .common import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def countsToDict(kmerSdf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Must be done in Spark because Pandas doesn't allow maptype columns
    ## (ID, Kmer, Percent, *)
    ##     => (ID, Kmer -> Percent, *)
    ##     => (ID, {Kmer -> Percent}, *)
    gCols = kmerSdf.schema.names
    gCols.remove(KMER_COL_NAME)
    gCols.remove(COUNT_COL_NAME)

    schema  = F.struct([KMER_COL_NAME, COUNT_COL_NAME])
    f       = F.map_from_entries(F.collect_list(schema))
    kmerSdf = kmerSdf.groupBy(gCols).agg(f.alias(KMER_COL_NAME))
    return kmerSdf

def countsToSortedList(kmerSdf):
    def _createSchema():
        ## Create a new schema
        colTypes = [
            T.StringType(), T.LongType(), T.StringType(),
            T.ArrayType(T.StringType()), T.ArrayType(T.FloatType())
        ]
        cols     = [T.StructField(c, t) for c, t in zip(KMERDF_COL_NAMES, colTypes)]
        schema   = T.StructType(cols)
        return schema

    def _addZeroCounts(kmerSdf):
        ## Collect Kmers in each sequence
        f = F.collect_list(KMER_COL_NAME)
        y = kmerSdf.groupby(KMERDF_SEQINFO_COL_NAMES).agg(f.alias(KMER_COL_NAME))

        ## Collect Kmers across all sequences
        f = F.collect_set(KMER_COL_NAME)
        x = kmerSdf.select(f.alias('oKmers'))

        ## Join the tables
        z = x.crossJoin(y)

        ## Find Kmers in some but not all sequences
        f = F.explode(F.array_except('oKmers', KMER_COL_NAME))
        g = F.lit(0)
        cols = [*KMERDF_SEQINFO_COL_NAMES, f.alias(KMER_COL_NAME), g.alias(COUNT_COL_NAME)]
        z = z.select(cols)

        ## Add them to the original table
        kmerSdf = kmerSdf.union(z)
        return kmerSdf

    def _sortCounts(kmerPdf):
        ## Convert the counts into a sorted list based on the Kmers
        kmerPdf = kmerPdf.sort_values(by=[KMER_COL_NAME])
        counts  = kmerPdf.groupby([*KMERDF_SEQINFO_COL_NAMES])[COUNT_COL_NAME] \
            .apply(list).reset_index()
        kmerPdf = kmerPdf.groupby([*KMERDF_SEQINFO_COL_NAMES])[KMER_COL_NAME] \
            .apply(list).reset_index()
        kmerPdf[COUNT_COL_NAME] = counts[COUNT_COL_NAME]
        return kmerPdf

    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Create schema
    schema = _createSchema()
    nParts = kmerSdf.rdd.getNumPartitions()

    ## Find Kmers that do not exist in each sequence (i.e., zero count Kmers)
    ## and add them to the table
    ## Repartitioning 'should' make groupbys run a bit faster
    kmerSdf = kmerSdf.repartition(nParts, KMERDF_SEQINFO_COL_NAMES)
    kmerSdf = _addZeroCounts(kmerSdf)

    ## Normalise Kmer columns
    kmerSdf = kmerSdf.repartition(nParts, SEQID_COL_NAME)
    kmerSdf = (
        kmerSdf
        .groupby(SEQID_COL_NAME)
        .applyInPandas(_sortCounts, schema=schema)
    )
    return kmerSdf

def sumCounts(kmerSdf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Sum the counts for each Kmer
    newKmerSdf = kmerSdf.groupby(KMER_COL_NAME) \
        .agg(F.sum(COUNT_COL_NAME).alias(COUNT_COL_NAME))

    ## Sum the lengths of each sequence
    seqLen = kmerSdf.select(SEQID_COL_NAME, SEQLEN_COL_NAME).distinct() \
        .groupby().sum(SEQLEN_COL_NAME).collect()[0][0]

    ## Format the rest of the columns and adjust the total counts
    newKmerSdf = newKmerSdf.withColumn(SEQID_COL_NAME, F.lit('Sum'))
    newKmerSdf = newKmerSdf.withColumn(FILE_COL_NAME, F.lit('Sum'))
    newKmerSdf = newKmerSdf.withColumn(SEQLEN_COL_NAME, F.lit(seqLen))
    return newKmerSdf

def averageCounts(kmerSdf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Sum the counts for each K-mer across all sequences and
    ## divide by the number of sequences
    newKmerSdf = sumCounts(kmerSdf)
    numSeqs    = kmerSdf.select(SEQID_COL_NAME).distinct().count()
    f          = countsToProbabilitiesUdf(COUNT_COL_NAME, F.lit(numSeqs))
    newKmerSdf = newKmerSdf.withColumn(COUNT_COL_NAME, f)
    return newKmerSdf

def clusterRevCompCounts(kmerSdf):
    def _f(df):
        ## Create a column containing the reverse complement of each sequence
        df['revComp'] = df[KMER_COL_NAME].apply(getReverseComplement)

        ## Sum up the counts of complementary K-mers, and double the counts
        ## for the total frequency of each K-mer.
        f = lambda x: '-'.join(sorted(x))
        df[KMER_COL_NAME] = df[[KMER_COL_NAME, 'revComp']].apply(f, axis=1)
        df = df.groupby([*KMERDF_SEQINFO_COL_NAMES, KMER_COL_NAME]).sum() * 2
        df = df.reset_index()
        return df

    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    kmerSdf = kmerSdf.groupby(SEQID_COL_NAME) \
        .applyInPandas(_f, schema=kmerSdf.schema)
    return kmerSdf

def clusterCountsByColumn(kmerSdf, cIdDf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    cKmerDf = kmerSdf.join(cIdDf, SEQID_COL_NAME, 'left') \
        .groupby(CID_COL_NAME) \
        .applyInPandas(sumCountsByColumn, schema=kmerSdf.schema)
    return cKmerDf

def sumCountsByColumn(key, kmerPdf):
    ## Sum the counts for each Kmer
    cols       = [KMER_COL_NAME, COUNT_COL_NAME]
    newKmerPdf = kmerPdf[cols].groupby(KMER_COL_NAME).sum()
    newKmerPdf = newKmerPdf.reset_index()

    ## Sum the lengths of each sequence
    cols   = [SEQID_COL_NAME, SEQLEN_COL_NAME]
    seqLen = kmerPdf[cols].drop_duplicates()[SEQLEN_COL_NAME].sum()

    ## Format the rest of the columns and adjust the total counts
    newKmerPdf[SEQID_COL_NAME]  = str(key[0])
    newKmerPdf[FILE_COL_NAME]   = str(key[0])
    newKmerPdf[SEQLEN_COL_NAME] = seqLen
    return newKmerPdf

def averageCountsByColumn(key, kmerPdf):
    ## Sum the counts for each K-mer across sequences
    ## and divide by the number of sequences
    newKmerPdf = sumCountsByColumn(key, kmerPdf)
    numSeqs    = len(kmerPdf[SEQID_COL_NAME].unique())
    newKmerPdf[COUNT_COL_NAME] = newKmerPdf[COUNT_COL_NAME] / numSeqs
    return newKmerPdf

#------------------- Protected Classes & Functions ------------#

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
