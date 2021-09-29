#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import numpy as np
import pandas as pd
from pyspark.sql import SparkSession

# Internal imports
from .common import *
from ... import transform

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getDistances(obsKmerDf, expKmerDf, **kwargs):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    metric = kwargs.pop('metric', 'd2Star')
    metric = _getMetric(metric)

    ## Create a single table containing the Obs and Exp counts
    ## This may takes some time due to the vectorisation
    kmerDf = _createObsExpDf(obsKmerDf, expKmerDf)

    ## Crossjoin the table to find all the pairs we want and
    ## ensure that we only need to process the upper triangle
    kmerDf = getUpperTriangle(kmerDf)

    ## Calculate the distance for each pair
    ## This currently takes a long time to run, but we should be able to
    ## make it faster if we use some udfs
    f = lambda x: (x[0], _toPdf(x[1:3], x[3:5]), x[6], _toPdf(x[7:9], x[9:11]))
    g = lambda x: (x[0], x[2], float(metric(x[1], x[3])))
    kmerRdd    = kmerDf.rdd.map(list).map(f).map(g)
    kmerDistDf = kmerRdd.toDF(['seqId_1', 'seqId_2', 'distance'])
    return kmerDistDf

#------------------- Private Classes & Functions ------------#

def _getMetric(metric):
    if (metric == 'd2Star'):
        metric = _d2Star

    elif (metric == 'd2S'):
        metric = _d2S

    return metric

def _createObsExpDf(obsKmerDf, expKmerDf):
    ## Convert the tables into vectors
    obsKmerDf = transform.vectoriseSdf(obsKmerDf) \
        .withColumnRenamed(KMER_COL_NAME, 'kmer_obs') \
        .withColumnRenamed(COUNT_COL_NAME, 'count_obs') \
        .drop(SEQLEN_COL_NAME, FILE_COL_NAME)
    expKmerDf = transform.vectoriseSdf(expKmerDf) \
        .withColumnRenamed(KMER_COL_NAME, 'kmer_exp') \
        .withColumnRenamed(COUNT_COL_NAME, 'count_exp') \
        .drop(SEQLEN_COL_NAME, FILE_COL_NAME)

    ## Join the tables
    kmerDf = obsKmerDf.join(expKmerDf, on=[SEQID_COL_NAME], how='left')
    return kmerDf

def _toPdf(obs, exp):
    obsDf = pd.DataFrame([obs[1].values], columns=obs[0], index=['Obs'])
    expDf = pd.DataFrame([exp[1].values], columns=exp[0], index=['Exp'])
    df    = pd.concat([obsDf, expDf]).fillna(0)
    return df

def _d2S(srcDf, tgtDf):
    srcDf.index = ["{}_1".format(x) for x in srcDf.index]
    tgtDf.index = ["{}_2".format(x) for x in tgtDf.index]
    df = pd.concat([srcDf, tgtDf]).fillna(0).T

    df['Diff_1'] = df['Obs_1'] - df['Exp_1']
    df['Diff_2'] = df['Obs_2'] - df['Exp_2']
    df['Sq_1']   = df['Diff_1'] * df['Diff_1']
    df['Sq_2']   = df['Diff_2'] * df['Diff_2']
    df['Denom']  = np.sqrt((df['Sq_1'] + df['Sq_2']));

    sumNum   = ((df['Diff_1'] * df['Diff_2']) / df['Denominator']).sum()
    sumDenom = (np.sqrt(df['Sq_1'] / df['Denominator']) * np.sqrt(df['Sq_2'] / df['Denominator'])).sum()
    d2S      = 0.5 * (1.0 - (sumNum / sumDenom))
    return d2S

def _d2Star(srcDf, tgtDf):
    srcDf.index = ["{}_1".format(x) for x in srcDf.index]
    tgtDf.index = ["{}_2".format(x) for x in tgtDf.index]
    df = pd.concat([srcDf, tgtDf]).fillna(0).T

    df['Diff_1']  = df['Obs_1'] - df['Exp_1']
    df['Diff_2']  = df['Obs_2'] - df['Exp_2']
    df['Num']     = (df['Diff_1'] * df['Diff_2']) / np.sqrt(df['Exp_1'] * df['Exp_2'])
    df['SqDiv_1'] = (df['Diff_1'] * df['Diff_1']) /  df['Exp_1']
    df['SqDiv_2'] = (df['Diff_2'] * df['Diff_2']) /  df['Exp_2']

    sumNum   = df['Num'].sum()
    sumDenom = np.sqrt(df['SqDiv_1'].sum()) * np.sqrt(df['SqDiv_2'].sum())
    d2Star   = 0.5 * (1.0 - (sumNum / sumDenom))
    return d2Star

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
