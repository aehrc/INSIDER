#!/bin/python

#------------------- Description & Notes --------------------#


#------------------- Dependencies ---------------------------#

# Standard library imports
import math
import hashlib

# External imports
import numpy as np
import pandas as pd
import pyspark.sql.functions as F
from pyspark.sql import SparkSession

# Internal imports
from ..common import *
from .. import transform
from ... import ml

#------------------- Constants ------------------------------#

## Maximum number of Seqs per Partition
CHUNK_SIZE = 1000

#------------------- Public Classes & Functions -------------#

def getLabels(kmerDf, params):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Broadcast the parameters
    sc = ss.sparkContext
    paramsBc = sc.broadcast(params)
    print(paramsBc.value)

    ## Iteratively cluster sequences until we've reached our condition:
    ## That all the data can fit into a single partition. This approach
    ## is similar to heirarchical clustering.
    ##
    ## We do this due to some problems associated with larger data sets
    ## (especially the TSNE). However, we may need to do something
    ## about parameter tuning for each iteration.
    clusterDf = _getClusterDf(kmerDf, paramsBc)
    clusterDf = clusterDf.withColumnRenamed('childId', CID_COL_NAME)
    return clusterDf

def getClusterSizes(kmerDf, cIdDf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Count the number of sequences in each cluster
    csDf = (
        kmerDf
        .join(F.broadcast(cIdDf), SEQID_COL_NAME, 'left')
        .select(SEQID_COL_NAME, CID_COL_NAME).distinct()
        .groupby(CID_COL_NAME).count().toPandas()
    )
    csDf = csDf.rename(columns={'count':'NumSeqs'})
    return csDf

def getClusterMeanCounts(kmerDf, cIdDf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    numKmersDf = (
        kmerDf
        .groupby(KMERDF_SEQINFO_COL_NAMES)
        .agg(F.sum(COUNT_COL_NAME).alias('numKmers'))
    )

    ## Count the number of Kmers in each cluster
    cmcDf = (
        numKmersDf
        .join(F.broadcast(cIdDf), SEQID_COL_NAME, 'left')
        .select(SEQID_COL_NAME, CID_COL_NAME, 'numKmers').distinct()
        .groupby(CID_COL_NAME).mean('numKmers').toPandas()
    )
    cmcDf = cmcDf.rename(columns={'avg(numKmers)':'AvgNumKmers'})
    return cmcDf

#------------------- Private Classes & Functions ------------#

def _getClusterDf(kmerDf, paramsBc):
    ## Repartition the DataFrame so that we can
    ## reduce and cluster more efficiently
    cKmerDf = _repartition(kmerDf)
    nParts  = cKmerDf.rdd.getNumPartitions()

    ## Identify clusters within each partition
    ## (ID, ClusterId)
    parentClusterDf = _getClusters(cKmerDf, paramsBc)

    ## If we can't find anymore clusters, then we're done
    ## New termination condition
    cCount = parentClusterDf.select(CID_COL_NAME).distinct().count()
    sCount = cKmerDf.select(SEQID_COL_NAME).distinct().count()
    if (cCount == sCount):
        return parentClusterDf

    # ## If we only have one partition, then we're done
    # ## Old termination condition
    # if (nParts == 1):
    #     return parentClusterDf

    ## If we have more than one partition, then cluster the counts
    ## and rerun the clustering
    else:
        ## Calculate the average Kmer frequencies for each cluster
        cKmerDf = transform.agg.clusterCountsByColumn(kmerDf, parentClusterDf)

        ## Identify clusters across each partition
        ## (ID, ..., ClusterId)
        childClusterDf = _getClusterDf(cKmerDf, paramsBc)

        ## Rename the columns before we join the tables
        ## (ID, ..., ClusterId) -> (ClusterId, ..., childId)
        childClusterDf = childClusterDf.withColumnRenamed(CID_COL_NAME, 'childId') \
                                       .withColumnRenamed(SEQID_COL_NAME, CID_COL_NAME)

        ## Join the parent and the child tables
        ## Parent = (ID, ClusterId); & Child = (ClusterId, childId)
        ##     newParent = (ID, ClusterId, childId)
        ##               = (ID, parentId, childId)

        ## Parent = (ID, ClusterId); Child = (ClusterId, parentId, childId)
        ##     newParent = (ID, ClusterId, parentId, childId)
        ##               = (ID, parentId, parentId, childId)

        ## Parent = (ID, ClusterId); Child = (ClusterId, parentId, parentId, childI)
        ##     newParent = (ID, ClusterId, parentId, parentId, childId)
        ##               = (ID, parentID, parentId, parentId, childId)
        parentClusterDf = parentClusterDf.join(childClusterDf, CID_COL_NAME, 'left') \
                                         .withColumnRenamed(CID_COL_NAME, 'parentId')
        return parentClusterDf

def _repartition(kmerDf):
    ## Repartition the DataFrame so that we can reduce and cluster more
    ## efficiently later
    nSeqs  = kmerDf.select(SEQID_COL_NAME).distinct().count()
    nParts = math.ceil(nSeqs / CHUNK_SIZE)
    kmerDf = kmerDf.repartition(nParts, SEQID_COL_NAME)
    print("REPARTITION\t{}_{}".format(nSeqs, nParts))
    return kmerDf

def _getClusters(kmerDf, paramsBc):
    intraRdd = _clusterPartitions(kmerDf, paramsBc)
    idxs     = intraRdd.zipWithIndex().values()

    ## Relabel the clusters so that they're different across each partition
    f = lambda x: x.to_numpy().tolist()
    g = lambda x: (x[1][0], [x[1][1], x[0]])
    h = lambda x: hashlib.sha1(repr(tuple(x)).encode()).hexdigest()
    intraDf = idxs.zip(intraRdd).flatMapValues(f).map(g).mapValues(h) \
                  .toDF([SEQID_COL_NAME, CID_COL_NAME])

    ## Persist the results to avoid inconsistencies due to multiple runs
    intraDf.persist()
    return intraDf

def _clusterPartitions(kmerDf, paramBc):
    ## Convert the Spark DataFrame into a Spark RDD of (K, V) pairs
    kmerPdfRdd = transform.toPdfRdd(kmerDf, randomise=True)
    kmerPdfRdd = kmerPdfRdd.map(transform.rotatePdf).map(transform.splitPdf)
    f = lambda x: len(x[0])
    pSizes = kmerPdfRdd.map(f).collect()
    print("PARTITIONS\t{}".format(str(pSizes)))

    ## There are several parameters that need to be tuned for this to work
    ## * Meta clustering
    ## * Consensus clustering
    rKwargs   = paramBc.value.get('rkwargs')
    cKwargs   = paramBc.value.get('ckwargs')
    conKwargs = paramBc.value.get('conkwargs')
    # f = lambda x: _clusterPdf(x, rKwargs, cKwargs, n=1)
    # f = lambda x: _clusterPdf(x, rKwargs, cKwargs, n=10)
    f = lambda x: _clusterPdf(x, rKwargs, cKwargs, n=100)
    g = lambda x: _getConsensusClusters(x, **conKwargs)
    pClusters = kmerPdfRdd.map(f).map(g)
    return pClusters

def _clusterPdf(kmerPdf, reduce_params, cluster_params, n=1):
    kmerId    = kmerPdf[0]
    kmerCount = kmerPdf[1]

    ## To reduce sequence length biases, convert the counts to probabilities
    kmerCount = kmerCount.divide(kmerCount.sum(axis=1), axis=0)

    ## Assign cluster labels to each sequence
    mClusters = [_reduce_n_cluster(kmerCount, reduce_params, cluster_params)
                 for i in range(0, n)]
    mClusters = [m.iloc[:, -1] for m in mClusters]
    mClusters = pd.concat(mClusters, axis=1)
    mClusters.index = kmerId[SEQID_COL_NAME]
    return mClusters

def _reduce_n_cluster(kmerCount, reduce_params, cluster_params):
    # kmerComps  = ml.feature.sklearnReduce(kmerCount, 'PCA',
    #     n_components=2, **reduce_params)
    kmerComps  = ml.feature.sklearnReduce(kmerCount, 'TSNE',
        n_components=2, n_jobs=-1, **reduce_params)
    kmerComps  = ml.feature.scale(kmerComps)

    ## DBSCAN seems to work best with TSNE
    ## Suggested parameters:
    ## algo='DBSCAN', eps=0.05, min_samples=1
    labels     = ml.cluster.assign(kmerComps, 'DBSCAN',
        min_samples=1, **cluster_params)
    compLabels = pd.concat([kmerComps, labels], axis=1)
    return compLabels

def _getConsensusClusters(mClusters, **kwargs):
    ## Suggested parameters:
    ## eps=0.5, min_samples=1
    cClusters = ml.cluster.assignConsensus(mClusters, **kwargs)
    cClusters.index = mClusters.index
    cClusters = cClusters.reset_index()
    return cClusters

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
