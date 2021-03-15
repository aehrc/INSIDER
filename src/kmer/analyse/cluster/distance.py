#!/bin/python

#------------------- Description & Notes --------------------#

'''
Let's look at the distribution of distances amongst Cas9 reads,
amongst non-Cas9 reads and between Cas9 and non-Cas9 reads. This should help
us in determining whether there is an obvious difference between Cas9 reads
and non-Cas9 reads, and therefore any potential in distinguishing between Cas9
and non-Cas9 reads.

I'd expect the distribution amongst Cas9 reads and amongst non-Cas9 reads to be
relatively similar (i.e., skewed towards 0). Meanwhile, I'd expect the
distribution between Cas9 and non-Cas9 reads to be relatively
different (i.e., either bimodal / more evenly distributed).

However, inspection of the plots revealed that the distribution of distances
in all 3 cases were essentially the same. Overall, normally distributed
with a slight skew towards 1. However, there might be a VERY SLIGHT shift in
distribution between the Cas9 and non-Cas9 reads (which may provide some hope).
Why are the distributions essentially the same?
* It could be because the reads actually come from the same
  sequence and therefore the same distribution.

This suggests that detecting any sort of differences between Cas9 reads
and non-Cas9 reads is REALLY difficult (in Yeast, and with 250bp reads).
'''

'''
Compared to DBSCAN and OPTICS (Density-based clustering), I feel that
hierarchical clustering makes a bit more sense given that the reads come
from the same sequence / genome. To better understand the differences,
let's compare the results of the three clustering algorithms.

After some parameter optimisation, DBSCAN returned clusters that seemed
meaningful. Inspection showed that some (2-3) clusters contained
only Cas9 reads, which is great. The size of these clusters ranged
between 10-40 reads, which raises the question of whether they will
provide enough statistical power later on?

Like DBSCAN, OPTICS returned some clusters, but each cluster contained
fewer reads than with DBSCAN. Because of this, we may not have enough
statistical power later on.

Finally, hierarchical clustering assigns every point to a cluster.
Inspection of hierarchical clusters revealed several clusters containing
Cas9 reads, suggesting that Cas9 reads do not exclusively cluster
together. However, there were several clusters containing only non-Cas9 reads.
'''

#------------------- Dependencies ---------------------------#

# Standard library imports
import sys

# External imports
import numpy as np
import pandas as pd
import pyspark.sql.functions as sparkF
import pyspark.sql.types as sparkT
from pyspark.sql import SparkSession

# Internal imports
from ...constants import *
from ... import ops
from ... import transform
from ... import distance
from .... import ml

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getLabels(oKmerDf, eKmerDf, algo, **kwargs):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Calculate the distance between Kmer frequencies
    ## This requires LOTS of resources since we're essentially doing
    ## an all vs all comparison. However, in order to use the
    ## distance matrix for clustering, we MUST load the distance matrix
    ## into memory for sklearn. The alternative is to use a
    ## non-sklearn clustering algorithm.
    print("Calculating d2S distance")
    kmerDist = _getDistanceMatrix(oKmerDf, eKmerDf)
    print(kmerDist.head())
    kwargs['metric'] = 'precomputed'

    ## Identify clusters based on the distance matrix
    ## Ideally, we should only rely on a single optimised clustering algorithm
    print("Clustering distances")
    print(kwargs)
    labels = _getClusters(kmerComp, algo, **kwargs)
    print(labels.head())

    # ## Fix up the table
    # args = ['{}:{}'.format(k, v) for k, v in kwargs.items()]
    # args = '_'.join(args)
    # colName = '{}_{}_{}'.format(labels.columns.tolist()[0], algo, args)
    # labels.columns = [colName]
    # labels.index   = kmerDist.index
    # labels = labels.reset_index().rename(columns={'index':'id'})
    return labels

#------------------- Private Classes & Functions ------------#

def _getDistanceMatrix(oKmerDf, eKmerDf):
    ## Normalise counts for the d2S distance calculation
    f = transform.counts.countsToCentralisedCounts
    oNormKmerDf = oKmerDf.groupby(oKmerDf.id) \
                         .cogroup(eKmerDf.groupby(eKmerDf.id)) \
                         .applyInPandas(f, schema=eKmerDf.schema)
    # oNormKmerDf = sampleRows(oNormKmerDf)

    ## This will load the distance matrix into memory
    # kmerDist = distance.calculate(oNormKmerDf, metric='d2S', n_jobs=-1)
    kmerDist = distance.calculate(oNormKmerDf, metric='euclidean', n_jobs=-1)
    np.fill_diagonal(kmerDist.to_numpy(), 0)
    # kmerDist = distance.scale(kmerDist)
    return kmerDist

def _getClusters(kmerDf, algo, **kwargs):
    labels = ml.cluster.assign(kmerDf, algo, **kwargs)
    return labels

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
