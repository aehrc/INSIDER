#!/bin/python

#------------------- Description & Notes --------------------#

'''
D2S and D2Star seem to be the 'standard' metric for measuring
dissimilarity between Kmer frequencies.

In terms of definition they are pretty similar and D2Star reportedly
performs better than D2S. However, D2S appears to be more widely
used because its easier to compute.
'''

'''
Spark allows us to do very large and scalable pairwise comparisons
However, doing anything useful with a super large square distance matrix
will be very problematic.

One thing with Spark is that it can't seem to handle wide data
very well (at least in the case of the DataFrames). This can get a
bit problematic when we want to do pairwise comparisons as the
square matrices become MASSIVE.

Alternatively, we could find out the distances for only the X closest points
since we're predominantly interested in points that are relatively close.
Therefore, we'll end up with a table of idxs instead of distances.
'''

#------------------- Dependencies ---------------------------#

# Standard library imports
import math

# External imports
import numpy as np
import pandas as pd
from pyspark.sql import SparkSession
from sklearn.metrics import pairwise_distances
from sklearn import preprocessing

# Internal imports
from ..constants import *
from .. import ops
from .. import transform
from .transform import matrixToTable

#------------------- Constants ------------------------------#

## Maximum number of Seqs per Partition
CHUNK_SIZE = 100000

## Number of Partitions multiplier for Spark
N_PARTS_MULTI = 6

#------------------- Public Classes & Functions -------------#

def calculate(kmerDf, **kwargs):
    ## Set up Spark to calculate pairwise distances
    ## The idea for this is from:
    ## https://medium.com/@rantav/large-scale-matrix-multiplication-with-pyspark-or-how-to-match-two-large-datasets-of-company-1be4b1b2871e
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    asMatrix  = kwargs.pop('asMatrix', True)
    metric    = kwargs.pop('metric', 'euclidean')

    (bc, rdd) = _setup(ss, kmerDf)
    metric    = _getMetric(metric)
    kmerDist  = _getKmerDist(rdd, bc, metric, **kwargs)
    if (not asMatrix):
        kmerDist = matrixToTable(kmerDist)

    return kmerDist

def findById(dist, id_xs, id_ys=None):
    p    = '|'.join(id_xs)
    cond = ((dist['id_x'].str.contains(p)))

    if (id_ys is not None):
        p       = '|'.join(id_ys)
        idyCond = ((dist['id_y'].str.contains(p)))
        cond    = (cond & idyCond)

    return dist[cond]

def compare(dist1, dist2):
    idColNames   = ['id_x', 'id_y']
    dist         = dist1.merge(dist2, left_on=idColNames, right_on=idColNames)
    distColNames = ['distance_x', 'distance_y']
    corr         = dist['distance_x'].corr(dist['distance_y'])
    return (dist, corr)

def scale(dist, scale=(0, 100)):
    scaler = preprocessing.MinMaxScaler(scale)
    c      = scaler.fit_transform(dist)
    c      = pd.DataFrame(c, columns=dist.columns, index=dist.index)
    return c

#------------------- Private Classes & Functions ------------#

def _setup(ss, kmerDf):
    sc  = ss.sparkContext

    ## Construct the RDD
    rdd = transform.toPdfRdd(kmerDf)
    rdd = rdd.map(transform.rotatePdf).map(transform.splitPdf) \
             .flatMap(ops.partitionByRows)
    rdd = _repartition(rdd, sc.defaultParallelism)
    print(rdd.getNumPartitions())

    ## Construct the broadcast variable
    ## Unfortunately, this will load all the Kmer counts into memory
    kmerCount = rdd.values().reduce(lambda x,y: pd.concat([x, y]))
    kmerCount = kmerCount.fillna(0).to_numpy().astype(np.float64)
    bc  = sc.broadcast(kmerCount)
    return (bc, rdd)

def _repartition(kmerPdfRdd, defParts):
    ## Partition sequences. This should be proportional to
    ## the number of Pdfs
    eParts = math.floor(kmerPdfRdd.count() / CHUNK_SIZE)
    nParts = defParts + (eParts * N_PARTS_MULTI)
    return kmerPdfRdd.repartition(nParts)

def _getKmerDist(rdd, bc, metric, **kwargs):
    ## Calculate distances
    ## This will load the distance matrix into memory
    kmerDists = _getSubKmerDists(rdd, bc, metric, **kwargs)
    kmerDist  = pd.concat(kmerDists)

    ## Rename the columns
    kmerDist.columns      = kmerDist.index
    kmerDist.columns.name = None
    kmerDist.index.name   = None
    return kmerDist

def _getSubKmerDists(rdd, bc, metric, **kwargs):
    ## Evaluate distances between vectors
    ## We encounter imprecise floating points (in particular, the diagonal
    ## is non-zero). I fixed them in the past by typecasting the values, but
    ## for some reason it's not working anymore. It'll be good to fix this up,
    ## but I suspect this won't significantly impact the overall results
    f = lambda x: pairwise_distances(x.to_numpy().astype(np.float64), bc.value,
        metric=metric, **kwargs)
    g = lambda x: pd.DataFrame(x[1], index=x[0][SEQID_COL_NAME])
    kmerDists = rdd.mapValues(f).map(g)

    for kmerDist in kmerDists.toLocalIterator():
        print(kmerDist)
        yield kmerDist

def _getMetric(metric):
    if (metric == 'D2'):
        metric = _D2

    elif (metric == 'D2S'):
        metric = _D2S

    elif (metric == 'D2Star'):
        metric = _D2Star

    elif (metric == 'd2S'):
        metric = _d2S

    elif (metric == 'd2Star'):
        metric = _D2Star

    return metric

def _D2(aKmerCount, bKmerCount):
    ## D2(X, Y) = sum(Xw * Yw)
    D2 = np.dot(aKmerCount, bKmerCount)
    return D2

def _d2(aKmerCount, bKmerCount):
    ## A = sqrt(sum(Xw^2))
    ## B = sqrt(sum(Yw^2))
    ## d2(X, Y) = 1/2 * (1 - (D2 / (A * B))
    D2 = _D2(aKmercount, bKmerCount)
    A  = math.sqrt(np.sum((aKmerCount ** 2)))
    B  = math.sqrt(np.sum((bKmerCount ** 2)))
    d2 = 1/2 * (1 - (D2 / (A * B)))

    ## Ensure that the distances range between 0 and 1 (as per the definition)
    d2 = 0 if d2 < 0 else d2
    d2 = 1 if d2 > 1 else d2
    return d2

def _D2S(aKmerCount, bKmerCount):
    ## Xw = Xw - N
    ##  N = (len(X) - k + 1) * P(w)
    ## Yw = Yw - M
    ##  M = (len(Y) - k + 1) * P(w)
    ## D2S(X, Y) = sum([Xw * Yw] / [sqrt((Xw^2 + Yw^2))])
    top    = np.multiply(aKmerCount, bKmerCount)
    bottom = np.add((aKmerCount ** 2), (bKmerCount ** 2))
    bottom = np.sqrt(bottom)

    D2S    = np.divide(top, bottom)
    D2S[np.isnan(D2S)] = 0    ## Replace all NaN's with 0's.
    D2S    = np.sum(D2S)
    return D2S

def _d2S(aKmerCount, bKmerCount):
    ##  N = (len(X) - k + 1) * P(w)
    ##  M = (len(Y) - k + 1) * P(w)
    ## Xw = Xw - N
    ## Yw = Yw - M
    ##  C = sqrt((Xw^2 + Yw^2))
    ##  A = sqrt(sum((Xw^2 / C)))
    ##  B = sqrt(sum((Yw^2 / C)))
    ## d2S(X, Y) = 1/2 * (1 - (D2S / (A * B)))
    D2S = _D2S(aKmerCount, bKmerCount)

    C   = np.add((aKmerCount ** 2), (bKmerCount ** 2))
    C   = np.sqrt(C)
    A   = np.divide((aKmerCount ** 2), C)
    A[np.isnan(A)] = 0    ## Replace all NaN's with 0's.
    A   = math.sqrt(np.sum(A))
    B   = np.divide((bKmerCount ** 2), C)
    B[np.isnan(B)] = 0    ## Replace all NaN's with 0's.
    B   = math.sqrt(np.sum(B))

    ## Ensure that the distances range between 0 and 1 (as per the definition)
    d2S = 1/2 * (1 - (D2S / (A * B)))
    d2S = 0 if d2S < 0 else d2S
    d2S = 1 if d2S > 1 else d2S
    return d2S

## I'm still not 100% sure whether this is the correct formula
## Assuming theres not too much difference between D2S and D2Star,
## Then its probably OK if we stick with the D2S
def _D2Star(aKmercount, bKmerCount):
    ##  N = (len(X) - k + 1) * P(w)
    ##  M = (len(Y) - k + 1) * P(w)
    ## Xw = (Xw - N) / sqrt(N)
    ## Yw = (Yw - M) / sqrt(M)
    ## D2Star(X, Y) = sum(Xw * Yw)
    return np.dot(aKmerCount, bKmerCount)

## We can't calculate this without passing in the P(W)
## So this function is currently not feasible!
def _d2Star(aKmerCount, bKmerCount):
    ## d2Star(X, Y) = 1/2 * (1 - (D2Star / (A * B)))
    ## A = sqrt(sum((Xw^2 / P(w))))
    ## B = sqrt(sum((Yw^2 / P(w))))
    D2Star = _D2Star(aKmerCount, bKmerCount)
    A = math.sqrt(np.sum(np.divide((aKmercount ** 2), C)))
    B = math.sqrt(np.sum(np.divide((bKmerCount ** 2), C)))
    d = 1/2 * (1 - (D2Star / (A * B)))
    return d

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
