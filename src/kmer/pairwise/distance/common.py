#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import numpy as np
import pandas as pd
import pyspark.sql.functions as F
import pyspark.sql.types as T
from pyspark.ml.feature import MinMaxScaler
from pyspark.ml.feature import VectorAssembler
from pyspark.ml import Pipeline
from scipy.spatial.distance import squareform

# Internal imports
from ..common import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def tableToSymMatrix(kmerDistDf):
    ## Construct a list of distances from a sorted dictionary
    rows  = kmerDistDf.to_dict('records')
    rows  = {tuple(sorted([r['seqId_1'], r['seqId_2']])):r['distance'] for r in rows}
    dists = [dist[1] for dist in sorted(rows.items())]

    ## Get all the IDs and construct the matrix
    ids = list(pd.concat([kmerDistDf['seqId_1'], kmerDistDf['seqId_2']]).unique())
    kmerDistDf = pd.DataFrame(squareform(dists), index=sorted(ids), columns=sorted(ids))
    return kmerDistDf

def matrixToSymMatrix(kmerDistDf, useMax=True):
    triu = np.triu(kmerDistDf)
    tril = np.tril(kmerDistDf)

    ## Compare the triangles and get the maximum value for each pair
    newTril = tril.T
    newTriu = np.maximum(triu, newTril) if useMax else np.minimum(triu, newTril)

    ## Reconstruct the matrix
    m = newTriu + newTriu.T
    np.fill_diagonal(m, 0)
    kmerDistDf = pd.DataFrame(m, columns=kmerDistDf.columns,
        index=kmerDistDf.index)

    return kmerDistDf

def scale(kmerDistDf, scale=(0, 100)):
    va   = VectorAssembler(inputCols=['distance'], outputCol="v")
    mms  = MinMaxScaler(min=scale[0], max=scale[1], inputCol="v", outputCol="s")
    pipe = Pipeline(stages=[va, mms])

    ## Perform scaling pipeline
    m          = pipe.fit(kmerDistDf)
    kmerDistDf = m.transform(kmerDistDf)

    ## Reorganise and cleanup
    f = lambda x: float(x[0])
    f = F.udf(f, T.FloatType())
    kmerDistDf = kmerDistDf.withColumn('distance', f('s')) \
        .select('seqId_1', 'seqId_2', 'distance')
    return kmerDistDf

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
