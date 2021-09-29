#!/bin/python

#------------------- Description & Notes --------------------#

'''
tSNE calculates the similarity between pairs of vectors in both the
high dimensional space and the low dimensional space. It then tries to optimise
these two similarity measures using a cost function (KL-divergence)

Similarity between pairs of vectors in high dimensional space is measured
by its density and probabilities using the Gaussian distribution. Vectors
that are similar have similar probabilities. The distribution varies according
to the perplexity.

In the same way, the similarity between pairs of vectors in low dimensional
space is measured using a Student t-distribution (i.e., Cauchy distribution).

The above probability distributions are then compared by measuring the
KL divergence. A relatively low KL divergence indicates that two distributions
are similar. However, unlike distance metrics, KL divergence is not a symmetric
measure.

Unlike PCA, the components of tSNE don't really have any meaning.
The important thing is the relative distance between points. But,
keep in mind they're probabilities and not distances.

Instead, the performance of tSNE is measured by the KL divergence
which measures the similarity between the two probability
distributions used for the dimensionality reduction

tSNE seems to represent both global and local differences well
Both MDS and PCA appear to represent global differences well. However,
they don't seem to represent the local differences very well.

tSNE is good for preserving local/relative structure and therefore generates
very informative visualisations of the heterogeneity in the data. However, it
does not guarantee that inter-class distances are preserved correctly, making
it difficult for us cluster and understand something about relatedness between
clusters.
'''

#------------------- Dependencies ---------------------------#

# Standard library imports
import os
import sys

# External imports
import pandas as pd
import holoviews as hv
from pyspark.sql import SparkSession
import pyspark.sql.functions as sparkFsql
import pyspark.ml.functions as sparkFml
from pyspark.ml.feature import PCA
from pyspark.ml.feature import VectorAssembler
from pyspark.ml import Pipeline
from sklearn import decomposition
from sklearn import preprocessing
from sklearn import manifold
from sklearn import random_projection

# Internal imports
from .. import kmer
from .constants import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def sklearnReduce(kmerCount, frMethod, **kwargs):
    algo = _getAlgorithm(frMethod, **kwargs)

    ## Perform feature reduction
    comps = algo.fit_transform(kmerCount)
    comps = _formatComponents(comps)
    return comps

def fastTsneReduce(kmerCount, **kwargs):
    ## We need a better way of doing this
    sys.path.append('../../software/FIt-SNE/')
    # p = os.path.join(os.environ['SCRATCH1DIR'], 'FIt-SNE')
    # sys.path.append(p)
    from fast_tsne import fast_tsne

    comps = fast_tsne(kmerCount, **kwargs)
    comps = _formatComponents(comps)
    return comps

def sparkPcaReduce(kmerCount, **kwargs):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    df = ss.createDataFrame(kmerCount)

    ## Convert kmerCount into Spark Vectors and setup the pipeline
    k    = kwargs.pop('k', 2)
    va   = VectorAssembler(inputCols=df.columns, outputCol=kmer.KMER_COL_NAME)
    pca  = PCA(inputCol=kmer.KMER_COL_NAME, outputCol=FLABEL_COL_NAME, k=k)
    pipe = Pipeline(stages=[va, pca])

    ## Perform PCA pipeline
    m     = pipe.fit(df)
    comps = m.transform(df)

    ## Gather the components
    ## This will load the reduced feature table into memory
    f     = lambda x: sparkFml.vector_to_array(x)
    g     = lambda x: [sparkFsql.col(x)[i] for i in range(0, k)]
    comps = comps.withColumn(FLABEL_COL_NAME, f(FLABEL_COL_NAME)) \
                 .select(g(FLABEL_COL_NAME)) \
                 .toPandas()

    colNames = [FLABEL_COL_NAME + str(i + 1) for i in range(0, comps.shape[1])]
    comps.columns = colNames
    return comps

def sparkIncrementalPcaReduce(kmerCount, **kwargs):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ipca = decomposition.IncrementalPCA(**kwargs)

    ## Fit-transform
    for x in kmerCount.toLocalIterator():
        ipca.partial_fit(x)
    comps = kmerCount.map(lambda x: ipca.transform(x))

    ## Scale values
    scaler = preprocessing.MinMaxScaler((0, 1))
    for x in comps.toLocalIterator():
        scaler.partial_fit(x)

    comps = comps.map(lambda x: scaler.transform(x)) \
        .map(_formatComponents)

    print("Explained Variance\t{}".format(ipca.explained_variance_))
    print("Explained Variance Ratio\t{}".format(ipca.explained_variance_ratio_))
    return comps

def visualisePcaPerformance(kmerCount, **kwargs):
    perf = _getPcaPerformance(kmerCount, **kwargs)

    ## Create the plot
    d  = hv.Dataset(perf, *perf.columns)
    c  = d.to(hv.Curve, *perf.columns)
    t  = d.to(hv.Table, *perf.columns)
    p  = (c + t).cols(2)
    return p

def visualiseTsnePerformance(kmerCount, **kwargs):
    perf = _getTsnePerformance(kmerCount, **kwargs)

    ## Create the plot
    d  = hv.Dataset(perf, *perf.columns)
    c  = d.to(hv.Curve, *perf.columns)
    t  = d.to(hv.Table, *perf.columns)
    p  = (c + t).cols(2)
    return p

def scale(prevComps, scale=(0, 1)):
    scaler    = preprocessing.MinMaxScaler(scale)
    currComps = scaler.fit_transform(prevComps.to_numpy())
    currComps = pd.DataFrame(currComps,
        columns=prevComps.columns, index=prevComps.index)
    return currComps

#------------------- Private Classes & Functions ------------#

def _getAlgorithm(frMethod, **kwargs):
    if (frMethod == 'PCA'):
        algo = decomposition.PCA(**kwargs)

    elif (frMethod == 'SVD'):
        algo = decomposition.TruncatedSVD(**kwargs)

    elif (frMethod == 'FA'):
        algo = decomposition.FactorAnalysis(**kwargs)

    elif (frMethod == 'IPCA'):
        algo = decomposition.IncrementalPCA(**kwargs)

    elif (frMethod == 'KPCA'):
        algo = decomposition.KernelPCA(**kwargs)

    elif (frMethod == 'TSNE'):
        ## Find the optimal tSNE configuration
        algo = manifold.TSNE(**kwargs)
        # from MulticoreTSNE import MulticoreTSNE
        # algo = MulticoreTSNE(**kwargs)

    elif (frMethod == 'MDS'):
        algo = manifold.MDS(**kwargs)

    elif (frMethod == 'LLE'):
        algo = manifold.LocallyLinearEmbedding(**kwargs)

    elif (frMethod == 'ISOMAP'):
        algo = manifold.Isomap(**kwargs)

    elif (frMethod == 'SPECE'):
        algo = manifold.SpectralEmbedding(**kwargs)

    elif (frMethod == 'RP'):
        algo = random_projection.SparseRandomProjection(**kwargs)

    elif (frMethod == 'UMAP'):
        ## Requires the umap package
        import wpca
        algo = umap.UMAP(**kwargs)

    elif (frMethod == 'WPCA'):
        ## Requires the wpca package
        import umap
        algo = wpca.WPCA(**kwargs)

    else:
        raise NotImplementedError('Unknown algorithm')

    return algo

def _formatComponents(c):
    colNames = [FLABEL_COL_NAME + str(i + 1) for i in range(0, c.shape[1])]
    c        = pd.DataFrame(c, columns=colNames)
    return c

def _getPcaPerformance(kmerCount, **kwargs):
    algo = _getAlgorithm('PCA', **kwargs)
    if (algo.n_components > kmerCount.shape[1]):
        raise ValueError('Cannot compute explained variance')

    algo.fit_transform(kmerCount)

    ## Components of PCA are a linear combination of
    ## the feature variances in the the input data
    colNames = ['Components', 'ExplainedVariance']
    perf     = algo.explained_variance_ratio_
    perf     = pd.DataFrame(enumerate(perf, 1), columns=colNames)
    return perf

def _getTsnePerformance(kmerCount, **kwargs):
    n_components = kwargs.get('n_components')
    if (n_components > kmerCount.shape[1]):
        raise ValueError('Cannot compute KL divergence')

    ## Use the exact method for better performance estimation
    algos = [_getAlgorithm('TSNE', n_components=n, method='exact')
             for n in range(1, n_components + 1)]
    [algo.fit_transform(kmerCount) for algo in algos]

    colNames = ['Components', 'KL divergence']
    perf     = [algo.kl_divergence_ for algo in algos]
    perf     = pd.DataFrame(enumerate(perf, 1), columns=colNames)
    return perf

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

## Previous feature code
# from sklearn.ensemble import ExtraTreesClassifier
# from sklearn.pipeline import make_pipeline

# def investigateAlgorithms(kmerId, kmerCount, **kwargs):
#     n_components = kwargs.get('n_components', DEFAULT_N_COMPONENTS)
#     random_state = kwargs.get('random_state', DEFAULT_RANDOM_STATE)

#     params = {'n_components':n_components, 'random_state':random_state}
#     algos  = (
#         ('PCA', decomposition.PCA(**params)),
#         ('Modified LLE', manifold.LocallyLinearEmbedding(method='modified',
#             **params)),
#         ('tSNE', manifold.TSNE(**params)),
#         ('Isomap', manifold.Isomap(n_components=n_components)),
#         ('MDS', manifold.MDS(**params)),
#         ('SE', manifold.SpectralEmbedding(**params)))

#     ## Visualise data and manually determine which algorithm will be good
#     for i, (name, algo) in enumerate(algos, 1):
#         comps    = _getComponents(algo, kmerCount)
#         kmerDf   = pd.concat([kmerId, comps], axis=1)
#         dataCols = list(filter(lambda x: FLABEL_COL_NAME in x, kmerDf.columns))

#         dataset  = hv.Dataset(kmerDf, dataCols)
#         scatter  = dataset.to(hv.Scatter, dataCols)
#         scatter.opts(opts.Scatter(show_legend=True))
#         plots[name] = scatter

#     plots = hv.HoloMap(plots, kdims='algo')
#     plots = plots.collate()
#     return plots

# def investigateTsneParameters(kmerCount, **kwargs):
#     perp   = kwargs.get('perplexity', [2, 5, 10, 50, 100, 300, 500])
#     params = getParameterGrid(
#         ('perplexity', perp))

#     ## Algorithm needs to be updated depending on their suitability
#     algos     = [manifold.TSNE(**p) for p in params]

#     ## Find the optimal algorithm
#     f         = lambda x: (x, _getTsneScore(x, kmerCount))
#     algos     = list(map(f, algos))
#     scores    = [(p, {'score':algo[1]}) for p, algo in zip(params, algos)]
#     scores    = [{**s[0], **s[1]} for s in scores]
#     paramPerf = pd.DataFrame(scores)

#     ## Create the plot
#     d  = hv.Dataset(paramPerf, *paramPerf.columns)
#     c  = d.to(hv.Curve, *paramPerf.columns)
#     t  = d.to(hv.Table, *paramPerf.columns)
#     p  = (c + t).cols(2)
#     return p

# def getImportance(kmerId, kmerCount):
#     ## Determine the importance of each feature
#     ## VarianceThreshold doesn't work very well
#     ## Especially since I got no idea what threshold to use
#     algo       = ExtraTreesClassifier(n_estimators=50)
#     fittedAlgo = algo.fit(kmerCount, kmerId[kmer.SEQID_COL_NAME])
#     ranks      = fittedAlgo.feature_importances_

#     ## Construct a dataframe of Kmers and their relative importance
#     rows = [(r, k) for r, k in zip(ranks, kmerCount.columns)]
#     df   = pd.DataFrame(rows, columns=['importance', 'kmer'])
#     return df

# def select(kmerId, kmerCount):
#     df   = getImportance(kmerId, kmerCount)
#     cond = (df['importance'] != 0)
#     cols = df[cond]['kmer'].values
#     kmerCount = kmerCount[cols]
#     return kmerCount

