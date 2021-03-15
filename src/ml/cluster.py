#!/bin/python

#------------------- Description & Notes --------------------#

'''
There are different approaches to model selection
See: https://stackoverflow.com/questions/39920862/model-selection-for-gaussianmixture-by-using-gridsearch

Unlike other clustering algorithms that tell us which data point belongs to
which cluster, Gaussian mixture is a density estimation algorithm that tells
us the probabilities that a data point belongs to each cluster.

Spectral Clustering looks really good on filtered data, but very inefficient!
DBSCAN is basically a faster version of Spectral Clustering, whereas OPTICS
is a more efficient version of DBSCAN. Unfortunately, I always encounter
a division by 0 error with OPTICS!

'''

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram
from sklearn import cluster
from sklearn import mixture

# Internal imports
from .constants import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def assign(kmerCount, algo, **kwargs):
    algos  = _getAlgorithms(algo, **kwargs)
    labels = _getLabels(algos, kmerCount)

    args = ['{}:{}'.format(k, v) for k, v in kwargs.items()]
    args = '_'.join(args)
    colName = '{}_{}_{}'.format(labels.columns.tolist()[0], algo, args)
    labels.columns = [colName]
    return labels

def visualiseHierarchicalClusters(kmerCount, **kwargs):
    algo  = _getAlgorithms(algo='AGGLO', **kwargs)[0]
    algo  = _fitPredict(algo, kmerCount)[0]
    d     = _getDendrogram(algo, truncate_mode='level', p=3)
    return d

#------------------- Private Classes & Functions ------------#

def _getAlgorithms(algo, **kwargs):
    if (algo == 'KMEANS'):
        algos = [cluster.KMeans(**kwargs)]

    elif (algo == 'DBSCAN'):
        algos = [cluster.DBSCAN(**kwargs)]

    elif (algo == 'OPTICS'):
        algos = [cluster.OPTICS(**kwargs)]

    elif (algo == 'SPEC'):
        algos = [cluster.SpectralClustering(**kwargs)]

    elif (algo == 'AGGLO'):
        algos = [cluster.AgglomerativeClustering(**kwargs)]

    elif (algo == 'GMM'):
        algos = [mixture.GaussianMixture(**kwargs)]

    else:
        raise NotImplementedError('Unknown algorithm')

    return algos

def _getLabels(algos, kmerCount):
    f  = lambda x: _fitPredict(x, kmerCount)
    ll = list(map(f, algos))
    ## Return the first prediction
    return ll[0][1]

def _fitPredict(algo, kmerCount):
    algo.fit(kmerCount)

    if (hasattr(algo, 'labels_')):
        l = algo.labels_.astype(np.int)

    else:
        l = algo.predict(kmerCount)
        # print(algo.cluster_centers_)

    l = np.reshape(l, (-1, 1))
    l = pd.DataFrame(l, columns=[CLABEL_COL_NAME])
    return (algo, l)

def _getDendrogram(model, **kwargs):
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)

    return dendrogram(linkage_matrix, **kwargs)

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------#

## Previous clustering code
# import holoviews as hv
# from holoviews import opts
# import matplotlib.pyplot as plt
# from sklearn.metrics import silhouette_score
# from sklearn.model_selection import GridSearchCV
# from .. import plot
# from .common import getParameterGrid
# from .common import getFeatureColumns

# ## Maximum number of clusters
# MAX_NUM_CLUSTERS = 50
# MAX_NUM_INIT     = 20

# def investigateOptimalAlgorithms(kmerId, kmerPca):
#     plot.setLibrary('bokeh')

#     pca   = getFeatureColumns(kmerPca)
#     plots = {}
#     algos = (
#         ('KMeans', cluster.KMeans()),
#         ('Affinity', cluster.AffinityPropagation()),
#         ('MeanShift', cluster.MeanShift()),
#         ('Spectral', cluster.SpectralClustering()),
#         ('Agglomerative', cluster.AgglomerativeClustering(linkage='average')),
#         ('Agglomerative', cluster.AgglomerativeClustering(linkage='ward')),
#         ('DBSCAN', cluster.DBSCAN()),
#         ('Gaussian', GaussianMixture()))

#     ## Visualise data and manually determine which algorithm will be good
#     for i, (name, algo) in enumerate(algos, 1):
#         labels  = _getLabels(algo, pca)
#         labels  = pd.DataFrame(labels, columns=[CLABEL_COL_NAME])
#         kmerDf  = pd.concat([kmerId, pca, labels], axis=1)

#         dataset  = hv.Dataset(kmerDf, dataCols)
#         scatter  = dataset.to(hv.Scatter, dataCols, groupby=CLABEL_COL_NAME).overlay()
#         scatter.opts(opts.Scatter(size=10, show_legend=True))
#         plots[name] = scatter

#     plots = hv.HoloMap(plots, kdims='algo')
#     plots = plots.collate()
#     return plots

# def _getScore(algo, pca):
#     cLabels = _getLabels(algo, pca)
#     if (isinstance(algo, cluster.KMeans)
#         or isinstance(algo, cluster.DBSCAN)
#         or isinstance(algo, cluster.SpectralClustering)):

#         score = silhouette_score(pca, cLabels)

#     elif (isinstance(algo, GaussianMixture)):
#         ## Use only one
#         score = algo.bic(pca)
#         # score = algo.aic(pca)

#     else:
#         raise NotImplementedError("Scoring for clustering algorithm not implemented.")

#     score = float(score)
#     return score

# def _getOptimalAlgorithm(pca, plotPerf=False):
#     ## Parameters need to be updated depending on the algorithm we use
#     params = getParameterGrid(
#         ('n_clusters', list(range(2, 10))),
#         # ('n_components', list(range(2, MAX_NUM_CLUSTERS))),
#         ('n_init', [MAX_NUM_INIT]))

#     ## Algorithm needs to be updated depending on their suitability
#     algos  = [cluster.KMeans(**p) for p in params]
#     # algos = ((p, OPTICS(**p)) for p in params)
#     # algos = ((p, cluster.AgglomerativeClustering(**p)) for p in params)
#     # algos = ((p, cluster.AgglomerativeClustering(**p)) for p in params)
#     # algos = ((p, GaussianMixture(**p)) for p in params)
#     # algos = ((p, cluster.DBSCAN(**p)) for p in params)
#     # algos = ((p, cluster.SpectralClustering(**p))) for p in params)

#     ## Find the optimal algorithm
#     f = lambda x: (x, _getScore(x, pca))
#     algos = list(map(f, algos))
#     if (plotPerf):
#         scores = [(p, {'score':algo[1]}) for p, algo in zip(params, algos)]
#         scores = [{**s[0], **s[1]} for s in scores]
#         sDf    = pd.DataFrame(scores)
#         sDf.plot(x='n_clusters', y='score')
#         plt.show()

#     algo = sorted(algos, key=lambda x: x[1], reverse=True)[0][0]
#     return algo

# class GaussianModel(Model):

#     def __init__(self, pca):
#         Model.__init__(self, pca)

#     def getHyperParameters(self):
#         ## The only hyperparameter we are interested
#         ## is the number of clusters to use
#         nClusters = range(2, self.MAX_NUM_CLUSTERS)
#         params    = list(product(eps, samples))
#         return params

#     def getMethod(self, param):
#         nClusters = param[0]  ## We know that nClusters is the only parameter

#         cMethod   = GaussianMixture(n_components=nClusters, n_init=10)
#         cMethod   = cMethod.fit(self.pca)
#         return cMethod

#     def calculateScores(self, cMethod):
#         ## Supposedly metrics for assessing the number of components to use
#         ## https://jakevdp.github.io/PythonDataScienceHandbook/05.12-gaussian-mixtures.html
#         BICScore  = cMethod.bic(self.pca)
#         AICScore  = cMethod.aic(self.pca)
#         nClusters = cMethod.n_components

#         scores   = pd.DataFrame([[nClusters, BICScore, AICScore]])
#         scores.columns = [NCLUSTERS_COL_NAME, BIC_COL_NAME, AIC_COL_NAME]
#         return scores

#     def groupScores(self, scores):
#         scores = scores.groupby([NCLUSTERS_COL_NAME], as_index=False).mean()
#         return scores

#     def plotPerformance(self):
#         self.addGradientColumn()

#         ## Plot BIC/AIC scores
#         y = [BIC_COL_NAME, AIC_COL_NAME]
#         self.perfScores.plot(x=NCLUSTERS_COL_NAME, y=y)
#         plt.show()

#         ## Plot BIC/AIC gradients
#         y = [BIC_GRADIENT_COL_NAME, AIC_GRADIENT_COL_NAME]
#         self.perfScores.plot(x=NCLUSTERS_COL_NAME, y=y)
#         plt.show()

#     def addGradientColumn(self):
#         BICScores = self.perfScores.loc[:, BIC_COL_NAME]
#         AICScores = self.perfScores.loc[:, AIC_COL_NAME]
#         self.perfScores[BIC_GRADIENT_COL_NAME] = np.gradient(BICScores)
#         self.perfScores[AIC_GRADIENT_COL_NAME] = np.gradient(AICScores)
#         return self.perfScores

#     def selectOptimalMethod(self):
#         ## Choose the number of clusters based on the BIC score
#         ## BIC seems a bit better / more popular than the AIC score

#         ## Approch 1:
#         # ## Lower BIC score = better
#         # minBicIdx = self.perfScores[BIC_COL_NAME].idxmin()
#         # nClusters = self.perfScores.loc[minBicIdx, NCLUSTERS_COL_NAME].astype(int)

#         ## Approach 2:
#         ## Automating the 'Elbow Rule'
#         ## https://www.datasciencecentral.com/profiles/blogs/how-to-automatically-determine-the-number-of-clusters-in-your-dat
#         bicColumn = self.perfScores.loc[:, BIC_COL_NAME]
#         delta1    = bicColumn.diff() * -1
#         delta2    = delta1.diff() * -1
#         strength  = delta2 - delta1
#         strength  = strength.shift(-1)

#         idx       = strength.idxmax()   ## Higher strength == better
#         nClusters = self.perfScores.loc[idx, NCLUSTERS_COL_NAME].astype(int)
#         cMethod   = self.getMethod((nClusters,))
#         print("Num. Clusters:\t" + str(nClusters))
#         return cMethod
