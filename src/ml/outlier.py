#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import numpy as np
import pandas as pd
from sklearn.covariance import EllipticEnvelope
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
from sklearn.svm import OneClassSVM

# Internal imports
from .constants import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def detect(kmerCount, algo, **kwargs):
    algos   = _getAlgorithms(algo, **kwargs)
    labels  = _getLabels(algos, kmerCount) 

    args = ['{}:{}'.format(k, v) for k, v in kwargs.items()]
    args = '_'.join(args)
    colName = '{}_{}_{}'.format(labels.columns.tolist()[0], algo, args)
    labels.columns = [colName]
    return labels

#------------------- Private Classes & Functions ------------#

def _getAlgorithms(algo, **kwargs):
    if (algo == 'EE'):
        algos = [EllipticEnvelope(**kwargs)]

    elif (algo == 'SVM'):
        algos = [OneClassSVM(**kwargs)]

    elif (algo == 'IF'):
        algos = [IsolationForest(**kwargs)]

    elif (algo == 'LOF'):
        algos = [LocalOutlierFactor(**kwargs)]

    else:
        raise NotImplementedError('Unknown algorithm')

    return algos

def _getLabels(algos, kmerCount):
    f  = lambda x: _fitPredict(x, kmerCount)
    ll = list(map(f, algos))
    ## Return the first prediction
    return ll[0][1]

def _fitPredict(algo, kmerCount):
    l = algo.fit_predict(kmerCount)
    l = np.reshape(l, (-1, 1))
    l = pd.DataFrame(l, columns=[OLABEL_COL_NAME])
    return (algo, l)

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------

# def investigateOptimalAlgorithms(kmerId, kmerPca):
#     plot.setLibrary('bokeh')

#     pca   = getMLColumns(kmerPca, FLABEL_COL_NAME)
#     plots = {}
#     algos = (
#         ('Elliptic', EllipticEnvelope()),
#         ('SVM', OneClassSVM()),
#         ('Forest', IsolationForest()),
#         ('Local', LocalOutlierFactor()))

#     ## Visualise data and manually determine which algorithm will be good
#     for i, (name, algo) in enumerate(algos, 1):
#         labels  = _getLabels(algo, pca)
#         kmerDf  = pd.concat([kmerId, pca, labels], axis=1)

#         dataset = hv.Dataset(kmerDf, dataCols)
#         scatter = dataset.to(hv.Scatter, dataCols, groupby=OLABEL_COL_NAME).overlay()
#         scatter.opts(opts.Scatter(size=10, show_legend=True))
#         plots[name] = scatter

#     plots = hv.HoloMap(plots, kdims='algo')
#     plots = plots.collate()
#     return plots
