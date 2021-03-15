#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import numpy as np
import pandas as pd
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity

# Internal imports
from ..util import spark
from .. import kmer
from .constants import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def run(rKmerId, rPca, *tPca):
    fNames = rKmerId[FILE_COL_NAME].unique().tolist()
    rPcas  = [getPcaSubset(rKmerId, rPca, f) for f in fNames]
    (rDen, bws, kdes) = fit(fNames, rPcas)

    if (len(tPca) != 0):
        tDen = calculateProbability(kdes, tPca[0])
        tDen = findMax(tDen)
        return (rDen, tDen, [DPROB_COL_NAME, DLABEL_COL_NAME], bws)

    return (rDen, [DPROB_COL_NAME], bws)

#------------------- Private Classes & Functions ------------#

def fit(fNames, pcas):
    ## Calculate the density on each data subset
    bws  = [(fName, selectBandwidth(pca)) for fName, pca in zip(fNames, pcas)]
    kdes = [getKde(pca, bw[1]) for pca, bw in zip(pcas, bws)]
    zs   = [getScore(pca, kde) for pca, kde in zip(pcas, kdes)]
    den  = np.concatenate(zs)
    return (den, bws, kdes)

def getPcaSubset(kmerId, pca, f):
    isFile       = (kmerId[FILE_COL_NAME] == f)
    kmerIdSubset = kmerId[isFile]
    pcaSubset    = pca[kmerIdSubset.index.values, :]
    return pcaSubset

def selectBandwidth(pca):
    ## 20-fold cross-validation
    grid = GridSearchCV(KernelDensity(kernel='gaussian'),
        {'bandwidth': np.linspace(0.1, 1.0, 30)}, cv=10)
    grid.fit(pca)
    return grid.best_params_['bandwidth']

def getKde(pca, bw):
    kde = KernelDensity(kernel='gaussian', bandwidth=bw)
    kde = kde.fit(pca)
    return kde

def getScore(pca, kde):
    z = np.exp(kde.score_samples(pca))
    z = np.reshape(z, (-1, 1))
    return z

def calculateProbability(kdes, pca):
    ## Calculate the score (i.e., probability) of new data points
    ## https://stackoverflow.com/questions/24681825/kernel-density-score-vs-score-samples-python-scikit
    scores  = [np.exp(k.score_samples(pca)) for k in kdes]
    scores  = [np.reshape(s, (-1, 1)) for s in scores]
    scores  = np.hstack((scores[0], scores[1]))
    probSum = np.sum(scores, axis=1)

    prob1 = np.divide(scores[:, 0], probSum, where=probSum != 0)
    prob1 = np.reshape(prob1, (-1, 1))
    prob2 = np.divide(scores[:, 1], probSum, where=probSum != 0)
    prob2 = np.reshape(prob2, (-1, 1))
    probs = np.hstack((prob1, prob2))
    return probs

def findMax(probs):
    idxs  = np.argmax(probs, axis=1)
    idxs  = np.reshape(idxs, (-1, 1))
    # tDens = np.hstack((probs, idxs)) ## Probabilities & Labels
    return idxs

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
