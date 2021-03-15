#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import functools
import random
import re

# External imports
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Internal imports

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getRandomRecord(featureDB, featureId):
    features = []

    if (featureId is None):
        featureId = _getRandomType(featureDB)

    ## Calculate interfeatures if we're looking for them
    ## the for first time
    if (featureId == 'inter_gene_gene'
        and featureId not in list(featureDB.featuretypes())):
            featureDB = _addInterfeatures(featureDB)    ## Experimental

    features = list(featureDB.features_of_type(featureId))

    randFeatureIdx = random.randint(0, len(features) - 1)
    randFeatureRec = features[randFeatureIdx]
    return randFeatureRec

def createSeqRecord(faRec, featureRec):
    if ('Name' in featureRec.attributes):
        fId = featureRec.attributes['Name'][0]

    elif ('ID' in featureRec.attributes):
        fId = featureRec.attributes['ID'][0]

    else:
        ## I don't know what we can use if they don't have any
        print(featureRec.attributes)
        raise NotImplementedError('No name detected')

    fType    = featureRec.featuretype
    startPos = featureRec.start
    endPos   = featureRec.stop
    strand   = featureRec.strand

    seq   = Seq(faRec[1][startPos:endPos])
    seq   = seq.reverse_complement() if (strand == '-') else seq
    seqId = "{}:{}:{}:{}:{}".format(faRec[0], fId, fType, startPos, endPos)
    return SeqRecord(seq, seqId, description='')

#------------------- Private Classes & Functions ------------#

def _getRandomType(featureDB):
    types = list(featureDB.featuretypes())

    randTypeIdx = random.randint(0, len(types) - 1)
    randType    = types[randTypeIdx]
    return randType

def _addInterfeatures(featureDB):
    ## GFFutils featureDB.interfeatures() does not work when
    ## there are features from different chromosomes or strand
    seqIdList  = _getSeqIds(featureDB)
    strandList = _getStrands(featureDB)
    features   = list(featureDB.features_of_type('gene'))

    ## There is an issue in that we are only looking at regions
    ## between gene features. Such regions MAY or MAY NOT contain
    ## features (i.e., non-coding exons, tRNAs etc...).
    ## In these cases, we need to either remove or shorten each feature
    for s in strandList:
        for seqId in seqIdList:
            f = lambda x: x.strand == s and x.seqid == seqId
            featuresSublist = filter(f, features)
            interfeatures   = featureDB.interfeatures(featuresSublist)
            interfeatures   = map(_updateInterfeatureName, interfeatures)
            featureDB.update(interfeatures)

    return featureDB

def _getSeqIds(featureDB):
    seqIdQuery  = 'SELECT DISTINCT seqid FROM features;'
    seqIdResult = featureDB.execute(seqIdQuery)

    f = lambda x, y: list(x) + list(y)
    seqIdList = list(functools.reduce(f, seqIdResult))
    return seqIdList

def _getStrands(featureDB):
    strandQuery  = 'SELECT DISTINCT strand FROM features \
                    WHERE strand != ".";'   ## Some features don't have a strand
    strandResult = featureDB.execute(strandQuery)

    f = lambda x, y: list(x) + list(y)
    strandList = list(functools.reduce(f, strandResult))
    return strandList

def _updateInterfeatureName(interfeature):
    f = lambda x: re.sub('^.*:', '', x)
    names = list(map(f, interfeature.attributes['ID']))

    name  = ['_'.join(names)]
    interfeature.attributes['Name'] = name
    return interfeature

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
