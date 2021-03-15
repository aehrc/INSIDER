#!/bin/python

#------------------- Description & Notes --------------------#

'''
If we know the whole host and foreign sequence, then we have 'should'
have access to the (population) distributions of the host and
foreign sequences. Therefore, we can assess how different the distributions
are by comparing the two distributions using the Chi-Squared test.
* Host Vs Foreign

In addition to the whole host and foreign sequence, if we know the
hybrid sequence, then we 'should' also have access to the (population)
distribution of the hybrid sequence. If the distributions for the host and
foreign distributions are substantially different, then the distribution
of the hybrid sequence could be bimodal. However, if the distributions
for the host and foreign sequences are not substantially unique, then the
distribution of the hybrid sequence could be similar to the distributions
of the foreign sequence and the host sequence and therefore difficult to
distinguish. Therefore, we can assess how different the distributions
are by comparing the three distributions using multiple Chi-Squared tests
* Host Vs Foreign
* Host Vs Hybrid
* Foreign Vs Hybrid

However, we are working under the assumption that we don't know anything
about the host and foreign sequence and hence, we don't know anything
about the hybrid sequence. Because of this, we don't have access to the
(population) distributions for the host, foreign or hybrid sequences. This is
quite reasonable, given that we would almost never have access to population
distributions anyway.

In the case of next-generation sequencing reads, we only know the
subsequences of the hybrid sequence. We can consider this as like
sampling from the (population) distribution of the hybrid sequence. Therefore,
we can only estimate the distribution of the hybrid sequence. In doing so,
we can only estimate the distribution of host and foreign sequences since
the distribution of the hybrid sequence is derived from the two distributions.

We can estimate the distributions of the host and foreign sequences by
considering each subsequence as an observation and averaging the subsequences
within a sample. Therefore, if we can (roughly) distinguish between
subsequences of the host sequence and subsequences of the foreign
sequence (via clustering), then we can assess how different the estimate
distributions are by comparing the two estimated distributions using
Hotelling's T-Squared test.
* Host Vs Foreign

Moreover, we can estimate the distribution of the hybrid sequence
by averaging all subsequences. Therefore, we can assess how different the
distributions are by comparing the three distributions using multiple
Hotelling's T-Squared tests.
* Host Vs Foreign
* Foreign Vs Average
* Host Vs Average
'''

#------------------- Dependencies ---------------------------#

# Standard library imports
import sys
import os
import functools

# External imports
import numpy as np
import pandas as pd
import holoviews as hv

# Internal imports
from src.kmer import plot as kmerplot
from src.kmer import analyse as kmeranalyse
from src.util import params
from src.util import spark
from src import io
from src import kmer
from src import ml

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

#------------------- Private Classes & Functions ------------#

def mainKmerFrequencies(obsDirs, cIdFile, oFile):
    print("Analysing frequencies")

    with spark.getSparkSession() as ss:
        with ss.sparkContext as sc:
            ## Instead of reading the Kmer frequencies back into memory
            ## we could extend the pipeline here so that all the computation
            ## and analysis is done in parallel. This would essentially be
            ## a merging of the two Dev scripts (along with a few
            ## other changes)

            ## Read Kmer frequencies
            print("Reading Kmers")
            kmerDf = io.kmer.read(*obsDirs)
            kmerDf = kmeranalyse.setup(kmerDf)

            ## Read cluster labels
            print("Reading cluster labels")
            cIdDf  = ss.read.csv(cIdFile, sep='\t', header=True)
            cIdDf  = cIdDf.select('id', 'ClusterId')

            ## Cluster Kmer frequencies
            ## We only need to do if we're going to look at cluster signatures
            print("Clustering Kmers")
            cKmerDf = kmer.transform.agg.clusterCountsByColumn(kmerDf, cIdDf)
            cKmerDf = cKmerDf.repartition(kmerDf.rdd.getNumPartitions(),
                kmer.SEQID_COL_NAME)
            cKmerDf.persist()
            cKmerDf.show()

            print("Calculating cluster info")
            csDf  = kmeranalyse.cluster.getClusterSizes(kmerDf, cIdDf)
            cmcDf = kmeranalyse.cluster.getClusterMeanCounts(kmerDf, cIdDf)

            ## Calculate Pvalue
            print("Calculating P-values")
            pvDf = kmeranalyse.stats.chi.getPvalues(cKmerDf, method='a')
            pvDf = kmeranalyse.stats.getSignificantPairs(pvDf)
            pvDf = kmeranalyse.stats.getZScore(pvDf, 'cramers_v')
            pvDf = pvDf.drop(columns=['id_y'])
            pvDf = pvDf.rename(columns={'id_x':'ClusterId','p_value':'chi_p_value'})

            ## Find outliers
            print("Finding outliers")
            oIdDf = kmeranalyse.outlier.getLabels(cKmerDf)
            oIdDf = oIdDf.rename(columns={'id':'ClusterId'})

            ## Generate PCA plots for visualisation
            print("Generating PCA")
            pcaDf = kmeranalyse.getPcaDf(cKmerDf)
            pcaDf = pcaDf.rename(columns={'id':'ClusterId'})

            ## Join the above tables. We shouldn't
            ## have any problems if every table is small
            print("Writing output")
            dfs = [cIdDf.toPandas(), csDf, cmcDf, oIdDf, pvDf]
            f   = lambda x, y: x.merge(y, on='ClusterId', how='left')
            df  = functools.reduce(f, dfs)
            df  = df.sort_values(by=['cramers_v']).reset_index(drop=True)

            ## Join the PCA table separately because we have an extra row for
            ## the estimated genome signature
            df  = df.merge(pcaDf, on='ClusterId', how='outer')
            cond = (df['ClusterId'] == 'Sum')
            df.loc[cond, [kmer.SEQID_COL_NAME, 'isSignificant']] = ['Sum', True]

            ## Sort and write the results table
            df.to_csv(oFile, sep='\t', index=False, na_rep='0')

def plotKmerFrequencies(iData, oFile):
    if (os.path.isdir(iData[0])):
        print("Plotting counts")
        with spark.getSparkSession() as ss:
            with ss.sparkContext as sc:
                ## Instead of reading the Kmer frequencies back into memory
                ## we could extend the pipeline here so that all the computation
                ## and analysis is done in parallel. This would essentially be
                ## a merging of the two Dev scripts (along with a few
                ## other changes)

                ## Read Kmer frequencies
                print("Reading Kmers")
                kmerDf = io.kmer.read(*obsDirs)
                kmerDf = kmeranalyse.setup(kmerDf)

                ## Get the estimated genome signature by summing the counts
                ## of each K-mer across all sequences
                eKmerDf    = kmer.transform.agg.sumCounts(kmerDf)
                eKmerCount = kmeranalyse.stats.getCountsPdf(eKmerDf)
                eKmerCount = eKmerCount.divide(eKmerCount.sum(axis=1), axis=0)
                print(eKmerCount)

                ## Get the counts for each sequence
                kmerCount = kmeranalyse.stats.getCountsPdf(kmerDf)
                kmerCount = kmerCount.divide(kmerCount.sum(axis=1), axis=0)
                print(kmerCount)

                ## If applicable, filter/select sequences
                # ids  = ['NODE_125_', 'NODE_123_', 'NODE_139_']
                ids  = ['NODE_128_', 'NODE_77_', 'NODE_109_']
                idxs = [idx for i in ids for idx in kmerCount.index if i in idx]
                kmerCount = kmerCount.loc[idxs, :]

                ## Make plots
                # kmerCount  = pd.concat([sKmerCount, kmerCount]).fillna(0)
                # kmerCount.to_csv('counts.txt', sep='\t')
                kmerplot.radar(kmerCount)
                # kmerplot.histogram(kmerCount)

                ## Alternatively transform the data into a single vector and plot
                # sKmerCount = pd.concat([sKmerCount] * kmerCount.shape[0])
                # kmerCount  = kmerCount.reset_index(drop=True).divide(sKmerCount.reset_index(drop=True))
                # print(kmerCount)

    else:
        print("Plotting PCA")
        for iFile in iData:
            kmerPca = pd.read_csv(iFile, sep='\\t')
            kmerPca['OutlierId'] = kmerPca['OutlierId'].apply(lambda x: 0 if x < 1 else 1)
            kmerPca = kmerPca.sort_values(by=['cramers_v', 'OutlierId'], ascending=False)
            print(kmerPca)
            s = kmerplot.pca(kmerPca)
            hv.save(s, oFile)

def main():
    argParser = params.AnalyseKmerArgParser()
    argParser.parse()

    if (argParser.cmd == 'plot'):
        argParser.parsePlotArgs()
        argParser.printArgs()

        plotKmerFrequencies(argParser.iData, argParser. oFile)

    elif (argParser.cmd == 'main'):
        argParser.parseMainArgs()
        argParser.printArgs()

        mainKmerFrequencies(argParser.obsDirs, argParser.cIdFile,
            argParser.oFile)

    print("DONE")

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
