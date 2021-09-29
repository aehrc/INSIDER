#!/bin/python

#------------------- Description & Notes --------------------#

'''
We can only visualise things if we load all the data into memory.
'''

#------------------- Dependencies ---------------------------#

# Standard library imports
import math

# External imports
import holoviews as hv
from holoviews import opts, dim
import matplotlib.pyplot as plt
import pandas as pd

# Internal imports

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def radar(kmerCount):
    ## Calculate the angle of each axis in the plot by dividing the plot
    ## by the number of variables
    nCols  = kmerCount.shape[1]
    angles = [n / float(nCols) * 2 * math.pi for n in range(nCols)]
    angles += angles[:1]    ## Complete the circle

    ## Initialise the spider plot, set the first axis to be on top,
    ## draw one axe per variable, add labels to each axes
    ax = plt.subplot(111, polar=True)
    ax.set_theta_offset(math.pi / 2)
    ax.set_theta_direction(-1)
    plt.xticks(angles[:-1], kmerCount.columns)

    ## Draw ylabels
    ax.set_rlabel_position(0)

    ## Plot each sequence
    for idx in kmerCount.index:
        v = kmerCount.loc[idx].tolist()
        v += v[:1]
        ax.plot(angles, v, linewidth=1, linestyle='solid', label=idx)
        ax.fill(angles, v, 'b', alpha=0.1)

    ## Add legend
    plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
    plt.show()

def histogram(kmerCount):
    ## Plot the counts
    kmerCount  = kmerCount.T
    kmerCount.plot.line(y=kmerCount.columns, legend=None)
    plt.show()

def pca(kmerPca, asBokeh=False):
    d = hv.Dataset(kmerPca, ['PCA1', 'PCA2'])
    s = d.to(hv.Scatter, ['PCA1', 'PCA2'], groupby='OutlierId').overlay()

    if (asBokeh):
        hv.extension('bokeh')
        s = s.opts(
            opts.Scatter(tools=['hover'], size=10, width=800, height=800,
                color='cramers_v', colorbar=True, cmap='bmy', logz=True,
                marker=dim('OutlierId').categorize({0:'circle', 1:'square'}),
                show_legend=True))

    else:
        hv.extension('matplotlib')
        s = s.opts(
                opts.Scatter(fig_size=300,
                    fontsize={'labels':20, 'xticks':10, 'yticks':10},
                    s=dim('OutlierId').categorize({0:25, 1:100}),
                    color='cramers_v', colorbar=True, cmap='bmy', logz=True,
                    marker=dim('OutlierId').categorize({0:'o', 1:'s'}),
                    show_legend=True))
    return s

#------------------- Protected Classes & Functions ------------#

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
