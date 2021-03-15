#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import itertools

# External imports
import holoviews as hv
from holoviews import dim         ## Requires python 3.7; not 3.5
from holoviews import opts

# Internal imports

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getPlots(plots, pDict):
    pNames = list(pDict.keys())
    pVals  = list(pDict.values())
    pVals  = itertools.product(*pVals)

    ## Create the map of plots
    plots = {a:b for a, b in zip(pVals, plots)}
    plots = hv.HoloMap(plots, kdims=pNames)
    plots = plots.collate()
    return plots

def setLibrary(library='bokeh'):
    ## Use Bokeh by default
    if (library == 'bokeh'):
        hv.extension('bokeh')
        # hv.archive.auto(filename_formatter="{obj:.7}")    ## For notebooks

        opts.defaults(
            opts.Curve(tools=['hover', 'tap'], padding=0.05,
                width=700, height=700),
            opts.Scatter(tools=['hover', 'tap'], padding=0.05,
                width=700, height=700),
            opts.HeatMap(tools=['hover', 'tap'], labelled=[], xrotation=45,
                width=700, height=700, colorbar=True, cmap=('Blues')),
            opts.Graph(tools=['hover', 'tap'], padding=0.05,
                width=700, height=700, xaxis=None, yaxis=None))

        ## The library that Bokeh uses to export to SVG is not longer supported
        ## and so cannot be exported to SVG

    elif (library == 'matplotlib'):
        hv.extension('matplotlib')
        hv.output(fig='svg')

        opts.defaults(
            opts.Curve(fig_size=300, padding=0.05),
            opts.Scatter(fig_size=300, padding=0.05),
            opts.HeatMap(fig_size=300, labelled=[], xrotation=45,
                colorbar=True, cmap=('Blues')),
            opts.Graph(fig_size=300, padding=0.05))

    elif (library == 'plotly'):
        hv.extension('plotly')
        opts.defaults(
            opts.Curve(padding=0.05, width=900, height=900),
            opts.Scatter(padding=0.05, width=900, height=900),
            opts.Scatter3D(padding=0.05, width=900, height=900))

    else:
        raise NotImplementedError("Unknown plotting library.")

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
