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

'''
We can calculate pairwise distances using the idea from the following link.
It seems to work well if we don't need the whole matrix (but not what we want)
https://medium.com/@rantav/large-scale-matrix-multiplication-with-pyspark-or-how-to-match-two-large-datasets-of-company-1be4b1b2871e
'''

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports

# Internal imports
from .adjusted import getDistances
from .traditional import getDistances

from .common import tableToSymMatrix
from .common import matrixToSymMatrix
from .common import scale

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
