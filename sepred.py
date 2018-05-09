#!/usr/bin/env python

"""%prog [options]

implementation and testing of correlation analysis found here:
Atias, N, Sharan, R (2011). An algorithmic framework for predicting side effects of drugs. J. Comput. Biol., 18, 3:207-18.

A class of objects (drugs) are described by $p$ attributes.
We wish to obtain a model to predict $m$ properties (side effects) of these objects,
given a training set of $n$ objects (drugs) and their associated
attributes (a pxn matrix R) and properties (an mxn matrix E).

           predicted
actual   | negative | positive
---------.----------|--------
negative |    a     |   b
positive |    c     |   d

recall:    r = d / (c + d)
precision: p = d / (b + d)
geometric mean:  sqrt(r p)

acc+ = d / (c + d)
acc- = a / (a + b)
g-mean: g = sqrt(acc+ acc-)
"""

import sys
sys.path.append('~/tools/sepred')
from optparse import OptionParser
import numpy as np
import sepred_utils as sru
import sepred_reader_binary_effects as srbe
import sepred_reader_effects as sre
import cPickle as pickle

parser = OptionParser(__doc__)

parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="verbose output")

parser.add_option("-a", "--attributes",
                  dest="attributesfile",
                  help="read in the attribute information for the dataset from ATTRIBUTES",
                  metavar="ATTRIBUTES",
                  default="nci60/nci_compounds_descriptors.csv")
                  # default="drugbank_approved_small_descp.txt")

parser.add_option("-e", "--effects",
                  dest="effectsfile",
                  help="read in the effect information for the dataset from EFFECTS",
                  metavar="EFFECTS",
                  default="nci60/GI50_RAW.txt")
                  # default="Drugs_Approved_Side_Effects_New.csv")

parser.add_option("-k", "--projected-size",
                  dest="projected_size",
                  help="size of projected space, an integer K",
                  metavar="K",
                  default="0")

(options, args) = parser.parse_args()

### treating effects as binary values gives terrible results, just treat
### effects as regular (continuous, non-binary) values going forward.
### binary effects code is candidate for removal from program
binary_effects = False
if binary_effects:
    efl = srbe.generate_effect_list(options.effectsfile, verbosity=options.verbose)
    attdat = srbe.read_attributes(options.attributesfile, verbosity=options.verbose)
    effdat = srbe.read_effects(options.effectsfile, efl, verbosity=options.verbose)
else:
    efl = sre.generate_effect_list(options.effectsfile, verbosity=options.verbose)
    (attdat, attnames) = sre.read_attributes(options.attributesfile, verbosity=options.verbose)
    effdat = sre.read_effects(options.effectsfile, efl, verbosity=options.verbose)

# for akey in attdat:
#     print akey.split()[0], akey.split()[0][0:3], akey.split()[0][3:]

##### brief summary of method in paper:
##### find the projection matrices WE(mxk) and WR(pxk) which project R and E into
##### a k-dimensional subspace such that
##### (WE' E) and (WR' R) are maximally correlated, i.e. where
##### R' WR WE' E is maximized

myccmodel = sru.CCModel()
myccmodel = sru.canonical_correlation(efl, attdat, effdat, options.projected_size, verbosity=options.verbose)
pickle.dump(myccmodel, open("ccmodel.pkl","wb"))
