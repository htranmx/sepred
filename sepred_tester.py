#!/usr/bin/env python
"""reads in attributes file (descriptors), makes predictions using modelfile (e.g. ccmodel.pkl), and compares with known values (nci60/GI50_RAW.txt). Prints out list of (known, predicted) pairs."""

import sys
import numpy as np
sys.path.append('~/tools/sepred')
import cPickle as pickle
import sepred_reader_effects as sre
from optparse import OptionParser

parser = OptionParser(__doc__)

parser.add_option("-a", "--attributes",
                  dest="attributesfile",
                  help="read in the attribute information for the dataset from ATTRIBUTES",
                  metavar="ATTRIBUTES",
                  default="nci60/nci_compounds_first100_of_2265.csv")

parser.add_option("-m", "--modelfile",
                  dest="modelfile",
                  help="use model in MODELFILE",
                  metavar="MODELFILE",
                  default="ccmodel.pkl")

(options, args) = parser.parse_args()

ccm = pickle.load(open(options.modelfile))
(newattdat, newattnames) = sre.read_attributes(options.attributesfile, verbosity=True)

### read in effects
expdata = 'nci60/GI50_RAW.txt'
efl = sre.generate_effect_list(expdata, verbosity=True)
expeffdat = sre.read_effects(expdata, efl, verbosity=True)
### note: effects are in -log10 format
### to convert back into format in GI50_RAW, do:
# print np.exp(- expeffdat['NSC9858'])

compounds = newattdat.keys()

# print ccm.binary_cut
# raise RuntimeError
# for ic in ['NSC9858']:
for ic in compounds:
    ### actual value
    actual = expeffdat[ic]
    ### predicted value
    predval = ccm.predict(newattdat[ic])
    if len(actual) != len(predval):
        print "Error: inconsistency in array lengths"
        raise RuntimeError
    else:
        for vv in range(len(actual)):
            print "{0:12.3g}{1:12.3g}".format(actual[vv], predval[vv])
            # print "{0:5d}{1:5d}".format(int(actual[vv] >= ccm.binary_cut), int(predval[vv] >= ccm.binary_cut))
