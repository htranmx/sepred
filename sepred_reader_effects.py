import numpy as np
import csv

# def read_effect_list(filename, verbosity=False):
#     """returns sorted list of side effects found in <filename>"""
#     with open(filename) as el:
#         ellines = el.readlines()
#     efl = set()
#     for ell in ellines:
#         if ell.strip():
#             efl.add(ell.strip())
#     efl = list(efl)
#     efl.sort()
#     if verbosity:
#         print "possible effects: {0}, from {1}".format(len(efl), filename)
#     return efl

def csvline_to_list(csvline):
    """converts comma separated, double-quoted fields into a list"""
    csvr = csv.reader([csvline], delimiter=',', quotechar='"')
    tmp = []
    for cl in csvr:
        tmp.append(cl)
    if len(tmp) == 1:
        return tmp[0]
    else:
        print "Error: csvline_to_list: found {0} elements (should have only 1)".format(len(tmp))
        raise RuntimeError

def read_attributes(filename, verbosity=False):
    """read attribues for each item in the dataset from <filename>"""
    attdat = {}
    nattr = -1
    with open(filename) as af:
        aflines = af.readlines()
    ### first line contains field names
    ### first field is drug name, 2nd is smiles string
    attnames = csvline_to_list(aflines[0])
    attnames = attnames[2:]
    nattr = len(attnames)
    if verbosity:
        print "Found {0} attributes:".format(nattr)
        print attnames
    for iafl in range(1, len(aflines)):
        afl = aflines[iafl]
        afls = csvline_to_list(afl)
        dbname = afls[0]
        if len(dbname.split()) > 1:
            dbname = dbname.split()[0]
        if dbname:
            if nattr != len(afls[2:]):
                print "read_attributes: Number of attributes not consistent {0}, {1}, attributes:".format(nattr, len(afls[2:]))
                print afls[2:]
                raise RuntimeError()
            else:
                tmpattdat = np.array(afls[2:], dtype='float')
                for itmp in range(len(tmpattdat)):
                    if not np.isfinite(tmpattdat[itmp]):
                        # print "setting {0}[{1}] to 0.0".format(dbname, itmp)
                        tmpattdat[itmp] = 0.0
                attdat[dbname] = tmpattdat
    if verbosity:
        print "(read_attributes): attribute dataset: {1} attributes for {0} items with drugbank ID, from {2}".format(len(attdat), nattr, filename)
    return attdat, attnames

# def effects_to_binary(search_effects, effectlist, verbosity=False):
#     """returns a binary array of length len(effectlist) where occurences of elements in
#     search_effects are 1, and others 0"""
#     retval = np.zeros(len(effectlist), dtype='int')
#     for ie in range(len(effectlist)):
#         if effectlist[ie] in search_effects:
#             retval[ie] = 1
#     return retval

# def binary_to_effects(binary_mask, effectlist, verbosity=False):
#     """returns a list of effects, given the binary mask generated by effects_to_binary_mask"""
#     retval = []
#     for ie in np.nonzero(binary_mask)[0]:
#         retval.append(effectlist[ie])
#     retval = list(set(retval))
#     retval.sort()
#     return retval

def read_effects(filename, effectlist, verbosity=False):
    """read effects for each item in dataset from <filename> and generates a mask by comparing
    to the list of effects <efl> (generated by read_effect_list).
    Returns effdat. To get the list of effects for <item>, do call
    binary_to_effects(effdat[item], effectlist)."""
    effdat = {}
    with open(filename) as ef:
        eflines = ef.readlines()
    #### one header line
    for iefl in range(1,len(eflines)):
        efl = eflines[iefl].strip()
        efls = efl.split('\t')
        nscno = efls[0]
        itemname = efls[1]
        dbname = 'NSC' + nscno.strip()
        effect = efls[2:]
        for ieit in range(len(effect)):
            try:
                dummy = float(effect[ieit])
            except ValueError:
                #### GI50 values of NA (resistant), set to 0.999 (a large number)
                effect[ieit] = '0.999'
        effect = -np.log10(np.array(effect, dtype='float'))
        if len(effect) != len(effectlist):
            print "Error: read_effects: found {0} effects, expecting {1}".format(len(effect), len(effectlist))
            print effect
            raise RuntimeError
        if dbname in effdat:
            print "Error: read_effects: {0} appears twice in {1}".format(dbname, filename)
            raise RuntimeError
        else:
            effdat[dbname] = effect
    # for item in effdat:
    #     effdat[item] = effects_to_binary(effdat[item], effectlist)
    if verbosity:
        print "Found {0} items with effects, in {1}".format(len(effdat), filename)
    return effdat

def generate_effect_list(filename, verbosity=False):
    """generates a list of unique effects found in the dataset"""
    efl = []
    with open(filename) as ef:
        eflines = ef.readlines()
    # first line has list of drug targets
    efl = eflines[0].split()
    if verbosity:
        print "{0} effects found in dataset from from {1}".format(len(efl), filename)
    return efl