import inspect
import numpy as np
from scipy import optimize


class CCModel:
    effectlist = []
    Rmean = np.array((), dtype='float')
    Emean = np.array((), dtype='float')
    WE = np.array((), dtype='float')
    WR = np.array((), dtype='float')
    CEE = np.array((), dtype='float')
    num_attributes = 0
    num_effects = 0
    binary_cut = 0.5
    def predict(self, attribute_vector, binary=False):
        """Input: unnormalized attribute vector.
        Output: unnormalized effect vector."""
        if len(attribute_vector) != self.num_attributes:
            print "Error: CCModel: Number of attributes ({0}) does not equal {1}".format(len(attribute_vector), self.num_attributes)
            raise RuntimeError
        else:
            norm_att_vec = attribute_vector - self.Rmean
            predef = np.dot(self.CEE, np.dot(self.WE, np.dot(self.WR.transpose(), norm_att_vec)))
        if binary:
            binary_effects = np.array(((predef + self.Emean) >= self.binary_cut), dtype='int')
            return binary_effects
        else:
            return predef + self.Emean
    def binary_to_effects(self, beff):
        if len(self.effectlist) != len(beff):
            print "CCModel: inconsistent effect list: {0}, {1}".format(len(self.effectlist), len(beff))
            raise RuntimeError
        else:
            retval = []
            for ii in np.nonzero(beff)[0]:
                retval.append(self.effectlist[ii])
            return retval
        


def canonical_correlation(effectlist, attribute_data, effect_data, projected_size, verbosity=False):
    retccmodel = CCModel()
    ### create list of items which are in both attribute_data and effect_data
    itemlist = []
    for it in  effect_data:
        if not (it in attribute_data):
            # print "{0} not found".format(it)
            continue
        else:
            itemlist.append(it)
    if verbosity:
        print "Found {0} items present in both attribute and effect data".format(len(itemlist))

    ##### create attribute matrix R(pxn) and effect E(mxn) matrices from the input data
    n = len(itemlist)
    p = len(attribute_data[attribute_data.keys()[0]])
    m = len(effectlist)
    
    R = np.zeros((p, n))
    E = np.zeros((m, n))
    for i in range(n):
        thisitem = itemlist[i]
        R[:,i] = attribute_data[thisitem]
        E[:,i] = effect_data[thisitem]
        if np.isnan(np.sum(R[:,i])):
            print "Error: {0} from {1}:".format(str(inspect.stack()[0][3]), str(inspect.stack()[1][3]))
            print "  Found not a number: (key {0}) {1}".format(thisitem, attribute_data[thisitem])
            raise RuntimeError

    #### normalize rows of E,R to have mean of 0
    Rmean = np.mean(R, axis=1)
    Emean = np.mean(E, axis=1)
    for i in range(n):
        R[:,i] = R[:,i] - Rmean
        E[:,i] = E[:,i] - Emean

    CEE = np.dot(E, E.transpose())
    CRR = np.dot(R, R.transpose())
    CRE = np.dot(R, E.transpose())
    CER = np.dot(E, R.transpose())

    try:
        CEEinv = np.linalg.inv(CEE)
    except np.linalg.LinAlgError:
        CEEinv = np.linalg.pinv(CEE)

    try:
        CRRinv = np.linalg.inv(CRR)
    except np.linalg.LinAlgError:
        CRRinv = np.linalg.pinv(CRR)

    cccc = np.dot(CEEinv, np.dot(CER, np.dot(CRRinv, CRE)))
    (lsq, we) = np.linalg.eig(cccc)
    isortlsq = np.argsort(-np.real(lsq))

    ### the subspace is of size 'k'
    k = int(projected_size)
    if k < 1:
        # estimate a value for k
        k = len(np.nonzero(np.real(lsq) > 0.25 * np.real(lsq[isortlsq[0]]))[0])
        # k = we.shape[0]
        if verbosity:
            print "k was set to {0}".format(k)

    if verbosity: print "eigenvalues:\n", np.real(lsq[isortlsq[0:k]])
    # print we.shape
    WE = np.zeros((m, k))
    WR = np.zeros((p, k))
    for i in range(k):
        thislsq = np.real(lsq[isortlsq[i]])
        thisl = np.sqrt(thislsq)
        thiswe = np.real(we[:,isortlsq[i]])
        WE[:,i] = thiswe
        WR[:,i] = np.dot(CRRinv, np.dot(CRE, thiswe)) / thisl

    if verbosity:
        errormat = np.dot(np.dot(WR.transpose(), R).transpose(),  np.dot(WE.transpose(), E))
        # print errormat.shape
        etrace = 0.0
        for i in range(errormat.shape[0]):
            etrace += errormat[i,i]
        print "Trace = {0}".format(etrace)

    #### populate CCModel
    retccmodel.effectlist = effectlist
    retccmodel.Rmean = Rmean
    retccmodel.Emean = Emean
    retccmodel.WE = WE
    retccmodel.WR = WR
    retccmodel.CEE = CEE
    retccmodel.num_attributes = p
    retccmodel.num_effects = m
    ### calculate optimum cutoff
    mean_effects = []
    for ekey in effect_data:
        mean_effects.append(np.mean(effect_data[ekey]))
    opterr = lambda x:  - ccmodel_gmean(attribute_data, effect_data, retccmodel, x)
    optcut = optimize.fmin(opterr, np.mean(mean_effects), xtol = 0.01)
    # optcut = optimize.fmin(opterr, -10.0, xtol = 0.01)
    print "Optimum cutoff: {0}".format(optcut[0])
    retccmodel.binary_cut = optcut[0]
    return retccmodel


def ccmodel_gmean(newattdat, effdat, myccmodel, cutoff):
    """newattdat --> attributes to be used for prediction
    effdat --> the correct effects
    cutoff"""
    if len(cutoff) == 1:
        cutoff = cutoff[0]
    npred = 0
    sumofpred = 0.0
    # cutoffs = np.arange(0, 5, 0.1)
    truepos = 0
    falspos = 0
    falsneg = 0
    trueneg = 0
    for item in newattdat:
        thisav = newattdat[item]
        epred = myccmodel.predict(thisav)
        epred_normalized = epred - myccmodel.Emean
        if item in effdat:
            npred += 1
            thiseff = effdat[item]
            ### compare predicted (epred) and actual (thiseff)
            # print myccmodel.Emean
            # print "Prediction for {0}".format(item)
            # print "   predicted: {0}".format(epred)
            # print "   actual:    {0}".format(thiseff)
            # for ik in range(len(cutoffs)):
            thiscut = cutoff
            dtruepos = len(np.nonzero(epred[np.nonzero(thiseff)[0]] >= thiscut)[0])
            dfalspos = len(np.nonzero(epred[np.nonzero(thiseff == 0)[0]] >= thiscut)[0])
            dfalsneg = len(np.nonzero(epred[np.nonzero(thiseff)[0]] < thiscut)[0])
            dtrueneg = len(np.nonzero(epred[np.nonzero(thiseff == 0)[0]] < thiscut)[0])
            if dtruepos + dfalspos + dfalsneg + dtrueneg != len(thiseff):
                print "Confusion matrix does not sum to total"
                raise RuntimeError
            truepos += dtruepos
            falspos += dfalspos
            falsneg += dfalsneg
            trueneg += dtrueneg            
            # thiseff_normalized = thiseff - myccmodel.Emean
            # sumofpred += np.dot(epred.transpose(), thiseff_normalized) / (np.linalg.norm(epred) * np.linalg.norm(thiseff_normalized))
            # print "---> {0} <---".format(item)
            # print np.mean(effdat[item]), np.std(effdat[item])
            # print np.mean(epred), np.std(epred)
    recall = float(truepos) / (falsneg + truepos)
    precis = float(truepos) / (falspos + truepos)
    gmean = np.sqrt(recall * precis)
    print "truepos {0}".format(truepos)
    print "falspos {0}".format(falspos)
    print "trueneg {0}".format(trueneg)
    print "falsneg {0}".format(falsneg)
    # for ii in range(len(cutoffs)):
    # print cutoff, gmean, recall, precis
    print "Cutoff = {0:6.2f}, gmean = {1:8.3f}, recall = {2:6.2f}, precision = {3:6.2f}".format(cutoff, gmean, recall, precis)
    # print truepos[ii] + falspos[ii] + falsneg[ii] + trueneg[ii]
    return gmean
    # return np.sqrt(recall * precis * precis)
