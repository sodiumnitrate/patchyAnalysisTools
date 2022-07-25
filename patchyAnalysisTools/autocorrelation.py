import numpy as np
from . import utils
from . import trajectory
import pdb

def autocorrelation(start,n_frames,A):
    '''
    Calculates autocorrelation for an observable A
    '''

    mean_A = np.mean(A[start:])
    autocorr = np.zeros(n_frames-start)
    for f in range(1,n_frames-start,1):
        integrand = 0
        count = 0
        # loop over all possible starting points
        for ff in range(start+1,n_frames-f):
            count += 1
            # calculate delta A at the starting point
            dAt = A[ff] - mean_A

            # 
            dAtt = A[ff:ff+f+1] - mean_A

            # sum over the product
            integrand += np.sum(dAt*dAtt)

        if count != 0:
            # normalize so that we have a proper average
            integrand /= count * f
            autocorr[f] = integrand

    assert(autocorr[1] != 0)
    return autocorr[1:] / autocorr[1]