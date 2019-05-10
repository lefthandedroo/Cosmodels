#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:52:36 2018

@author: BallBlueMeercat
"""

import numpy as np
from datasim import magn

def lnlike(theta, data, sigma, firstderivs_key, names):
    '''
    Takes in:
        theta = numpy.ndarray, guessed values for parameters;
        data = dictionary, {'mag':mag, 'zpicks':zpicks};
        sigma = int/float, error on data;
        firstderivs_key = string, model being tested;
        names = list of strings.
    Returns:
        float, likelihood for firstderivs_key model with theta parameters.
    '''
    mag = data['mag']
    model_mag, model_da = magn(names, theta, data, firstderivs_key) # mag, but with theta params
    var = sigma**2
    likelihood = -0.5*np.sum((mag-model_mag)**2.0 /var +np.log(2.0*np.pi*var))
    return likelihood

def lnprior(th, key, int_in):
    '''
    Takes in:
        th = numpy.ndarray, guessed values for parameters;
        key = string, model being tested;
        int_in = index of interaction terms, integer or None.
    Returns:
        0.0 if all conditions on theta values are met;
        -np.inf if theta values are outside of prior.
    '''
    if key == 'rainbow' or key == 'niagara' or key == 'kanangra' or 'waterfall' or key == 'stepfall' or key == 'exotic' or key == 'gammaxxz':
        int_lim_min, int_lim_max = -1, 1
    elif key == 'late_intxde':
        int_lim_min, int_lim_max = -2, 0.1
    elif key == 'heaviside_late_int':
        int_lim_min, int_lim_max = -1.45, 0.1
    elif key == 'late_int':
        int_lim_min, int_lim_max = -15, 0.1
    elif key == 'expgamma':
        int_lim_min, int_lim_max = -0.1, 1.5
    elif key == 'txgamma':
        int_lim_min, int_lim_max = -0.5, 0.1
    elif key == 'zxgamma':
        int_lim_min, int_lim_max = -10, 0.1
    elif key == 'zxxgamma':
        int_lim_min, int_lim_max = -0.1, 12
    elif key == 'rdecay':
        int_lim_min, int_lim_max = -10, 1
    elif key == 'interacting':
        int_lim_min, int_lim_max = -1.5, 0.1
    else:
        int_lim_min, int_lim_max = -10, 10

    Mcorr_min, Mcorr_max = -20, -18
    assert len(th) > 0, f'len(th) = 0, theta = {th}'

    fluids = th[1:int_in]
    int_terms = th[int_in:]

    if Mcorr_min < th[0] < Mcorr_max:
        if max(fluids) == 1 or  max(fluids) < 1:
            if min(fluids) == 0 or min(fluids) > 0:
                if len(int_terms) == 0:
                    return 0.0
                elif min(int_terms) > int_lim_min:
                    if max(int_terms) < int_lim_max:
                        return 0.0
    return -np.inf


def lnprob(theta, data, sigma, firstderivs_key, names, int_in):
    lp = lnprior(theta, firstderivs_key, int_in)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, data, sigma, firstderivs_key, names)