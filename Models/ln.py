#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:52:36 2018

@author: BallBlueMeercat
"""

import numpy as np
from datasim import magn

def lnlike(theta, zpicks, mag, sigma, firstderivs_key, ndim):
#    print('@@@@ lnlike has been called')
#    
#    print('lnline -- theta', theta)
#    print('type theta', type(theta))
#    print('theta.ndim', theta.ndim)
    
    params = {}
    if ndim == 1:
        params = {'m':theta}
    elif ndim == 2:
        params = {'m':theta[0],'gamma':theta[1]}
    elif ndim == 3:
        params= {'m':theta[0],'gamma':theta[1],'de':theta[2]}
    
    model = magn(params, zpicks, firstderivs_key)
    var = sigma**2
    return -0.5*np.sum((mag-model)**2 /var +0.5*np.log(2*np.pi*var))
#    inv_var = 1.0/(sigma**2)
#    return -0.5*(np.sum((mag-model)**2*inv_var - np.log(inv_var)))

def lnprior(theta, ndim):
#    print(' lnprior has been called')   
#    print('lnprior speaking: theta = ',theta)
#    print('type(theta)',str(type(theta)))
    
    if ndim == 1:
        m = theta
        if 0 < m < 1 or m == 1:
            return 0.0
    elif ndim == 2:
        m, gamma = theta
        if (0 < m < 1 or m == 1) and abs(gamma) < 10:
            return 0.0
    elif ndim == 3:
        m, gamma, de = theta
        if (0 < m < 1 or m == 1) and abs(gamma) < 10 and (0 < de < 1):
            return 0.0
        
    return -np.inf

def lnprob(theta, zpicks, mag, sigma, firstderivs_key, ndim):
#    print('@@@@@ lnprob has been called')
    lp = lnprior(theta, ndim)
    if not np.isfinite(lp):
        return -np.inf
    
    return lp + lnlike(theta, zpicks, mag, sigma, firstderivs_key, ndim)