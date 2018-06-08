#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:52:36 2018

@author: BallBlueMeercat
"""

import numpy as np
from datasim import magn

def lnlike(theta, zpicks, mag, sigma, firstderivs_key):
#    print('@@@@ lnlike has been called')
#    
#    print('lnline -- theta', theta)
#    print('type theta', type(theta))
#    print('theta.ndim', theta.ndim)
    
    params = {}
    for i in range(len(theta)):
        if i == 0:
            params['m'] = theta[i]
        elif i == 1:
            params['gamma'] = theta[i]
        elif i == 2:
            params['de'] = theta[i]
    
    model = magn(params, zpicks, firstderivs_key)
    inv_sigma2 = 1.0/(sigma**2)
    return -0.5*(np.sum((mag-model)**2*inv_sigma2 - np.log(inv_sigma2)))

def lnprior(theta):
#    print(' lnprior has been called')   
#    print('lnprior speaking: theta = ',theta)
    
    if len(theta) == 1:
        m = theta
        if 0 < m < 1 or m == 1:
            return 0.0
    elif len(theta) == 2:
        m, gamma = theta
        if (0 < m < 1 or m == 1) and abs(gamma) < 0.1:
            return 0.0
    elif len(theta) == 3:
        m, gamma, de = theta
        if (0 < m < 1 or m == 1) and abs(gamma) < 0.1 and (0 < de < 1):
            return 0.0
        
    return -np.inf

def lnprob(theta, zpicks, mag, sigma, firstderivs_key):
#    print('@@@@@ lnprob has been called')
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    
    return lp + lnlike(theta, zpicks, mag, sigma, firstderivs_key)