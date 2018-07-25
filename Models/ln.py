#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:52:36 2018

@author: BallBlueMeercat
"""

import numpy as np
from datasim import magn

def lnlike(theta, zpicks, mag, sigma, firstderivs_key, ndim):
    
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

#def lnprior(theta, ndim):
##    print(' lnprior has been called')   
##    print('lnprior speaking: theta = ',theta)
##    print('type(theta)',str(type(theta)))
#    
#    if ndim == 1:
#        m = theta
#        if 0 < m < 1 or m == 1:
#            return 0.0
#    elif ndim == 2:
#        m, gamma = theta
#        if (0 < m < 1 or m == 1) and -10 < gamma < 10:
#            return 0.0
#    elif ndim == 3:
#        m, gamma, de = theta
#        if (0 < m < 1 or m == 1) and abs(gamma) < 100 and (0 < de < 1):
#            return 0.0
#        
#    return -np.inf


def lnprior(theta, firstderivs_key):
    
    if firstderivs_key == 'LCDM':
        m = theta
        if 0 < m < 1 or m == 1:
            return 0.0
    elif firstderivs_key == 'late_int':
        m, gamma = theta
        if (0 < m < 1 or m == 1) and -1.45 < gamma < 0.2:
            return 0.0       
    elif firstderivs_key == 'rdecay':
        m, gamma = theta
        if (0 < m < 1 or m == 1) and -10 < gamma < 0:
            return 0.0
    elif firstderivs_key == 'interacting':
        m, gamma = theta
        if (0 < m < 1 or m == 1) and abs(gamma) < 1.45:
            return 0.0
    elif firstderivs_key == 'expgamma':
        m, gamma = theta
        if (0 < m < 1 or m == 1) and abs(gamma) < 25:
            return 0.0
    elif firstderivs_key == 'zxxgamma' or firstderivs_key == 'gammaxxz':
        m, gamma = theta
        if (0 < m < 1 or m == 1) and 0 < gamma:
            return 0.0        
    else:
        m, gamma = theta
        if (0 < m < 1 or m == 1) and abs(gamma) < 10:
            return 0.0
        
    return -np.inf


def lnprob(theta, zpicks, mag, sigma, firstderivs_key, ndim):

    lp = lnprior(theta, firstderivs_key)
    if not np.isfinite(lp):
        return -np.inf
    
    return lp + lnlike(theta, zpicks, mag, sigma, firstderivs_key, ndim)