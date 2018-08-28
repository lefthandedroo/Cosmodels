#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:52:36 2018

@author: BallBlueMeercat
"""

import numpy as np
from datasim import magn

#def lnlike(theta, data, sigma, firstderivs_key, ndim):
#    
##    print('type data',type(data))
#    mag = data['mag']
##    print('type mag',type(mag))
#    
#    params = {}
#    if ndim == 1:
#        params = {'m':theta}
#    elif ndim == 2:
#        params = {'m':theta[0],'gamma':theta[1]}
#    elif ndim == 3:
#        params= {'m':theta[0],'gamma':theta[1],'alpha':theta[2], 'beta':theta[3]}
#    
#    model = magn(params, data, firstderivs_key)
#    var = sigma**2
#    return -0.5*np.sum((mag-model)**2 /var +0.5*np.log(2*np.pi*var))

def lnlike(theta, data, sigma, firstderivs_key, ndim):
    
#    print('type data',type(data))
    mag = data['mag']
#    print('type mag',type(mag))
        
    params = {}
    if ndim == 1:
        params = {'m':theta}
    elif ndim == 2:
        params = {'m':theta[0],'gamma':theta[1]}
    elif ndim == 3:
        params= {'m':theta[0],'alpha':theta[1], 'beta':theta[2]}
    elif ndim == 4:
        params= {'m':theta[0],'gamma':theta[1],'alpha':theta[2], 'beta':theta[3]}
    
    model = magn(params, data, firstderivs_key)
    var = sigma**2
    return -0.5*np.sum((mag-model)**2 /var +0.5*np.log(2*np.pi*var))
    
#def lnlike(theta, zpicks, mag, sigma, firstderivs_key, ndim):
#    
#    params = {}
#    if ndim == 1:
#        params = {'m':theta}
#    elif ndim == 2:
#        if firstderivs_key == 'LCDM':
#            params = {'m':theta[0],'M':theta[1]}
#        else: 
#            params = {'m':theta[0],'gamma':theta[1]}
#    elif ndim == 3:
#        params= {'m':theta[0],'gamma':theta[1],'M':theta[2]}
#    
#    if params['M']:
#        mag -= params.get('M')
#    model = magn(params, zpicks, firstderivs_key)
#    var = sigma**2
#    return -0.5*np.sum((mag-model)**2 /var +0.5*np.log(2*np.pi*var))

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


#def lnprior(theta, key):
#    
#    if key == 'LCDM':
#        m = theta
#        if 0 < m < 1 or m == 1:
#            return 0.0
#    elif key == 'late_int' or 'heaviside_late_int' or 'late_intxde':
#        m, gamma = theta
#        if (0 < m < 1 or m == 1) and -1.45 < gamma < 0.2:
#            return 0.0       
#    elif key == 'rdecay':
#        m, gamma = theta
#        if (0 < m < 1 or m == 1) and -10 < gamma < 0:
#            return 0.0
#    elif key == 'interacting':
#        m, gamma = theta
#        if (0 < m < 1 or m == 1) and abs(gamma) < 1.45:
#            return 0.0
#    elif key == 'expgamma':
#        m, gamma = theta
#        if (0 < m < 1 or m == 1) and abs(gamma) < 25:
#            return 0.0
#    elif key == 'zxxgamma' or 'gammaxxz':
#        m, gamma = theta
#        if (0 < m < 1 or m == 1) and 0 < gamma < 10:
#            return 0.0        
#    else:
#        m, gamma = theta
#        if (0 < m < 1 or m == 1) and abs(gamma) < 10:
#            return 0.0
#        
#    return -np.inf
    
def lnprior(theta, key):
    
    if key == 'LCDM':
        m, alpha, beta = theta
        if 0 < m < 1 or m == 1 and -10 < alpha < 10 and -10 < beta < 10:
            return 0.0
    elif key == 'late_int' or 'heaviside_late_int' or 'late_intxde':
        m, gamma, alpha, beta = theta
        if (0 < m < 1 or m == 1) and -1.45 < gamma < 0.2 and -10 < alpha < 10 and -10 < beta < 10:
            return 0.0       
    elif key == 'rdecay':
        m, gamma, alpha, beta = theta
        if (0 < m < 1 or m == 1) and -10 < gamma < 0 and -10 < alpha < 10 and -10 < beta < 10:
            return 0.0
    elif key == 'interacting':
        m, gamma, alpha, beta = theta
        if (0 < m < 1 or m == 1) and abs(gamma) < 1.45 and -10 < alpha < 10 and -10 < beta < 10:
            return 0.0
    elif key == 'expgamma':
        m, gamma, alpha, beta = theta
        if (0 < m < 1 or m == 1) and abs(gamma) < 25 and -10 < alpha < 10 and -10 < beta < 10:
            return 0.0
    elif key == 'zxxgamma' or 'gammaxxz':
        m, gamma, alpha, beta = theta
        if (0 < m < 1 or m == 1) and 0 < gamma < 10 and -10 < alpha < 10 and -10 < beta < 10:
            return 0.0        
    else:
        m, gamma, alpha, beta = theta
        if (0 < m < 1 or m == 1) and abs(gamma) < 10 and -10 < alpha < 10 and -10 < beta < 10:
            return 0.0
        
    return -np.inf


def lnprob(theta, data, sigma, firstderivs_key, ndim):

    lp = lnprior(theta, firstderivs_key)
    if not np.isfinite(lp):
        return -np.inf
    
    return lp + lnlike(theta, data, sigma, firstderivs_key, ndim)