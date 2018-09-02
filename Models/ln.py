#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:52:36 2018

@author: BallBlueMeercat
"""

import numpy as np
from datasim import magn

#def lnlike(theta, data, sigma, firstderivs_key, ndim):
#    '''
#    Finding matter density m, interaction gamma.
#    '''
#    mag = data['mag']
#    
#    params = {}
#    if ndim == 1:
#        params = {'m':theta}
#    elif ndim == 2:
#        params = {'m':theta[0],'gamma':theta[1]}
#    
#    model = magn(params, data, firstderivs_key)
#    var = sigma**2
#    return -0.5*np.sum((mag-model)**2 /var +0.5*np.log(2*np.pi*var))


def lnlike(theta, data, sigma, firstderivs_key, ndim):
    '''
    Finding matter density m, corrected absolute mag M, interaction gamma.
    '''    
    mag = data['mag']
    
    params = {}
    if ndim == 2:
        params = {'m':theta[0], 'M':theta[1]}
    elif ndim == 3:
        params = {'m':theta[0],'M':theta[1], 'gamma':theta[2]}
    
    model = magn(params, data, firstderivs_key)
    var = sigma**2
    return -0.5*np.sum((mag-model)**2 /var +0.5*np.log(2*np.pi*var))


#def lnlike(theta, data, sigma, firstderivs_key, ndim):
#    '''
#    Finding matter density m, alpha, beta, interaction gamma.
#    '''
#    mag = data['mag']
#        
#    params = {}
#    if ndim == 1:
#        params = {'m':theta}
#    elif ndim == 2:
#        params = {'m':theta[0],'gamma':theta[1]}
#    elif ndim == 3:
#        params= {'m':theta[0],'alpha':theta[1], 'beta':theta[2]}
#    elif ndim == 4:
#        params= {'m':theta[0],'gamma':theta[1],'alpha':theta[2], 'beta':theta[3]}
#    
#    model = magn(params, data, firstderivs_key)
#    var = sigma**2
#    return -0.5*np.sum((mag-model)**2 /var +0.5*np.log(2*np.pi*var))










#def lnprior(theta, key):
#    '''
#    Finding matter density m, interaction gamma.
#    '''
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
    '''
    Finding matter density m, corrected absolute mag M, interaction gamma.
    '''  
    
    Mmin = -20
    
    Mmax = -18
    
    if key == 'LCDM':
        m, M = theta
        if (0 < m < 1 or m == 1) and Mmin < M < Mmax:
            return 0.0
    elif key == 'late_int' or 'heaviside_late_int' or 'late_intxde':
        m, M, gamma = theta
        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and -1.45 < gamma < 0.2:
            return 0.0       
    elif key == 'rdecay':
        m, M, gamma = theta
        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and -10 < gamma < 0 :
            return 0.0
    elif key == 'interacting':
        m, M, gamma = theta
        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(gamma) < 1.45:
            return 0.0
    elif key == 'expgamma':
        m, M, gamma = theta
        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(gamma) < 25 :
            return 0.0
    elif key == 'zxxgamma' or 'gammaxxz':
        m, M, gamma = theta
        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and 0 < gamma < 10:
            return 0.0        
    else:
        m, M, gamma = theta
        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(gamma) < 10:
            return 0.0
        
    return -np.inf
    

#def lnprior(theta, key):
#    '''
#    Finding matter density m, alpha, beta, interaction gamma.
#    '''  
#    
#    if key == 'LCDM':
#        m, alpha, beta = theta
#        if 0 < m < 1 or m == 1 and -10 < alpha < 10 and -10 < beta < 10:
#            return 0.0
#    elif key == 'late_int' or 'heaviside_late_int' or 'late_intxde':
#        m, gamma, alpha, beta = theta
#        if (0 < m < 1 or m == 1) and -1.45 < gamma < 0.2 and -10 < alpha < 10 and -10 < beta < 10:
#            return 0.0       
#    elif key == 'rdecay':
#        m, gamma, alpha, beta = theta
#        if (0 < m < 1 or m == 1) and -10 < gamma < 0 and -10 < alpha < 10 and -10 < beta < 10:
#            return 0.0
#    elif key == 'interacting':
#        m, gamma, alpha, beta = theta
#        if (0 < m < 1 or m == 1) and abs(gamma) < 1.45 and -10 < alpha < 10 and -10 < beta < 10:
#            return 0.0
#    elif key == 'expgamma':
#        m, gamma, alpha, beta = theta
#        if (0 < m < 1 or m == 1) and abs(gamma) < 25 and -10 < alpha < 10 and -10 < beta < 10:
#            return 0.0
#    elif key == 'zxxgamma' or 'gammaxxz':
#        m, gamma, alpha, beta = theta
#        if (0 < m < 1 or m == 1) and 0 < gamma < 10 and -10 < alpha < 10 and -10 < beta < 10:
#            return 0.0        
#    else:
#        m, gamma, alpha, beta = theta
#        if (0 < m < 1 or m == 1) and abs(gamma) < 10 and -10 < alpha < 10 and -10 < beta < 10:
#            return 0.0
#        
#    return -np.inf


def lnprob(theta, data, sigma, firstderivs_key, ndim):

    lp = lnprior(theta, firstderivs_key)
    if not np.isfinite(lp):
        return -np.inf
    
    return lp + lnlike(theta, data, sigma, firstderivs_key, ndim)