#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:52:36 2018

@author: BallBlueMeercat
"""

import numpy as np
from datasim import magn

def lnlike(theta, data, sigma, firstderivs_key, ndim, params):
    '''
    Finding matter density m, absolute M, alpha, beta, interaction gamma.
    '''
    mag = data['mag']
    
    for i in range(ndim):
        dic = params[i]
        for key in dic:
            dic[key] = theta[i]
    
    model = magn(params, data, firstderivs_key)
    var = sigma**2
    return -0.5*np.sum((mag-model)**2 /var +0.5*np.log(2*np.pi*var))

def lnprior(theta, key):
    '''
    Finding matter density m, absolute M, alpha, beta, interaction gamma.
    '''      
    Mcorr_min, Mcorr_max = -20, -18
#    alpha_lim = 5
#    beta_lim = 5
     
    if (0 < theta[0] < 1 or theta[0] == 1):
        if Mcorr_min < theta[1] < Mcorr_max:
#            if abs(theta[2]) < alpha_lim:
#                if abs(theta[3]) < beta_lim:
                    
            if key == 'exotic':
                if -2 < theta[2] < 0.1 and -1.5 < abs(theta[3]) < 3.5:
                    return 0.0
                
            elif key == 'late_intxde':
                if -2 < theta[2] < 0.1:
                    return 0.0
                
            elif key == 'heaviside_late_int':
                if -1.45 < theta[2] < 0.1:
                    return 0.0
                
            elif key == 'late_int':
                if -15 < theta[2] < 0.1:
                    return 0.0
                
            elif key == 'expgamma':
                if -0.1 < theta[2] < 1.5:
                    return 0.0
                
            elif key == 'txgamma':
                if -0.5 < theta[2] < 0.1:
                    return 0.0
                
            elif key == 'zxgamma':
                if -10 < theta[2] < 0.1:
                    return 0.0
                
            elif key == 'zxxgamma':
                if -0.1 < theta[2] < 12:
                    return 0.0
                
            elif key == 'gammaxxz':
                if -1 < theta[2] < 1:
                    return 0.0
                
            elif key == 'rdecay':
                if -2 < theta[2] < 0.1:
                    return 0.0
                
            elif key == 'interacting':
                if -1.5 < theta[2] < 0.1:
                    return 0.0
                
            elif key == 'LCDM':
                return 0.0
            
            else:
                if abs(theta[4]) < 10:
                    return 0.0
        
    return -np.inf

def lnprob(theta, data, sigma, firstderivs_key, ndim, test_params):

    lp = lnprior(theta, firstderivs_key)
    if not np.isfinite(lp):
        return -np.inf
    
    return lp + lnlike(theta, data, sigma, firstderivs_key, ndim, test_params)