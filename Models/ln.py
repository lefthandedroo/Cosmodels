#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:52:36 2018

@author: BallBlueMeercat
"""

import numpy as np
import datasim

def lnlike(theta, zpicks, mag, sigma):
#    print('@@@@ lnlike has been called')
    model = datasim.mag(theta, zpicks)
    inv_sigma2 = 1.0/(sigma**2)
    return -0.5*(np.sum((mag-model)**2*inv_sigma2 - np.log(inv_sigma2)))

def lnprior(theta):
#    print(' lnprior has been called')   
#    print('lnprior speaking: theta = ',theta)
    m = theta
    if 0 < m < 1 or m == 1:
        return 0.0
    return -np.inf

def lnprob(theta, zpicks, mag, sigma):
#    print('@@@@@ lnprob has been called')
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    
    return lp + lnlike(theta, zpicks, mag, sigma)