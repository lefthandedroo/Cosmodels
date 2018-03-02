#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:51:18 2018

@author: BallBlueMeercat

put constraints on the parameters for msim.msim
"""

import msim
import numpy as np
#import lnprior

def lnlike(theta, zpicks, mag, sigma):
    gamma, m, de = theta
#    if not (-0.1 < gamma < 0.1 and 0.299 < m < 0.301 and 0.699 < de < 0.701):
#        print('@@@@ lnlike has been called')
#        print('theta = ',theta)
#        return -np.inf
    
#    lp = lnprior.lnprior(theta)
#    if not np.isfinite(lp):
#        return -np.inf
    
    model = msim.msim(gamma, m, de, zpicks)
    inv_sigma2 = 1.0/(sigma**2)
    return -0.5*(np.sum((mag-model)**2*inv_sigma2 - np.log(inv_sigma2)))