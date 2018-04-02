#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:51:18 2018

@author: BallBlueMeercat

put constraints on the parameters for msim.msim
"""

import msim
import lnprior
import numpy as np

def lnlike(theta, zpicks, mag, noise):
    gamma, m, de = theta
#    print('@@@@ lnlike has been called')
    
    lp = lnprior.lnprior(theta)
    if not np.isfinite(lp):
        print('bad theta passed to msim from lnlike')
        print('theta = ',theta)

    model = msim.msim(gamma, m, de, zpicks)
    inv_sigma2 = 1.0/(noise**2)
    return -0.5*(np.sum((mag-model)**2*inv_sigma2 - np.log(inv_sigma2)))