#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:52:36 2018

@author: BallBlueMeercat
"""

import lnlike
import lnprior
import numpy as np

def lnprob(theta, zpicks, mag, sigma):
#    print('@@@@@ lnprob has been called')
    lp = lnprior.lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    
    return lp + lnlike.lnlike(theta, zpicks, mag, sigma)
