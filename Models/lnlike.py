#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:51:18 2018

@author: BallBlueMeercat
"""

import msim
import numpy as np

def lnlike(theta, n, p, zpicks, mag, sigma):
#    print('@@@@ lnlike has been called')
    lamb, m, de = theta
    model = msim.msim(lamb, m, de, n, p, zpicks)
    inv_sigma2 = 1.0/(sigma**2)
    return -0.5*(np.sum((mag-model)**2*inv_sigma2 - np.log(inv_sigma2)))