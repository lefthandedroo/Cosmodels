#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:54:17 2018

@author: BallBlueMeercat
"""
import numpy as np

gamma_prior = 0.1

def lnprior(theta):
#    print(' lnprior has been called')
    gamma, m = theta
#    if -3 < gamma < 2.5 and -3.5 < m < 3 and -2 < de < 4.5:
    if -gamma_prior < gamma < gamma_prior and (0 < m < 1.1 or m == 0):
        return 0.0
    return -np.inf