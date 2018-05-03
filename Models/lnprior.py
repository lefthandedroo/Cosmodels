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
#    print('lnprior speaking: theta = ',theta)
    m = theta
    if 0 < m < 1 or m == 1:
        return 0.0
    return -np.inf