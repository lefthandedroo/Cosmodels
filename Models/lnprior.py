#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:54:17 2018

@author: BallBlueMeercat
"""
import numpy as np

def lnprior(theta):
#    print(' lnprior has been called')
    lamb, m, de = theta
    if (-0.001 < lamb < 0.001 and 0 < m < 1.0 and 0.0 < de < 1.0):
        return 0.0
    return -np.inf