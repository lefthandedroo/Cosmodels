#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:54:17 2018

@author: BallBlueMeercat
"""
import numpy as np

def lnprior(theta):
#    print(' lnprior has been called')
    gamma, m, de = theta
#    if -0.1 < gamma < 0.1 and 0.299 < m < 0.301 and 0.699 < de < 0.701:
    flat = m + de
#    if -3 < gamma < 2.5 and -3.5 < m < 3 and -2 < de < 4.5 and flat == 1:
    if -3 < gamma < 2.5 and -3.5 < m < 3 and -2 < de < 4.5 and 0.5 < flat < 1.5:

        return 0.0
    return -np.inf