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
    if -0.1 < gamma < 0.1 and 0.299 < m < 0.301 and 0.699 < de < 0.701:
#        print('lnprior good theta ',theta)
        return 0.0
#    print('lnprior bad theta ',theta)
    return -np.inf