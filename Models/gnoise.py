#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:33:56 2018

@author: BallBlueMeercat
"""
import numpy as np

def gnoise(mag, sigma, mu, n):
    """
    calculates and adds random p% Gaussian noise to each mag datapoint
    """
#    print('-gnoise has been called')
    noise = np.random.normal(mu,sigma,n)
#    print('noise from inside gnoise is = ', noise)
    mag = mag + noise
    return mag, noise