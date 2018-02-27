#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:33:56 2018

@author: BallBlueMeercat
"""
import numpy as np

def gnoise(mag, mu, sigma, n):
    """
   Returns:
       mag = mag, each point offset by unique Gaussian noise;
       noise = Gaussian noise.
    """
    print('-gnoise has been called')
    noise = np.random.normal(mu,sigma,n)
    print(type(noise))
#    print('noise from inside gnoise is = ', noise)
    mag = mag + noise
    return mag, noise