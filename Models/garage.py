#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""

import numpy as np


def lnlike(theta, x, y, sigma):
        a, b = theta
        model = a * x + b
        inv_sigma2 = 1.0/(sigma**2)
        return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2))) 
    
a = 1
b = 2

theta = a, b

N = 10
mu = 0
sigma= 0.5

x = np.random.rand(N)*4                 # picking random points on x-axis
yerr = np.random.normal(mu,sigma,N)     # Gaussian noise
y = a * x + b 
y += yerr                               # data, offset in y with noise

outcome = lnlike(theta, x, y, sigma)