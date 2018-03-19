#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""

import numpy as np

# Choose the "true" parameters.
m_true = -0.9594
b_true = 4.294
f_true = 0.534

# Generate some synthetic data from the model.
N = 50
x = np.sort(10*np.random.rand(N))
yerr = 0.1+0.5*np.random.rand(N)
print('yerr = ',yerr)
y = m_true*x+b_true
print('y = ',y)
y += np.abs(f_true*y) * np.random.randn(N)
print('y += np.abs(f_true*y) * np.random.randn(N) = ',y)
y += yerr * np.random.randn(N)
print('y += yerr * np.random.randn(N)',y)