#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""

import numpy as np


m = [-22.68201682720258,
 -32.05185120861957,
 -27.30200425753162,
 -24.132781241304805,
 -23.650626679759398,
 -27.434931605155437,
 -23.763156312261422,
 -26.602708576455566,
 -28.393389697398934,
 -24.537029076423305]

#m = np.array(m)

n = 10
p = 10

m_mu = abs(np.sum(m)/n)         # mean of simulated apparent magnitudes m
sigma = m_mu /2 /100 * p       # standard deviation
mu = 0                       # mean of noise distribution

noise = np.random.normal(mu,sigma,n)
print('noise',noise)
m = m + noise
print('m',m)
