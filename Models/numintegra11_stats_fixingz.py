#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main code
"""
import numpy as np
import time

import timer
import gnoise
import msim
import zpicks
import stats

# Script timer.
timet0 = time.time()

# Parameters:

# Model specific parameters.  
gamma_true = 0  # Interaction term, rate at which DE decays into matter.
m_true = 0.3    # (= e_m(t)/e_crit(t0) at t=t0).
de_true = 0.7   # (de = e_de(t)/e_crit(t0) at t=t0).

# Number of datapoints to be simulated.
n = 1000 #10, 1000

# Statistical parameters:
mu = 0          # mean
sigma = 0.1     # standard deviation

# Code:

# Picking redshifts to investigate.
zmin, zmax = 0.001, 2
zpicks = zpicks.zpicks(zmin, zmax, n)

# Generating apparent magnitues mag at redshift z < zmax (calculated from
# luminosity distances given by LambdaCMD with parameters stated above.
model = msim.msim(gamma_true, m_true, de_true, zpicks)
model = np.asarray(model)
mag = gnoise.gnoise(model, mu, sigma, n)

# Stats timer.
times0 = time.time()

gamma, m, de, slnprob, pos = stats.stats(gamma_true, m_true, de_true, zpicks, mag, sigma)

# Time taken by stats. 
times1=time.time()
timer.timer('stats', times0, times1)

# Time taken by the script. 
timet1=time.time()
timer.timer('script', timet0, timet1)
