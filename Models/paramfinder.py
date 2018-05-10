#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main code
"""
import numpy as np
import time
import os

import timer
import gnoise
import zmsim
import stats

def paramfinder(npoints, nsteps, sigma, mu, m_true):
    # Script timer.
    timet0 = time.time()
    
    directory = str(int(time.time()))
    
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    # Generating redshifts to simulate mag.
    import zpicks # DO MOVE ABOVE FUNCTION
    zpicks = zpicks.zpicks(0.005, 2, npoints)
    
    # Generating apparent magnitues mag at redshifts z=0 to z < zmax.
    model = zmsim.zmsim(m_true, zpicks)
    model = np.asarray(model)
    mag = gnoise.gnoise(model, mu, sigma)
    
    # emcee parameter search.
    propert, sampler = stats.stats(m_true, zpicks, mag, sigma, nsteps)
    
    
#    # Saving sampler to file.
#    import outputsave
#    outputsave.samplersave(sampler)
    
    # Time taken by script. 
    timet1=time.time()
    timer.timer('script', timet0, timet1)
    
    return propert, sampler