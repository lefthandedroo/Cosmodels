#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main code
"""
import numpy as np
import time

import tools
import datasim
import stats

def paramfinder(npoints, nsteps, sigma, mu, m_true, directory):
    # Script timer.
    timet0 = time.time()
        
    # Generating redshifts to simulate mag.
    import zpicks # DO MOVE ABOVE FUNCTION
    zpicks = zpicks.zpicks(0.005, 2, npoints)
    
    # Generating apparent magnitues mag at redshifts z=0 to z < zmax.
    model = datasim.mag(m_true, zpicks)
    model = np.asarray(model)
    mag = datasim.gnoise(model, mu, sigma)
    
    # emcee parameter search.
    propert, sampler = stats.stats(m_true, zpicks, mag, sigma, nsteps, directory)
    
    # Time taken by script. 
    timet1=time.time()
    tools.timer('script', timet0, timet1)
    
    return propert, sampler