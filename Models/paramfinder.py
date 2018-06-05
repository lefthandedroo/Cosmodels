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

def paramfinder(npoints, nsteps, sigma, mu, params, save_path, firstderivs_key):
    # Script timer.
    timet0 = time.time()
        
    # Generating redshifts to simulate mag.
    import zpicks # DO NOT MOVE THIS LINE UP
    zpicks = zpicks.zpicks(0.005, 2, npoints)
    
    # Generating apparent magnitues mag at redshifts z=0 to z < zmax.
    model = datasim.mag(params, zpicks, 'LCDM')
    model = np.asarray(model)
    mag = datasim.gnoise(model, mu, sigma)
    
    # emcee parameter search.
    propert, sampler = stats.stats(
            params, zpicks, mag, sigma, nsteps, save_path, firstderivs_key)
    
    # Time taken by script. 
    timet1=time.time()
    tools.timer('script', timet0, timet1)
    
    return propert, sampler