#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main code
"""
import numpy as np
import time
import os

import tools
import datasim
import stats

def paramfinder(npoints, nsteps, sigma, mu, m_true, directory):
    # Script timer.
    timet0 = time.time()
<<<<<<< HEAD
    
    directory = str(int(time.time()))
    
    if not os.path.exists(directory):
        os.makedirs(directory)
=======
>>>>>>> 8f82b3f0fa02bbdb918431688fd1b124eb201c6d
        
    # Generating redshifts to simulate mag.
    import zpicks # DO MOVE ABOVE FUNCTION
    zpicks = zpicks.zpicks(0.005, 2, npoints)
    
    # Generating apparent magnitues mag at redshifts z=0 to z < zmax.
    model = datasim.mag(m_true, zpicks)
    model = np.asarray(model)
    mag = datasim.gnoise(model, mu, sigma)
    
    # emcee parameter search.
<<<<<<< HEAD
    propert, sampler = stats.stats(m_true, zpicks, mag, sigma, nsteps)
    
    
#    # Saving sampler to file.
#    import outputsave
#    outputsave.samplersave(sampler)
=======
    propert, sampler = stats.stats(m_true, zpicks, mag, sigma, nsteps, directory)
>>>>>>> 8f82b3f0fa02bbdb918431688fd1b124eb201c6d
    
    # Time taken by script. 
    timet1=time.time()
    tools.timer('script', timet0, timet1)
    
    return propert, sampler