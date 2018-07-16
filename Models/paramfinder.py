#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main code
"""

import time
import tools
import stats

def paramfinder(npoints, nsteps, sigma, mu, params, 
                zpicks, mag, save_path, firstderivs_key):
    # Script timer.
    timet0 = time.time()
    
    # emcee parameter search.
    propert, sampler = stats.stats(params, zpicks, mag, sigma, nsteps, 
                                   save_path, firstderivs_key)
    
    # Time taken by script. 
    timet1=time.time()
    tools.timer('script', timet0, timet1)
    
    return propert, sampler