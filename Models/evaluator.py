#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 21:45:40 2018

@author: BallBlueMeercat
"""
from pylab import figure, plot, xlabel, ylabel, title, show
import matplotlib.pyplot as pl


import paramfinder

# Model parameteres:  
gamma_true = 0  # Interaction term, rate at which DE decays into matter.
m_true = 0.3    # (= e_m(t)/e_crit(t0) at t=t0).
de_true = 0.7   # (de = e_de(t)/e_crit(t0) at t=t0).

def evaluator():
#    sigma = 0.01
#    nsteps = 100
    
    gamma = []
    m = []
    de = []
    numpoints = []
    
    n = 1
    while n < 10000000:
        n = n * 10
        gammabest, mbest, debest = paramfinder.paramfinder(n)
        gamma.append(gammabest)
        m.append(mbest)
        de.append(debest)
        numpoints.append(n)
    
    figure()
    xlabel('parameter value')
    ylabel('number of datapoints used')
    pl.axhline(gamma_true, color='red')
    plot(gamma, numpoints, '.')
    show()
    
    figure()
    xlabel('parameter value')
    ylabel('number of datapoints used')
    pl.axhline(m_true, color='red')
    plot(m, numpoints, '.')
    show()
    
    figure()
    xlabel('parameter value')
    ylabel('number of datapoints used')
    pl.axhline(de_true, color='red')
    plot(de, numpoints, '.')
    show()