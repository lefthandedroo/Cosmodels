#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 12:40:48 2018

@author: BallBlueMeercat


Cheat sheet:
    
    # Walker steps.
    m = sampler.flatchain[:,0]
    slnprob = sampler.flatlnprobability
    
    figure()
    xlabel('parameter value')
    ylabel('step number')
    plot(m, slnprob, '.', color='red')
    title('slnprob for m')
    show()
    
    # Chains.    
    figure()
    pl.title('flatChains with m_true in red')
    pl.plot(sampler.flatchain[:,0].T, '-', color='k', alpha=0.3)
    pl.axhline(m_true, color='red')
    pl.show
"""

import pickle
import os
from pylab import figure, plot, xlabel, ylabel, title, show

def walkersteps(filename):
    
    if os.path.exists(filename):
        with open(filename,'rb') as rfp: 
            results = pickle.load(rfp)
    
    ttdd = [item[0] for item in results]
    print(ttdd)
    
    pick = input('Enter ID of entry you want to plot: ')
    
    if not pick:
        m = results[-1][1]
        slnprob = results[-1][2]
    else:
        pick = int(pick)
        m = results[pick][1]
        slnprob = results[pick][2]
    
    figure()
    xlabel('parameter value')
    ylabel('step number')
    plot(m, slnprob, '.', color='red')
    title('slnprob for m')
    show()

    return m, slnprob

walkersteps('resultlist.p')


from pylab import grid, scatter

def modelcheck(t, mag, zpicks, dlpc, dl, gamma, ombar_m0, ombar_de0, a, ombar_m, ombar_de):
    
    # Plotting selected results:
    # a and a_dot vs time.

    while True:
        figure()
        xlabel('redshift $z$')
        ylabel('a')
        grid(True)
        plot(zpicks, a, 'r', lw=1)
        title('IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s'
              %(ombar_m0, ombar_de0, gamma))
        show()
        break

    # e_dashm and e_dashde vs time.
    while True:
        # e_dashm
        figure()
        xlabel('redshift $z$')
        ylabel('$\epsilon_m \'$')
        lw = 1
        plot(zpicks, ombar_m, 'g', linewidth=lw)
        title('$\epsilon_m \'$ evolution, IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s'
              %(ombar_m0, ombar_de0, gamma))
        show()

        # e_dashde
        figure()
        xlabel('redshift $z$')
        ylabel('$\epsilon_{DE} \'$')
        lw = 1
        plot(zpicks, ombar_m, 'm', linewidth=lw)
        title('$\epsilon_{DE} \'$ evolution, IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s'
              %(ombar_m0, ombar_de0, gamma))
        show()
        break


    # Luminosity distance vs redshift.
    while True:
        figure()
        xlabel('redshift $z$')
        ylabel('$d_L$*($H_0$/c)')
        grid(True)
        plot(zpicks, dl, 'tab:green', lw=1)
        title('$d_L$, IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s'
              %(ombar_m0, ombar_de0, gamma))
        show()
        break   

    while False:
        figure()
        xlabel('redshift $z$')
        ylabel('pc')
        grid(True)
        plot(zpicks, dlpc, 'tab:green', lw=1)
        title('$d_L$, IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s'
              %(ombar_m0, ombar_de0, gamma))
        show()
        break       

    dlgpc = dlpc /10**9
    while False:
        figure()
        xlabel('redshift $z$')
        ylabel('Gpc')
        grid(True)
        plot(zpicks, dlgpc, 'tab:green', lw=1)
        title('$d_L$, IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s'
              %(ombar_m0, ombar_de0, gamma))
        show()
        break


    while True:
        # Redshift vs time.
        figure()
        xlabel('redshift $z$')
        ylabel('$z$')
#        axis([0,-0.8,0,5])
        grid(True)
        plot(t, zpicks, 'tab:pink', lw=1)
        title('Redshift, IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s'
              %(ombar_m0, ombar_de0, gamma))
        show()
        break


    while True:
        figure()
        xlabel('redshift')
        ylabel('magnitude')
        title('Mag simulated with msim parameters')
        scatter(zpicks, mag, lw='3', c='r')
        show()
        break

    return