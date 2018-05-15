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

from pylab import figure, plot, xlabel, ylabel, title, show, grid, scatter
import numpy as np

from results import load

def onepercent(directory):
    directory = str(directory)
    
    results = []

    vc = load(directory, 'vc.p')
    sd = load(directory, 'sd.p')
    mean  = load(directory, 'mean.p')
    sigma = load(directory, 'sigma.p')
    npoints = load(directory, 'npoints.p')
    
#    vc = [0.01, 0, 0.02, 3, 5, 4, 6]
#    sd = [20, 10, 30, 40, 60, 50, 70]
#    mean = [0, 1, 2, 3, 4, 5, 6]
#    sigma = [0, 1, 2, 3, 4, 5, 6] 
#    npoints = [0, 1, 2, 3, 4, 5, 6]

    results.append(vc)
    results.append(sd)
    results.append(mean)
    results.append(sigma)
    results.append(npoints)
    
#    print('vc',len(vc))
#    print('sd',len(sd))
#    print('mean',len(mean))
#    print('sigma',len(sigma))
#    print('npoints',len(npoints))
    
    results = np.column_stack(results)
    results = results[results[:,0].argsort()]
    
    vc, sd, mean, sigma, npoints = np.hsplit(results,5)
#    vc = vc.flatten()
    
    print('vc')
    print(vc)
    
    si = np.where(vc < 0.08)
    print('si')
    print(si)
    
    print('vc[si]')
    print(vc[si])
    
    figure()
    xlabel('datazet size')
    ylabel('sigma of noise added to data')
    title('runs where m was found within 1%')
    scatter(npoints[si], sigma[si])
    show()
    
    return results, vc, sd, mean, sigma, npoints

results, vc, sd, mean, sigma, npoints = onepercent(1526341352)


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