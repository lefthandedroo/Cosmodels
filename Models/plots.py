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
import matplotlib.pyplot as plt
import numpy as np
import os

from results import load

def onepercent():
    
    direclist = []
    for d in os.walk('./results/'):
        direclist.append(d[0])
    direclist.pop(0)
    
    vc = []
    sigma =[]
    npoints = []
    for directory in direclist:
            vc += load(directory, 'vc.p')
            sigma += load(directory, 'sigma.p')
            npoints += load(directory, 'npoints.p')

    vc = np.asarray(vc)
    sigma = np.asarray(sigma)
    npoints = np.asarray(npoints)
      
    pi = np.where(vc < 1)   # Indicies of rows with vc < 1%.
    pinpoints = npoints[pi]
    pisigma = sigma[pi]
    
    sinpoints = []  # Single results, removing doubles of nsteps.
    sisigma = []
    
    for i in range(len(pinpoints)):
    
        if pinpoints[i] in sinpoints:
           index = np.where(sinpoints == pinpoints[i])
           index = int(index[0])
           if pisigma[i] > sisigma[index]:
              sisigma[index] = pisigma[i]

        else:
            sinpoints.append(pinpoints[i])
            sisigma.append(pisigma[i])

    
#    ind = np.ones((10,), bool)
#    ind[n] = False
#    A1 = A[ind,:]
    figure()
    xlabel('datazet size')
    ylabel('sigma of noise added to data')
    title('noisiest runs where m was found within 1%')       
    plt.scatter(sinpoints, sisigma, c='m', label='1% sd on m')
    plt.scatter(npoints, sigma, c='c', marker='x', label='all runs')
    plt.legend()
    plt.show()
    
#    figure()
#    xlabel('datazet size')
#    ylabel('sigma of noise added to data')
#    title('noisiest runs where m was found within 1%')
#    scatter(sinpoints, sisigma, c='m')
#    show()
#    
#    figure()
#    xlabel('datazet size')
#    ylabel('sigma of noise added to data')
#    title('all runs')
#    scatter(npoints, sigma, c='c')
#    show()
    
    return vc, sigma, npoints

vc, sigma, npoints = onepercent()


def modelcheck(t, mag, zpicks, dlpc, dl, 
               gamma, ombar_m0, ombar_de0, a, ombar_m, ombar_de):
    
    # Plotting selected results:
    # a and a_dot vs time.

    while True:
        figure()
        xlabel('redshift $z$')
        ylabel('a')
        grid(True)
        plot(zpicks, a, 'r', lw=1)
        title('IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0}',
              ' \'$ =%s, $\gamma$ = %s'%(ombar_m0, ombar_de0, gamma))
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
        title('$\epsilon_m \'$ evolution, IC: $\epsilon_{m0}',
              ' ,\'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s'
              %(ombar_m0, ombar_de0, gamma))
        show()

        # e_dashde
        figure()
        xlabel('redshift $z$')
        ylabel('$\epsilon_{DE} \'$')
        lw = 1
        plot(zpicks, ombar_m, 'm', linewidth=lw)
        title('$\epsilon_{DE} \'$ evolution, IC: $\epsilon_{m0}',
              ' \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s'
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
        title('$d_L$, IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0}',
              ' \'$ =%s, $\gamma$ = %s'
              %(ombar_m0, ombar_de0, gamma))
        show()
        break   

    while False:
        figure()
        xlabel('redshift $z$')
        ylabel('pc')
        grid(True)
        plot(zpicks, dlpc, 'tab:green', lw=1)
        title('$d_L$, IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0}',
              ' \'$ =%s, $\gamma$ = %s'
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
        title('$d_L$, IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0}',
              ' \'$ =%s, $\gamma$ = %s'
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
        title('Redshift, IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0}',
              ' \'$ =%s, $\gamma$ = %s'
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