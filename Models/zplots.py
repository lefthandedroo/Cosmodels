#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 12:40:48 2018

@author: BallBlueMeercat
"""
from pylab import figure, plot, xlabel, ylabel, grid, title, show, scatter


def zplots(t, mag, zpicks, dlpc, dl, gamma, e_dash0m, e_dash0de, a, e_dashm, e_dashde):

#    f = open('results.pckl', 'rb')
#    obj = pickle.load(f)
#    f.close()
#    
    
    
    # Plotting selected results:
    # a and a_dot vs time.

    while True:
        figure()
        xlabel('redshift $z$')
        ylabel('a')
        grid(True)
        plot(zpicks, a, 'r', lw=1)
        title('IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s'
              %(e_dash0m, e_dash0de, gamma))
        show()
        break

    # e_dashm and e_dashde vs time.
    while True:
        # e_dashm
        figure()
        xlabel('redshift $z$')
        ylabel('$\epsilon_m \'$')
        lw = 1
        plot(zpicks, e_dashm, 'g', linewidth=lw)
        title('$\epsilon_m \'$ evolution, IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s'
              %(e_dash0m, e_dash0de, gamma))
        show()

        # e_dashde
        figure()
        xlabel('redshift $z$')
        ylabel('$\epsilon_{DE} \'$')
        lw = 1
        plot(zpicks, e_dashde, 'm', linewidth=lw)
        title('$\epsilon_{DE} \'$ evolution, IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s'
              %(e_dash0m, e_dash0de, gamma))
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
              %(e_dash0m, e_dash0de, gamma))
        show()
        break   

    while False:
        figure()
        xlabel('redshift $z$')
        ylabel('pc')
        grid(True)
        plot(zpicks, dlpc, 'tab:green', lw=1)
        title('$d_L$, IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s'
              %(e_dash0m, e_dash0de, gamma))
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
              %(e_dash0m, e_dash0de, gamma))
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
              %(e_dash0m, e_dash0de, gamma))
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