#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 12:40:48 2018

@author: BallBlueMeercat
"""
from matplotlib.font_manager import FontProperties
from pylab import figure, plot, xlabel, ylabel, grid, legend, title, show, scatter

def plots(mag, zpicks, z, dlpc, dl, gamma, e_dash0m, e_dash0de, t, a, a_dot, t_cut, a_cut, a_dotcut, e_dashm, e_dashde):

    # Plotting selected results:
    # a and a_dot vs time.
    while True:
        figure()
        xlabel('time in $H_0^{-1}$')
        grid(True)
        plot(t_cut, a_cut, 'r', t_cut, a_dotcut, 'b', lw=1)
        legend((r'$a$', r'$\.a$'), prop=FontProperties(size=16))
        title('IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s [cut]'
              %(e_dash0m, e_dash0de, gamma))
        show()
        break


    # e_dashm and e_dashde vs time.
    while True:
        # e_dashm
        figure()
        xlabel('t in 1/H0')
        ylabel('$\epsilon_m \'$')
        lw = 1
        plot(t_cut, e_dashm, 'g', linewidth=lw)
        title('$\epsilon_m \'$ evolution, IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s [cut]'
              %(e_dash0m, e_dash0de, gamma))
        show()

        # e_dashde
        figure()
        xlabel('t in 1/H0')
        ylabel('$\epsilon_{DE} \'$')
        lw = 1
        plot(t_cut, e_dashde, 'm', linewidth=lw)
        title('$\epsilon_{DE} \'$ evolution, IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s [cut]'
              %(e_dash0m, e_dash0de, gamma))
        show()
        break


    # Luminosity distance vs redshift.
    while True:
        figure()
        xlabel('redshift $z$')
        ylabel('$d_L$*($H_0$/c)')
        grid(True)
        plot(z, dl, 'tab:green', lw=1)
        title('$d_L$, IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s [cut]'
              %(e_dash0m, e_dash0de, gamma))
        show()
        break   

    while False:
        figure()
        xlabel('redshift $z$')
        ylabel('pc')
        grid(True)
        plot(z, dlpc, 'tab:green', lw=1)
        title('$d_L$, IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s [cut]'
              %(e_dash0m, e_dash0de, gamma))
        show()
        break       

    dlgpc = dlpc /10**9
    while False:
        figure()
        xlabel('redshift $z$')
        ylabel('Gpc')
        grid(True)
        plot(z, dlgpc, 'tab:green', lw=1)
        title('$d_L$, IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s [cut]'
              %(e_dash0m, e_dash0de, gamma))
        show()
        break


    while True:
        # Redshift vs time.
        figure()
        xlabel('time in $H_0^{-1}$')
        ylabel('$z$')
#        axis([0,-0.8,0,5])
        grid(True)
        plot(t_cut, z, 'tab:pink', lw=1)
        title('Redshift, IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s [cut]'
              %(e_dash0m, e_dash0de, gamma))
        show()
        break


    # Complete results with blow up resulting from a approaching big bang.
    while True:  
        figure()
        xlabel('time in $H_0^{-1}$')
        grid(True)
        plot(t, a, 'r', lw=1)
        plot(t, a_dot, 'b', lw=1)
        legend((r'$a$', r'$\.a$', r'$\'\epsilon$'), prop=FontProperties(size=16))
        title('IC: $\epsilon_{m0} \'$ = %s, $\epsilon_{DE0} \'$ =%s, $\gamma$ = %s [uncut]'
              %(e_dash0m, e_dash0de, gamma))
        show()
        break

    while True:
        figure()
        title('Mag simulated with msim parameters')
        scatter(zpicks, mag, lw='3', c='r')
        show()
        break

    return