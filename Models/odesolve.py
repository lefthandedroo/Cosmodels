#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:38:48 2018

@author: BallBlueMeercat
"""

import firstderivs
import numpy as np
from scipy.integrate import odeint

import matplotlib.pyplot as pl
from matplotlib.font_manager import FontProperties
from pylab import figure, plot, xlabel, grid, legend, title


# Standard cosmological parameters.
H0 = 1       # Hubble parameter at t=now
tH = 1.0/H0  # Hubble time
# Eq of state parameters for known fluids:
w_r = 1/3     # radiation
w_m = 0.0     # matter
w_de = -1.0   # cosmological constant (dark energy?)

c_over_H0 = 4167 * 10**6    # c/H0 in parsecs


def odesolve(lamb,m,de):
    """
    Takes in:
        lamb = e_lamb(t)/ec(t0) at t=t0;
        m = e_m(t)/ec(t0) at t=t0;
        de = e_de(t)/ec(t0) at t=t0.
    Returns: 
        dlmpc = luminosity distance in Mpc;
        z = redshift under 2.
    
    """
#    print('@@ odesolve has been called')
    # Last value for a before results are considered close enough to z = 2.
    a_d = 0.25
    
    # Time (in 1/H0) to integrate until.  If this time isn't long enough for a to 
    # decrease to a_d then stoptime will be extended by time until a_d is reached.
    # 0.665 matter only, 0.96 standard m+de
    time = 0.9
    
    
    # Initial conditions at time = t0.
    a0 = 1.0        # scale factor
    a_dot0 = 1.0    # speed of expansion
    e_dash0m = m    # e_m(t)/ec(t0)
    e_dash0de = de  # e_de(t)/ec(t0)
    z0 = 0
    dl0 = 0
    
    # ODE solver parameters:
    abserr = 1.0e-8
    relerr = 1.0e-6
    numpoints = 1000000
    
    stoptime = -time # Integrating back in time as time now is t0.
    
    while True:
        # Create time samples for the ODE solver.
        t = [stoptime * tH * float(i) / (numpoints - 1) for i in range(numpoints)]
#        print('time is : ',t[0])
        # Pack up the initial conditions and eq of state parameters.
        v0 = [a0, a_dot0, e_dash0m, e_dash0de, z0, dl0]
        w = [w_m, w_de]
        
        # Call the ODE solver. maxstep=5000000 added later to try and avoid 
        # ODEintWarning: Excess work done on this call (perhaps wrong Dfun type).
        vsol = odeint(firstderivs.firstderivs, v0, t, args=(w,lamb,), atol=abserr, rtol=relerr, mxstep=5000000)
        # vsol type is:  <class 'numpy.ndarray'>
                
        # Remove unwanted results which are too close to big bang from the plot.
        # Separate results into their own arrays:
        a = vsol[:,0]
        a_dot = vsol[:,1]
#        e_dashm = vsol[:,2]
#        e_dashde = vsol[:,3]
        z = vsol[:,4]
        dl = vsol[:,5] * (1+z)
        dlmpc = dl * c_over_H0    # dl in Mega parsecs (= vsol[dl] * c/H0)
    
        
        # Find where results start to get strange (smaller than a_d):
        blowups = np.where(a < a_d)    # Tuple with indecies of a so
                                       # small that other results blow up.                             
        blowups = np.asarray(blowups)  # Converting to np array.
    
        if blowups.any():              # Check if instances of a < a_d exist.   
            blowup = blowups[0,0]
        else:                          # If no instance of a < a_d was found
            stoptime -= time           # then integrate further back in time.
            continue
        
        
        # Remove the values after the index of first instance of a < a_d.
        t_cut = np.asarray(t)
        
        t_cut = t_cut[:blowup]
        a = a[:blowup]
        a_dot = a_dot[:blowup]
        z = z[:blowup]
        dl = dl[:blowup]
        dlmpc = dlmpc[:blowup]
                
        
        # Age of the universe.
        age = t_cut[np.argmin(t_cut)]
        age = -round(age, 2)

        # Plotting selected results:
        # a and a_dot vs time.
        while False:
            figure()
            xlabel('time in $H_0^{-1}$')
            grid(True)
            plot(t_cut, a, 'r', t_cut, a_dot, 'b', lw=1)
            legend((r'$a$', r'$\.a$'), prop=FontProperties(size=16))
            title('Cut results for $\omega$ = %s, $\lambda$ = %s, age = %s $H_0^{-1}$'
                  %(w,lamb,age))
            break
        
        # Luminosity distance dl vs redshift.
        while False:
            figure()
            xlabel('redshift $z$')
            grid(True)
            plot(z, dl, 'tab:green', lw=1)
            title('$D_L$ vs redshift for $\omega$ = %s, $\lambda$ = %s,'
                  ' age = %s $H_0^{-1}$'%(w,lamb,age))
            break
        
        while False:
            # Redshift vs time.
            figure()
            xlabel('time in $H_0^{-1}$')
            pl.axis([0,-0.1,0,5])
            grid(True)
            plot(t_cut, z, 'tab:pink', lw=1)
            title('Redshift evolution for $\omega$ = %s, $\lambda$ = %s,'
              ' age = %s $H_0^{-1}$'%(w,lamb,age))
            break
            
        break
            
    
#    # Complete results with blow up resulting from a approaching big bang.
    while False:  
        figure()
        xlabel('time in $H_0^{-1}$')
        grid(True)
        
        # Plotting complete results.
        plot(t, vsol[:,0], 'r', lw=1)
        plot(t, vsol[:,1], 'b', lw=1)
        
        legend((r'$a$', r'$\.a$', r'$\'\epsilon$'), prop=FontProperties(size=16))
        title('Complete results for $\omega$ = %s'%(w))
        break

    return z, dlmpc