7#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:38:48 2018

@author: BallBlueMeercat
"""
from scipy.integrate import odeint
import firstderivs_cython as f
import numpy as np
import matplotlib.pyplot as plt

firstderivs_functions = {
        'stepfall':f.stepfall,
        'waterfall':f.waterfall,
        'rainbow':f.rainbow,
        'exotic':f.exotic,
        'late_intxde':f.late_intxde,
        'heaviside_late_int':f.heaviside_late_int,
        'late_int':f.late_int,
        'expgamma':f.expgamma,
        'txgamma':f.txgamma,
        'zxgamma':f.zxgamma,
        'gamma_over_z':f.gamma_over_z,
        'zxxgamma':f.zxxgamma,
        'gammaxxz':f.gammaxxz,
        'rdecay_m':f.rdecay_m,
        'rdecay_de':f.rdecay_de,
        'rdecay_mxde':f.rdecay_mxde,
        'rdecay':f.rdecay,
        'interacting':f.interacting,
        'LCDM':f.LCDM
        }

def zodesolve(names, values, zpicks, model, plot_key):
    """
    Takes in:
        names = list of strings, names of parameters to be fitted;
        values = np.array, values of parameters;
        zpicks = np.ndarray of redshifts ;
        model = string, name of model being tested.

    """

    # Inserting 0 at the front of redshifts to use initial conditions.
    zpicks = np.insert(zpicks, 0, 0.0)

    # Standard cosmological parameters.
    H0 = 1.0
    c = 1.0
    c_over_H0 = 4167 * 10**6    # c/H0 in parsecs

    # Initial conditions at z = 0 (now).
    t0 = 0.0              # time
    a0 = 1.0              # scale factor
    z0 = 0.0              # redshift
    dl0 = 0.0             # luminosity distance
    rho_c0 = H0**2        # critical density

    # Pack up the initial conditions and interaction terms.
    int_terms = []

    if model == 'waterfall':
        index = 6
    elif model == 'stepfall':
        index = 4
    elif model == 'exotic':
        index = 3
    elif model == 'LCDM':
        index = len(values)
    else:
        index = 2

    int_terms = values[index:]
    fluids = values[1:index]
    ombar_de0 = rho_c0/rho_c0 - np.sum(fluids)

    t0a0 = np.array([t0, a0])
    de0z0dl0 = np.array([ombar_de0, z0, dl0])
    # Remember that you lost precision when concatenating arr over using a list.
    v0 = np.concatenate((t0a0, fluids, de0z0dl0))

    # Extracting the parsed mode of interaction.
    firstderivs_function = firstderivs_functions.get(model,0)
    if firstderivs_function == 0:
        raise ValueError("zodesolve doesn't have this firstderivs_key at the top")

    # Call the ODE solver.
    vsol = odeint(firstderivs_function, v0, zpicks, args=(int_terms,H0),
                  atol=1.0e-8, rtol=1.0e-6)
    z = vsol[1:,-2]
    dl = vsol[1:,-1] * (1+z)  # in units of dl*(H0/c)
    dlpc = dl * c_over_H0    # dl in parsecs (= vsol[dl] * c/H0)

    plot_var = {}
    if plot_key:
        # Separate results into their own arrays:
        plot_var['t'] = vsol[1:,0]
        plot_var['a'] = vsol[1:,1]
        
        # Collecting fluids and their names for plotting:
        fluid_arr = np.zeros(((index), (len(zpicks)-1)))
        fluid_names = []
        for i in range((index-1)):
            fluid_names.append(names[i+1])
            fluid_arr[i] = vsol[1:,(i+2)]
        fluid_names.append('de_ombar')
        fluid_arr[-1] = vsol[1:,-3]
        plot_var['fluid_names'] = fluid_names
        plot_var['fluid_arr'] = fluid_arr
        
        plot_var['z'] = vsol[1:,-2]
        plot_var['dl'] = vsol[1:,-1] * (1+z)
        plot_var['int_terms'] = int_terms
        
        da = dl * (1.0+z)**(-2.0)
        plot_var['da'] = da        
#        plt.figure()
#        plt.title('Angular diameter distance evolution')
#        plt.xlabel('z')
#        plt.ylabel(r'$ \left( \frac{H_0}{c} \right) d_A $', fontsize=15, labelpad=10)
#        plt.plot(z, da)
        
        Hz = H0 * (np.sum(fluid_arr, axis=0))**(0.5)
        plot_var['Hz'] = Hz
        
        da = dlpc/10**6 * (1.0+z)**(-2.0)
        dv = ((1+z)**2 * da**2 * c*z/Hz)**(1/3)        
        plot_var['dv'] = dv
#        plt.figure()
#        plt.title(r'$D_v$ evolution')
#        plt.xlabel('z')
#        plt.ylabel(r'$ D_v (z)$ [Mpc]')
#        plt.plot(z, dv)
        
        plt.show()

    return dlpc, plot_var