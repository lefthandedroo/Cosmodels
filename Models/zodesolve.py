7#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:38:48 2018

@author: BallBlueMeercat
"""
from scipy.integrate import odeint
import firstderivs_cython as f
import numpy as np

firstderivs_functions = {
        'dlLCDM':f.dlLCDM,
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
        zpicks = list of redshifts ;
        model = string, name of model being tested.

    """

    # Inserting 0 at the front of redshifts to use initial conditions.
    zpicks = [0.0] + zpicks

    # Standard cosmological parameters.
    H0 = 1.0
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
        for i in range(index,len(values)):
            plot_var[names[i]] = values[i]
        # Separate results into their own arrays:
        plot_var['t'] = vsol[1:,0]
        plot_var['a'] = vsol[1:,1]
        plot_var['ombar_m'] = vsol[1:,2]
        if model == 'exotic' or model == 'waterfall':
            plot_var['ombar_r'] = vsol[1:,3]
            plot_var['ombar_r0'] = values[2]
            if model == 'waterfall':
                plot_var['a_ombar'] = vsol[1:,4]
                plot_var['b_ombar'] = vsol[1:,5]
                plot_var['c_ombar'] = vsol[1:,6]
        plot_var['ombar_de0'] = ombar_de0
        plot_var['ombar_de'] = vsol[1:,-3]
        plot_var['z'] = vsol[1:,-2]
        plot_var['dl'] = vsol[1:,-1] * (1+z)


    return dlpc, plot_var