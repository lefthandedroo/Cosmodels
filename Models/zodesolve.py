7#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:38:48 2018

@author: BallBlueMeercat
"""
from scipy.integrate import odeint
import firstderivs_cython as f
#import firstderivs as f


firstderivs_functions = {
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

def zodesolve(params, zpicks, firstderivs_key):
    """
    Takes in:
        params = list of dictionaries {string:value} of names and 
        starting values of parameters to be emcee fitted:
            [{'m':int/float} = e_m(t)/ec(t0) at t=t0;
            {'Mcorr':int/float} = Corrected magnitude for SN;
            {'gamma':int/float}] = interaction constant;
        zpicks = list of redshifts ;
        firstderivs_key = string, name of IVCDM model to use for model mag.
    
    """
    plot_var = {}
    
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
    ombar_m0 = params[0].get('matter', 0.0)    # e_m(z)/ec(z=0)
    plot_var['ombar_m0'] = ombar_m0
    if firstderivs_key != 'LCDM':
        ombar_r0 = params[2].get('radiation', 0.0)
        
        if firstderivs_key == 'rainbow':
            ombar_a = params[3].get('ombar_a', 0)
            ombar_b = params[4].get('ombar_b', 0)
            ombar_c = params[5].get('ombar_c', 0)
            ombar_d = params[6].get('ombar_d', 0)
            ombar_e = params[7].get('ombar_e', 0)
            ombar_f = params[8].get('ombar_f', 0)
            ombar_g = params[9].get('ombar_g', 0)
            ombar_h = params[10].get('ombar_h', 0)
            ombar_i = params[11].get('ombar_i', 0)
            ombar_de0 = rho_c0/rho_c0 -ombar_m0 -ombar_r0 -ombar_a -ombar_b 
            ombar_de0 = ombar_de0 -ombar_c -ombar_d -ombar_e -ombar_f -ombar_g
            ombar_de0 = ombar_de0 -ombar_h - ombar_i
        elif firstderivs_key == 'waterfall':
            ombar_a = params[3].get('ombar_a', 0)
            ombar_b = params[4].get('ombar_b', 0)
            ombar_c = params[5].get('ombar_c', 0)
            ombar_de0 = rho_c0/rho_c0 -ombar_m0 -ombar_r0 -ombar_a -ombar_b -ombar_c
    else:
        ombar_de0 = rho_c0/rho_c0 - ombar_m0

    # Packing up interaction terms:
    int_terms = []
    if firstderivs_key == 'rainbow':
         for i in range(12,len(params)):
            for key in params[i]:
                int_terms.append(params[i].get(key,0))
                plot_var[key] = params[i][key]     
    elif firstderivs_key == 'waterfall':
         for i in range(6,len(params)):
            for key in params[i]:
                int_terms.append(params[i].get(key,0))
                plot_var[key] = params[i][key]         
    else:
        for i in range(2,len(params)):
            for key in params[i]:
                int_terms.append(params[i].get(key,0))
                plot_var[key] = params[i][key]
        
    # Extracting the parsed mode of interaction.
    firstderivs_function = firstderivs_functions.get(firstderivs_key,0)
    if firstderivs_function == 0:
        raise ValueError("zodesolve doesn't have this firstderivs_key at the top")
    
    # Pack up the initial conditions and eq of state parameters.
    if firstderivs_key == 'exotic':
        v0 = [t0, a0, ombar_m0, ombar_r0, ombar_de0, z0, dl0]
    elif firstderivs_key == 'rainbow':
        v0 = [t0, a0, ombar_m0, ombar_r0, ombar_de0, z0, dl0, 
              ombar_a, ombar_b, ombar_c, ombar_d, ombar_e, ombar_f, ombar_g, 
              ombar_h, ombar_i]
    elif firstderivs_key == 'waterfall':
        v0 = [t0, a0, ombar_m0, ombar_r0, ombar_a, ombar_b, ombar_c, 
              ombar_de0, z0, dl0]        
    else:
        v0 = [t0, a0, ombar_m0, ombar_de0, z0, dl0]
    
    # Call the ODE solver. 
    vsol = odeint(firstderivs_function, v0, zpicks, args=(int_terms,H0), 
                  atol=1.0e-8, rtol=1.0e-6)      

    # Separate results into their own arrays:
    plot_var['t'] = vsol[1:,0]
    plot_var['a'] = vsol[1:,1]
    plot_var['ombar_m'] = vsol[1:,2]        
    if firstderivs_key == 'exotic' or firstderivs_key == 'rainbow' or firstderivs_key == 'waterfall':  
        plot_var['ombar_r'] = vsol[1:,3]
        plot_var['ombar_r0'] = ombar_r0

    plot_var['ombar_de0'] = ombar_de0
    plot_var['ombar_de'] = vsol[1:,-3]          
    plot_var['z'] = vsol[1:,-2]    
    plot_var['dl'] = vsol[1:,-1] * (1+plot_var['z'])  # in units of dl*(H0/c)
    dlpc = plot_var['dl'] * c_over_H0    # dl in parsecs (= vsol[dl] * c/H0)
   
    
    return dlpc, plot_var