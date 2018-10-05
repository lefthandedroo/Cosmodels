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
        params = dictionary w/ key:value
            'gamma':int/float = interaction constant;
            'm':int/float = e_m(t)/ec(t0) at t=t0;
        zpicks = list of redshifts ;
        firstderivs_key = string, name of IVCDM model to use for model mag.
    
    """
#    print('@@ zodesolve has been called')


    # Inserting 0 at the front of redshifts to allow initial conditions.
    zpicks = [0.0] + zpicks
    
    # Standard cosmological parameters.
    H0 = 1.0
    c_over_H0 = 4167 * 10**6    # c/H0 in parsecs
    
    # Initial conditions at z = 0 (now).
    t0 = 0.0              # time
    a0 = 1.0            # scale factor
    z0 = 0.0              # redshift
    dl0 = 0.0             # luminosity distance
    rho_c0 = H0**2      # critical density
    ombar_m0 = params.get('m', 0)                        # e_m(z)/ec(z=0)
    gamma = params.get('gamma',0)
    zeta = params.get('zeta', 0) 
    ombar_de0 = params.get('de', rho_c0/rho_c0 -ombar_m0) # e_de(z)/ec(z=0)
    ombar_r0 = 0.0
    
    # ODE solver parameters:
    abserr = 1.0e-8
    relerr = 1.0e-6
    
        
    # Extracting the parsed mode of interaction.
    firstderivs_function = firstderivs_functions.get(firstderivs_key,0)
    if firstderivs_function == 0:
        print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
        print("firstderivs_functions dict didn't have the key zodeosolve asked for")
        print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
        
    if firstderivs_key == 'exotic':
        # Pack up the initial conditions and eq of state parameters.
        v0 = [t0, a0, ombar_m0, ombar_r0, ombar_de0, z0, dl0]
        
        # Call the ODE solver. 
        vsol = odeint(firstderivs_function, v0, zpicks, args=(gamma,zeta,H0), 
                      atol=abserr, rtol=relerr)        

        # Separate results into their own arrays:
        t = vsol[1:,0]
        a = vsol[1:,1]
        ombar_m = vsol[1:,2]
        ombar_r = vsol[1:,3]
        ombar_de = vsol[1:,4]
        z = vsol[1:,5]    
        dl = vsol[1:,6] * (1+z)  # in units of dl*(H0/c)
        dlpc = dl * c_over_H0    # dl in parsecs (= vsol[dl] * c/H0)
        
#        import matplotlib.pyplot as plt
#        plt.figure()
#        plt.xlabel('redshift $z$')
#        plt.ylabel(r'$\bar \Omega $')
#        plt.grid(True)
#        plt.plot(z, ombar_m, label=r'$\bar \Omega_{m}$', 
#                 color='xkcd:coral', linestyle=':')
#        plt.plot(z, ombar_de, label=r'$\bar \Omega_{DE}$', 
#                 color='xkcd:aquamarine')
#        plt.plot(z, ombar_r, label=r'$\bar \Omega_{r}$', 
#                 color='xkcd:black')
#        plt.legend()
#        plt.title(r'$\bar \Omega_{r}$ evolution, model = %s, $\gamma$ = %s'
#              %(firstderivs_key, gamma))        
        
        plot_var = t, dlpc, dl, a, ombar_m, ombar_r, gamma, zeta, ombar_de, ombar_m0, ombar_r0, ombar_de0
        
    else:
        # Pack up the initial conditions and eq of state parameters.
        v0 = [t0, a0, ombar_m0, ombar_de0, z0, dl0] 
    
        # Call the ODE solver. 
        vsol = odeint(firstderivs_function, v0, zpicks, args=(gamma,H0), 
                      atol=abserr, rtol=relerr)
            
        # Separate results into their own arrays:
        t = vsol[1:,0]
        a = vsol[1:,1]
        ombar_m = vsol[1:,2]
        ombar_de = vsol[1:,3]
        z = vsol[1:,4]    
        dl = vsol[1:,5] * (1+z)  # in units of dl*(H0/c)
        dlpc = dl * c_over_H0    # dl in parsecs (= vsol[dl] * c/H0)
        
        plot_var = t, dlpc, dl, a, ombar_m, gamma, ombar_de, ombar_m0, ombar_de0
    
    return dlpc, plot_var

#def odesolve(params, zpicks, firstderivs_key):
#    """
#    Takes in:
#        gamma = interaction constant;
#        m = e_m(t)/ec(t0) at t=t0;
##        de = e_de(t)/ec(t0) at t=t0.
#    Returns: 
#        z = numpoints number of redshifts zmin<z<zmax;
#        dlpc = luminosity distance in pc.
#    
#    """
##    print('@@ zodesolve has been called')
#
#    # Inserting 0 at the front of redshifts to allow initial conditions.
#    zpicks = [0.0] + zpicks
#    
#    # Standard cosmological parameters.
#    H0 = 1
#    c_over_H0 = 4167 * 10**6    # c/H0 in parsecs
#    
#    # Initial conditions at z = 0 (now).
#    dl0 = 0             # luminosity distance
#    rho_c0 = H0**2      # critical density
#    ombar_m0 = params.get('m', 0)                        # e_m(z)/ec(z=0)
#    gamma = params.get('gamma',0)
#    ombar_de0 = params.get('de', rho_c0/rho_c0 -ombar_m0)# e_de(z)/ec(z=0)
#    
#    # ODE solver parameters:
#    abserr = 1.0e-8
#    relerr = 1.0e-6
#    
#    # Pack up the initial conditions and eq of state parameters.
#    v0 = [ombar_m0, ombar_de0, dl0]
#    
#    # Extracting the parsed mode of interaction.
#    firstderivs_function = firstderivs_functions.get(firstderivs_key,0)
#    if firstderivs_function == 0:
#        print("firstderivs_functions dict didn't have the key zodeosolve asked for")
#    
#    # Call the ODE solver.
#    vsol = odeint(firstderivs_function, v0, zpicks, args=(gamma,H0), 
#                  atol=abserr, rtol=relerr)
#            
#    # Separate results into their own arrays:
#    z = np.asarray(zpicks)
#    z = z[1:]
#    ombar_m = vsol[1:,0]
#    ombar_de = vsol[1:,1]
#    dl = vsol[1:,2] * (1+z)  # in units of dl*(H0/c)
#    dlpc = dl * c_over_H0    # dl in parsecs (= vsol[dl] * c/H0)
#    
#    plot_var = dlpc, dl, ombar_m, gamma, ombar_de, ombar_m0, ombar_de0
#    
#    return dlpc, plot_var