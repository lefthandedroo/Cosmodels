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
        'heaviside_sudden':f.heaviside_sudden,
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
        'LCDM':f.LCDM,
        'rLCDM':f.rLCDM
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

    if model == 'rainbow':
        int_in = 12
    elif model == 'waterfall':
        int_in = 6
    elif model == 'stepfall':
        int_in = 4
    elif model == 'exotic':
        int_in = 3
    elif model == 'LCDM' or model == 'rLCDM':
        int_in = len(values)
    else:
        int_in = 2

    int_terms = values[int_in:]
    fluids = values[1:int_in]
    ombar_de0 = rho_c0/rho_c0 - np.sum(fluids)

    t0a0 = np.array([t0, a0])
    de0z0dl0 = np.array([ombar_de0, z0, dl0])

    # Remember that you lost precision when concatenating arr over using a list.
    v0 = np.concatenate((t0a0, fluids, de0z0dl0))

    # Extracting the parsed mode of interaction.
    firstderivs_function = firstderivs_functions.get(model,0)
    assert firstderivs_function != 0, "zodesolve doesn't have this firstderivs_key at the top"

    # Call the ODE solver.
    vsol = odeint(firstderivs_function, v0, zpicks, args=(int_terms,H0), mxstep=5000000, atol=1.0e-8, rtol=1.0e-6)
    z = vsol[1:,-2]
    dl = vsol[1:,-1] * (1+z)        # in units of dl*(H0/c)
    da = dl * (1.0+z)**(-2.0)       # in units of dl*(H0/c)
    dlpc = dl * c_over_H0           # dl in parsecs (= vsol[dl] * c/H0)
    dapc = dlpc * (1.0+z)**(-2.0)   # in units of pc
    dapc = dapc / 10**3             # in units of kpc
#    theta = dl/da
#    print('theta',theta[-1])
#    print(model,'redshift = ',z[-1],'da =',da[-1])
#    plt.figure()
#    plt.title('Angular diameter distance vs redshift')
#    plt.xlabel('z')
#    plt.ylabel('D_A')
#    plt.plot(z, da, label='da')
#    plt.legend()
#
#    plt.figure()
#    plt.title('Angular diameter vs redshift')
#    plt.xlabel('z')
#    plt.ylabel(r'$\theta$')
#    plt.plot(z, theta, label='theta')
#    plt.legend()
#    plt.show()

    plot_var = {}
    if plot_key:
        # Separate results into their own arrays:
        plot_var['t'] = vsol[1:,0]
        plot_var['a'] = vsol[1:,1]

        # Collecting fluids and their names for plotting:
        fluid_arr = np.zeros(((int_in), (len(zpicks)-1)))
        fluid_names = []
        for i in range((int_in-1)):
            fluid_names.append(names[i+1])
            fluid_arr[i] = vsol[1:,(i+2)]
        fluid_names.append('de_ombar')
        fluid_arr[-1] = vsol[1:,-3]
        plot_var['fluid_names'] = fluid_names
        plot_var['fluid_arr'] = fluid_arr

        plot_var['z'] = z
        plot_var['dl'] = dl # in units of dl*(H0/c)
        plot_var['int_terms'] = int_terms

        plot_var['da'] = da

#        plt.figure()
#        plt.title('Angular diameter distance evolution')
#        plt.xlabel('z')
#        plt.ylabel(r'$ \left( \frac{H_0}{c} \right) d_A $', fontsize=15, labelpad=10)
#        plt.plot(z, da)

        Hz = H0 * (np.sum(fluid_arr, axis=0))**(0.5)
        plot_var['Hz'] = Hz

        daMpc = dlpc/10**6 * (1.0+z)**(-2.0) # in units of dl in Mpc*(H0/c)
        dV = (daMpc**2 * c*z/Hz)**(1/3) # combines radial and transverse dilation
        plot_var['dV'] = dV
#        plt.figure()
#        plt.title(r'$d_V$ evolution')
#        plt.xlabel('z')
#        plt.ylabel(r'$ d_V (z)$ [Mpc]')
#        plt.grid(True)
#        plt.plot(z, dV)


#        Dv = ((1+z)**2 * daMpc**2 * c*z/Hz)**(1/3)
#        plt.figure()
#        plt.title(r'$D_v$ evolution')
#        plt.xlabel('z')
#        plt.ylabel(r'$ D_v (z)$ [Mpc]')
#        plt.grid(True)
#        plt.plot(z, Dv)

#        # Calculating the sound horizon
#        ombar_m = vsol[1:,2]
#        ombar_baryon =  ombar_m*0.04 #0.0125
##        ombar_baryon = values[(fluid_in-1)]
#        s = 44.5 * np.log(9.83 / ombar_m) / (1 +10 * ombar_baryon**(3/4))**(1/2)
#        plt.figure()
#        plt.title('Sound horizon')
#        plt.xlabel(r'$z$')
#        plt.ylabel('Physical length in Mpc')
#        plt.grid(True)
#        plt.plot(z, s, label=r'$s_H$')
#        plt.legend()

        plt.show()
    return dlpc, da, plot_var