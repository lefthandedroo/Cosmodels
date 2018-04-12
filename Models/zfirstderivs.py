#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:23:25 2018

@author: BallBlueMeercat
"""
# Eq of state parameters for known fluids:
w_r = 1/3     # radiation
w_m = 0.0     # matter
w_de = -1.0   # cosmological constant (dark energy?)

def zfirstderivs(v, t, gamma):
    """
    Takes in:
        v = values at z=0;
        t = list of times to be integrated over;
        w = omega parameters;
        gamma = interaction constant.
                
    Returns a function with:    a_dot, a_dotdot, 
                                e'_dotm, e'_dotde, 
                                omegam_dot, omegade_dot,
                                z_dot,
                                dl_dot
    ready to be integrated with odeint.
    Uses same lambda for all fluids.
    """
    (t, a, e_dashm, e_dashde, z, dl) = v #omegam, omegade, z, dl) = v
    
    Hz = (e_dashm+e_dashde)**(1/2)

    
    # fist derivatives of functions I want to find:
    f = [# dt/dz (= time)
        -1/(1+z)/Hz,
            
        # d(a)/dz (= f.d. of scale factor)
         -(1+z)**(-2),
         
         # d(e'_m)/dz   (= f.d. of density(t) / crit density(t0))
         3*e_dashm /(1+z) - gamma/(1+z)/Hz,
         
         # d(e'_de)/dz
         gamma/(1+z)/Hz,
         
         # d(z)/dz (= f.d. of redshift)
         1,
         
         # d(dl)/dz (= f.d. of luminosty distance)
         1/Hz] # H + Hdz*(1+z)
        
    return f