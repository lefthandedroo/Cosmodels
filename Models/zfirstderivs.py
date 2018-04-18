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
        t = list of redshifts to integrate over;
        gamma = interaction term.
                
    Returns a function f =     [dt/dz, d(a)/dz, 
                                d(e'_m)/dz, d(e'_de)/dz, 
                                d(z)/dz,
                                d(dl)/dz]
    """
    (t, a, e_dashm, e_dashde, z, dl) = v #omegam, omegade, z, dl) = v
    
    Hz = (e_dashm+e_dashde)**(1/2)
    
    import numpy as np
    if np.isnan(Hz):
        print("z = %s, Hz = %s, e'_m = %s, e'_de = %s"%(z, Hz, e_dashm, e_dashde))
    
    # fist derivatives of functions I want to find:
    f = [# dt/dz (= f.d wrt z of time)
        -1/(1+z)/Hz,
            
        # d(a)/dz (= f.d wrt z of scale factor)
         -(1+z)**(-2),
         
         # d(e'_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         3*e_dashm /(1+z) - gamma/(1+z)/Hz,
         
         # d(e'_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         gamma/(1+z)/Hz,
         
         # d(z)/dz (= f.d wrt z of redshift)
         1,
         
         # d(dl)/dz (= f.d wrt z of luminosty distance)
         1/Hz] # H + Hdz*(1+z)
        
    return f