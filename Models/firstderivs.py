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

def firstderivs(v, t, gamma):
    """
    Takes in:
        v = values at t=0;
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
    (a, a_dot, e_dashm, e_dashde, z, dl) = v #omegam, omegade, z, dl) = v
    
    # fist derivatives of functions I want to find:
    f = [# a_dot (=scale factor)
         a_dot,
         
         # a_dotdot
         (-a/2) * (e_dashm * (1+3*w_m) + e_dashde * (1+3*w_de)), 
         
         # e'_dotm (=density(t) / crit density(t0))
#         -3 * (a_dot/a) * e_dashm * (1 + w_m -gamma/3 * a/a_dot * e_dashde/e_dashm), # old version
#         -3 * (a_dot/a) * (e_dashm + e_dashm * w_m + e_dashm * gamma/3 * e_dashde),
         gamma -3 * (a_dot/a) * e_dashm,
         
         # e'_dotde
#        -3 * (a_dot/a) * e_dashde * (1 + w_de +lamb/3 * a/a_dot), # old version
         -gamma,

         # z_dot (=redshift)
         -a_dot/a**2,
         
         # dl_dot (=luminosty distance)
         -1/a]
        
    return f