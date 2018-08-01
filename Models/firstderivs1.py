#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:23:25 2018

@author: BallBlueMeercat
"""
import numpy as np

# Eq of state parameters for known fluids:
w_r = 1/3     # radiation
w_m = 0.0     # matter
w_de = -1.0   # cosmological constant (dark energy?)

def expgamma(v, t, gamma, H0):
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
    (t, a, ombar_m, ombar_de, z, dl) = v

    Hz = H0 * (ombar_m + ombar_de)**(1/2)
        
    if np.isnan(Hz):
        print('expgamma')
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))

    irate = (1-np.exp(gamma))*(1-ombar_de/(ombar_de+ombar_m)) /(1+z)/Hz

    
    # first derivatives of functions I want to find:
    f = [# dt/dz (= f.d wrt z of time)
        -1/((1+z) * Hz),
            
        # d(a)/dz (= f.d wrt z of scale factor)
         -(1+z)**(-2),
         
         # d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         3*ombar_m /(1+z) - irate,
         
         # d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         irate,
         
         # d(z)/dz (= f.d wrt z of redshift)
         1,
         
         # d(dl)/dz (= f.d wrt z of luminosty distance)
         1/Hz] # H + Hdz*(1+z)
        
    return f


def txgamma(v, t, gamma, H0):
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
    (t, a, ombar_m, ombar_de, z, dl) = v #omegam, omegade, z, dl) = v

    Hz = H0 * (ombar_m + ombar_de)**(1/2)
        
    if np.isnan(Hz):
        print('txgamma')
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))
        
    irate = (gamma/(-t+0.0001))*(1-ombar_de/(ombar_de+ombar_m)) /(1+z)/Hz

    # first derivatives of functions I want to find:
    f = [# dt/dz (= f.d wrt z of time)
        -1/((1+z) * Hz),
            
        # d(a)/dz (= f.d wrt z of scale factor)
         -(1+z)**(-2),
         
         # d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         3*ombar_m /(1+z) - irate,
         
         # d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         irate,
         
         # d(z)/dz (= f.d wrt z of redshift)
         1,
         
         # d(dl)/dz (= f.d wrt z of luminosty distance)
         1/Hz] # H + Hdz*(1+z)
        
    return f

def zxgamma(v, t, gamma, H0):
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
    (t, a, ombar_m, ombar_de, z, dl) = v #omegam, omegade, z, dl) = v

    Hz = H0 * (ombar_m + ombar_de)**(1/2)
        
    if np.isnan(Hz):
        print('zxgamma')
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))

    irate = z*gamma*(1-ombar_de/(ombar_de+ombar_m)) /(1+z)/Hz
    
    # first derivatives of functions I want to find:
    f = [# dt/dz (= f.d wrt z of time)
        -1/((1+z) * Hz),
            
        # d(a)/dz (= f.d wrt z of scale factor)
         -(1+z)**(-2),
         
         # d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         3*ombar_m /(1+z) - irate,
         
         # d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         irate,
         
         # d(z)/dz (= f.d wrt z of redshift)
         1,
         
         # d(dl)/dz (= f.d wrt z of luminosty distance)
         1/Hz] # H + Hdz*(1+z)
        
    return f


def gamma_over_z(v, t, gamma, H0):
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
    (t, a, ombar_m, ombar_de, z, dl) = v #omegam, omegade, z, dl) = v

    Hz = H0 * (ombar_m + ombar_de)**(1/2)
        
    if np.isnan(Hz):
        print('gamma_over_z')        
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))

    irate = gamma/(z + 0.01)*(1-ombar_de/(ombar_de+ombar_m)) /(1+z)/Hz
    
    # first derivatives of functions I want to find:
    f = [# dt/dz (= f.d wrt z of time)
        -1/((1+z) * Hz),
            
        # d(a)/dz (= f.d wrt z of scale factor)
         -(1+z)**(-2),
         
         # d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         3*ombar_m /(1+z) - irate,
         
         # d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         irate,
         
         # d(z)/dz (= f.d wrt z of redshift)
         1,
         
         # d(dl)/dz (= f.d wrt z of luminosty distance)
         1/Hz] # H + Hdz*(1+z)
        
    return f


def zxxgamma(v, t, gamma, H0):
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
    (t, a, ombar_m, ombar_de, z, dl) = v #omegam, omegade, z, dl) = v

    Hz = H0 * (ombar_m + ombar_de)**(1/2)
        
    if np.isnan(Hz):
        print('zxxgamma')        
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))

    irate = (z**gamma)*(1-ombar_de/(ombar_de+ombar_m)) /(1+z)/Hz
    
    # first derivatives of functions I want to find:
    f = [# dt/dz (= f.d wrt z of time)
        -1/((1+z) * Hz),
            
        # d(a)/dz (= f.d wrt z of scale factor)
         -(1+z)**(-2),
         
         # d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         3*ombar_m /(1+z) - irate,
         
         # d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         irate,
         
         # d(z)/dz (= f.d wrt z of redshift)
         1,
         
         # d(dl)/dz (= f.d wrt z of luminosty distance)
         1/Hz] # H + Hdz*(1+z)
        
    return f


def gammaxxz(v, t, gamma, H0):
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
    (t, a, ombar_m, ombar_de, z, dl) = v #omegam, omegade, z, dl) = v

    Hz = H0 * (ombar_m + ombar_de)**(1/2)
        
    if np.isnan(Hz):
        print('gammaxxz')                
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))

    irate = (gamma**z)*(1-ombar_de/(ombar_de+ombar_m)) /(1+z)/Hz
    
    # first derivatives of functions I want to find:
    f = [# dt/dz (= f.d wrt z of time)
        -1/((1+z) * Hz),
            
        # d(a)/dz (= f.d wrt z of scale factor)
         -(1+z)**(-2),
         
         # d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         3*ombar_m /(1+z) - irate,
         
         # d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         irate,
         
         # d(z)/dz (= f.d wrt z of redshift)
         1,
         
         # d(dl)/dz (= f.d wrt z of luminosty distance)
         1/Hz] # H + Hdz*(1+z)
        
    return f


def rdecay_m(v, t, gamma, H0):
    """
    Takes in:
        v = values at z=0;
        t = list of redshifts to integrate over;
        gamma = interaction term;
        H0 = Hubble constant ar z=0.
                
    Returns a function f =     [dt/dz, d(a)/dz, 
                                d(e'_m)/dz, d(e'_de)/dz, 
                                d(z)/dz,
                                d(dl)/dz]
    """
    (t, a, ombar_m, ombar_de, z, dl) = v #omegam, omegade, z, dl) = v

    Hz = H0 * (ombar_m + ombar_de)**(1/2)
        
    if np.isnan(Hz):
        print('rdecay_m')                        
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))
    
    # rate of ombar change with redshift
    irate = gamma*(1-ombar_m/(ombar_de+ombar_m)) /(1+z)/Hz
    
    # first derivatives of functions I want to find:
    f = [# dt/dz (= f.d wrt z of time)
        -1/((1+z) * Hz),
            
        # d(a)/dz (= f.d wrt z of scale factor)
         -(1+z)**(-2),
         
         # d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         3*ombar_m /(1+z) - irate,
         
         # d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         irate,
         
         # d(z)/dz (= f.d wrt z of redshift)
         1,
         
         # d(dl)/dz (= f.d wrt z of luminosty distance)
         1/Hz] # H + Hdz*(1+z)
    
    return f


def rdecay_de(v, t, gamma, H0):
    """
    Takes in:
        v = values at z=0;
        t = list of redshifts to integrate over;
        gamma = interaction term;
        H0 = Hubble constant ar z=0.
                
    Returns a function f =     [dt/dz, d(a)/dz, 
                                d(e'_m)/dz, d(e'_de)/dz, 
                                d(z)/dz,
                                d(dl)/dz]
    """
    (t, a, ombar_m, ombar_de, z, dl) = v #omegam, omegade, z, dl) = v

    Hz = H0 * (ombar_m + ombar_de)**(1/2)
        
    if np.isnan(Hz):
        print('rdecay_de')                                
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))
    
    # rate of ombar change with redshift
    irate = gamma*(1-ombar_de/(ombar_de+ombar_m)) /(1+z)/Hz
    
    # first derivatives of functions I want to find:
    f = [# dt/dz (= f.d wrt z of time)
        -1/((1+z) * Hz),
            
        # d(a)/dz (= f.d wrt z of scale factor)
         -(1+z)**(-2),
         
         # d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         3*ombar_m /(1+z) - irate,
         
         # d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         irate,
         
         # d(z)/dz (= f.d wrt z of redshift)
         1,
         
         # d(dl)/dz (= f.d wrt z of luminosty distance)
         1/Hz] # H + Hdz*(1+z)
        
    return f


def rdecay_mxde(v, t, gamma, H0):
    """
    Takes in:
        v = values at z=0;
        t = list of redshifts to integrate over;
        gamma = interaction term;
        H0 = Hubble constant ar z=0.
                
    Returns a function f =     [dt/dz, d(a)/dz, 
                                d(e'_m)/dz, d(e'_de)/dz, 
                                d(z)/dz,
                                d(dl)/dz]
    """
    (t, a, ombar_m, ombar_de, z, dl) = v #omegam, omegade, z, dl) = v

    Hz = H0 * (ombar_m + ombar_de)**(1/2)
        
    if np.isnan(Hz):
        print('rdecay_mxde')                                        
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))
    
    # rate of ombar change with redshift
    irate = gamma*ombar_de*ombar_m /(1+z)/Hz
    
    # first derivatives of functions I want to find:
    f = [# dt/dz (= f.d wrt z of time)
        -1/((1+z) * Hz),
            
        # d(a)/dz (= f.d wrt z of scale factor)
         -(1+z)**(-2),
         
         # d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         3*ombar_m /(1+z) - irate,
         
         # d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         irate,
         
         # d(z)/dz (= f.d wrt z of redshift)
         1,
         
         # d(dl)/dz (= f.d wrt z of luminosty distance)
         1/Hz] # H + Hdz*(1+z)
        
    return f
    

def rdecay(v, t, gamma, H0):
    """
    Takes in:
        v = values at z=0;
        t = list of redshifts to integrate over;
        gamma = interaction term;
        H0 = Hubble constant at z=0.
                
    Returns a function f =     [dt/dz, d(a)/dz, 
                                d(e'_m)/dz, d(e'_de)/dz, 
                                d(z)/dz,
                                d(dl)/dz]
    """
    (t, a, ombar_m, ombar_de, z, dl) = v #omegam, omegade, z, dl) = v

    Hz = H0 * (ombar_m + ombar_de)**(1/2)
        
    if np.isnan(Hz):
        print('rdecay')                                        
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))
    
    # rate of ombar change with redshift
    irate = gamma*ombar_de /(1+z)/Hz
    
    # first derivatives of functions I want to find:
    f = [# dt/dz (= f.d wrt z of time)
        -1/((1+z) * Hz),
            
        # d(a)/dz (= f.d wrt z of scale factor)
         -(1+z)**(-2),
         
         # d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         3*ombar_m /(1+z) - irate,
         
         # d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         irate,
         
         # d(z)/dz (= f.d wrt z of redshift)
         1,
         
         # d(dl)/dz (= f.d wrt z of luminosty distance)
         1/Hz] # H + Hdz*(1+z)
        
    return f


def interacting(v, t, gamma, H0):
    """
    UNPHYSICAL FOR |gamma| > 0.1 BEFORE z = 2
    
    Takes in:
        v = values at z=0;
        t = list of redshifts to integrate over;
        gamma = interaction term.
                
    Returns a function f =     [dt/dz, d(a)/dz, 
                                d(e'_m)/dz, d(e'_de)/dz, 
                                d(z)/dz,
                                d(dl)/dz]
    """
    (t, a, ombar_m, ombar_de, z, dl) = v #omegam, omegade, z, dl) = v

    Hz = H0 * (ombar_m + ombar_de)**(1/2)
        
    if np.isnan(Hz):
        print('interacting')                                        
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))

    # first derivatives of functions I want to find:
    f = [# dt/dz (= f.d wrt z of time)
        -1/((1+z) * Hz),
            
        # d(a)/dz (= f.d wrt z of scale factor)
         -(1+z)**(-2),
         
         # d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         3*ombar_m /(1+z) - gamma/(1+z)/Hz,
         
         # d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         gamma/(1+z)/Hz,
         
         # d(z)/dz (= f.d wrt z of redshift)
         1,
         
         # d(dl)/dz (= f.d wrt z of luminosty distance)
         1/Hz] # H + Hdz*(1+z)
        
    return f


def LCDM(v, t, H0):
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
    (t, a, ombar_m, ombar_de, z, dl) = v #omegam, omegade, z, dl) = v

    Hz = H0 * (ombar_m + ombar_de)**(1/2)
        
    if np.isnan(Hz):
        print('LCDM')                                        
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, ombar_m, ombar_de))

    # first derivatives of functions I want to find:
    f = [# dt/dz (= f.d wrt z of time)
        -1/((1+z) * Hz),
            
        # d(a)/dz (= f.d wrt z of scale factor)
         -(1+z)**(-2),
         
         # d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         3*ombar_m /(1+z),
         
         # d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         0,
         
         # d(z)/dz (= f.d wrt z of redshift)
         1,
         
         # d(dl)/dz (= f.d wrt z of luminosty distance)
         1/Hz] # H + Hdz*(1+z)
        
    return f


#def LCDM(v, z, H0):
#    """
#    Takes in:
#        v = values at z=0;
#        t = list of redshifts to integrate over;
#        gamma = interaction term.
#                
#    Returns a function f =     [dt/dz, d(a)/dz, 
#                                d(e'_m)/dz, d(e'_de)/dz, 
#                                d(z)/dz,
#                                d(dl)/dz]
#    """
#    (ombar_m, ombar_de, dl) = v #omegam, omegade, z, dl) = v
#
#    Hz = H0 * (ombar_m + ombar_de)**(1/2)
#        
#    import numpy as np
#    if np.isnan(Hz):
#        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
#              %(z, Hz, ombar_m, ombar_de))
#
#    # fist derivatives of functions I want to find:
#    f = [# d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
#         3*ombar_m /(1+z),
#         
#         # d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
#         0,
#         
#         # d(dl)/dz (= f.d wrt z of luminosty distance)
#         1/Hz] # H + Hdz*(1+z)
#        
#    return f