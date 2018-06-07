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

def edecay(v, t, gamma, H0):
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
        
    import numpy as np
    if np.isnan(Hz):
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))

    # fist derivatives of functions I want to find:
    f = [# dt/dz (= f.d wrt z of time)
        -1/((1+z) * Hz),
            
        # d(a)/dz (= f.d wrt z of scale factor)
         -(1+z)**(-2),
         
#         # d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
#         3*ombar_m /(1+z) - gamma*np.exp(t)/(1+z)/Hz,
#         
#         # d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
#         gamma*ombar_de*ombar_m/(1+z)/Hz,
         
         # d(z)/dz (= f.d wrt z of redshift)
         1,
         
         # d(dl)/dz (= f.d wrt z of luminosty distance)
         1/Hz] # H + Hdz*(1+z)
        
    return f

def Hdecay(v, t, gamma, H0):
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
        
    import numpy as np
    if np.isnan(Hz):
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))

    # fist derivatives of functions I want to find:
    f = [# dt/dz (= f.d wrt z of time)
        -1/((1+z) * Hz),
            
        # d(a)/dz (= f.d wrt z of scale factor)
         -(1+z)**(-2),
         
#         # d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
#         3*ombar_m /(1+z) - gamma*ombar_de*ombar_m /(1+z)/Hz,
#         
#         # d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
#         gamma*ombar_de*ombar_m/(1+z)/Hz,
         
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
        gamma = interaction term.
                
    Returns a function f =     [dt/dz, d(a)/dz, 
                                d(e'_m)/dz, d(e'_de)/dz, 
                                d(z)/dz,
                                d(dl)/dz]
    """
    (t, a, ombar_m, ombar_de, z, dl) = v #omegam, omegade, z, dl) = v

    Hz = H0 * (ombar_m + ombar_de)**(1/2)
        
    import numpy as np
    if np.isnan(Hz):
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))

    # fist derivatives of functions I want to find:
    f = [# dt/dz (= f.d wrt z of time)
        -1/((1+z) * Hz),
            
        # d(a)/dz (= f.d wrt z of scale factor)
         -(1+z)**(-2),
         
         # d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         3*ombar_m /(1+z) - gamma*ombar_de*ombar_m /(1+z)/Hz,
         
         # d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         gamma*ombar_de*ombar_m/(1+z)/Hz,
         
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
        gamma = interaction term.
                
    Returns a function f =     [dt/dz, d(a)/dz, 
                                d(e'_m)/dz, d(e'_de)/dz, 
                                d(z)/dz,
                                d(dl)/dz]
    """
    (t, a, ombar_m, ombar_de, z, dl) = v #omegam, omegade, z, dl) = v

    Hz = H0 * (ombar_m + ombar_de)**(1/2)
        
    import numpy as np
    if np.isnan(Hz):
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))

    # fist derivatives of functions I want to find:
    f = [# dt/dz (= f.d wrt z of time)
        -1/((1+z) * Hz),
            
        # d(a)/dz (= f.d wrt z of scale factor)
         -(1+z)**(-2),
         
         # d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         3*ombar_m /(1+z) - gamma*ombar_de /(1+z)/Hz,
         
         # d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         gamma*ombar_de*ombar_m/(1+z)/Hz,
         
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
        gamma = interaction term.
                
    Returns a function f =     [dt/dz, d(a)/dz, 
                                d(e'_m)/dz, d(e'_de)/dz, 
                                d(z)/dz,
                                d(dl)/dz]
    """
    (t, a, ombar_m, ombar_de, z, dl) = v #omegam, omegade, z, dl) = v

    Hz = H0 * (ombar_m + ombar_de)**(1/2)
        
    import numpy as np
    if np.isnan(Hz):
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))

    # fist derivatives of functions I want to find:
    f = [# dt/dz (= f.d wrt z of time)
        -1/((1+z) * Hz),
            
        # d(a)/dz (= f.d wrt z of scale factor)
         -(1+z)**(-2),
         
         # d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         3*ombar_m /(1+z) - gamma*ombar_m /(1+z)/Hz,
         
         # d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         gamma*ombar_de*ombar_m/(1+z)/Hz,
         
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
        
    import numpy as np
    if np.isnan(Hz):
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))

    # fist derivatives of functions I want to find:
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
        
    import numpy as np
    if np.isnan(Hz):
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, ombar_m, ombar_de))

    # fist derivatives of functions I want to find:
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