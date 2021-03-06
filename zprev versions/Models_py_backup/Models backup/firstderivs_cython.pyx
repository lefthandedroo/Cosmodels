#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:23:25 2018

@author: BallBlueMeercat
"""
import numpy as np
from cpython cimport array
import array
import math

# 'exotic'
# 'late_intxde'
# 'heaviside_late_int'
# 'late_int'
# 'expgamma'
# 'txgamma'
# 'zxgamma'
# 'gamma_over_z'
# 'zxxgamma'
# 'gammaxxz'
# 'rdecay_m'
# 'rdecay_de'
# 'rdecay_mxde'
# 'rdecay'                        
# 'interacting'
# 'LCDM'

# Eq of state parameters for known fluids:
cdef double w_r = 1/3     # radiation
cdef double w_m = 0.0     # matter
cdef double w_de = -1.0   # cosmological constant (dark energy?)

def exotic(double[:] v, redshifts, double gamma, double zeta, double H0):
    """    
    matter decays into radiation, that decays into dark energy
    Takes in:
        v = values at z=0;
        t = list of redshifts to integrate over;
        gamma = float interaction term, matter decaying into radioation;
        zeta = float interaction term, radioation decaying into dark energy.
                
    Returns a function f =     [dt/dz, d(a)/dz, 
                                d(e'_m)/dz, d(e'_de)/dz, 
                                d(z)/dz,
                                d(dl)/dz]
    """
    cdef double t = v[0]
    cdef double a = v[1]
    cdef double ombar_m = v[2]
    cdef double ombar_r = v[3]
    cdef double ombar_de = v[4]
    cdef double z = v[5]
    cdef double dl = v[6]
        
    cdef double Hz = H0 * (ombar_m + ombar_de + ombar_r)**(0.5)
        
    if math.isnan(Hz):
        print('exotic')
        print('z = %s, Hz = %s, gamma = %s, zeta = %s'% (z, Hz, gamma, zeta))
        print('ombar_m = %s, ombar_r = %s, ombar_de = %s'
              % (ombar_m, ombar_r, ombar_de))
    
    cdef double dtdz = -1.0/((1.0+z) * Hz)
    cdef double dadz = -(1.0+z)**(-2.0)
    cdef double domdz = 3.0*ombar_m /(1.0+z) +gamma * ombar_m /(1.0+z) /Hz
    cdef double dordz = 4.0*ombar_r /(1.0+z) -gamma * ombar_m /(1.0+z) /Hz +zeta * ombar_r /(1.0+z) /Hz
    cdef double dodedz = -zeta * ombar_r /(1.0+z) /Hz
    cdef double ddldz = 1.0/Hz
    
    # first derivatives of functions I want to find:
    f = [dtdz,# dt/dz (= f.d wrt z of time)
         dadz,# d(a)/dz (= f.d wrt z of scale factor)
         domdz,# d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         dordz,# d(ombar_r)/dz   (= f.d wrt z of density_r(t) / crit density(t0))
         dodedz,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)
    
    return f


def late_intxde(double[:] v, redshifts, double gamma, double H0):
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
    cdef double t = v[0]
    cdef double a = v[1]
    cdef double ombar_m = v[2]
    cdef double ombar_de = v[3]
    cdef double z = v[4]
    cdef double dl = v[5]
        
    cdef double Hz = H0 * (ombar_m + ombar_de)**(0.5)
        
    if math.isnan(Hz):
        print('late_intxde')
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))

    irate = 0.0

    if z < 0.9:
        irate = ombar_de * gamma/(1.0+z)/Hz        
    
    cdef double dtdz = -1.0/((1.0+z) * Hz)
    cdef double dadz = -(1.0+z)**(-2.0)
    cdef double domdz = 3.0*ombar_m /(1.0+z) - irate
    cdef double ddldz = 1.0/Hz
    
    # first derivatives of functions I want to find:
    f = [dtdz,# dt/dz (= f.d wrt z of time)
         dadz,# d(a)/dz (= f.d wrt z of scale factor)
         domdz,# d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)
    
    return f

def heaviside_late_int(double[:] v, redshifts, double gamma, double H0):
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
    cdef double t = v[0]
    cdef double a = v[1]
    cdef double ombar_m = v[2]
    cdef double ombar_de = v[3]
    cdef double z = v[4]
    cdef double dl = v[5]
        
    cdef double Hz = H0 * (ombar_m + ombar_de)**(0.5)
        
    if math.isnan(Hz):
        print('late_int')
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))
    
    irate = gamma/(1.0+z)/Hz * np.heaviside(0.9-z, 1)     
    
    cdef double dtdz = -1.0/((1.0+z) * Hz)
    cdef double dadz = -(1.0+z)**(-2.0)
    cdef double domdz = 3.0*ombar_m /(1.0+z) - irate
    cdef double ddldz = 1.0/Hz
    
    # first derivatives of functions I want to find:
    f = [dtdz,# dt/dz (= f.d wrt z of time)
         dadz,# d(a)/dz (= f.d wrt z of scale factor)
         domdz,# d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)
    
    return f

def late_int(double[:] v, redshifts, double gamma, double H0):
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
    cdef double t = v[0]
    cdef double a = v[1]
    cdef double ombar_m = v[2]
    cdef double ombar_de = v[3]
    cdef double z = v[4]
    cdef double dl = v[5]
        
    cdef double Hz = H0 * (ombar_m + ombar_de)**(0.5)
        
    if math.isnan(Hz):
        print('late_int')
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))

    irate = 0.0

    if z < 0.9:
        irate = gamma/(1.0+z)/Hz        
    
    cdef double dtdz = -1.0/((1.0+z) * Hz)
    cdef double dadz = -(1.0+z)**(-2.0)
    cdef double domdz = 3.0*ombar_m /(1.0+z) - irate
    cdef double ddldz = 1.0/Hz
    
    # first derivatives of functions I want to find:
    f = [dtdz,# dt/dz (= f.d wrt z of time)
         dadz,# d(a)/dz (= f.d wrt z of scale factor)
         domdz,# d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)
    
    return f

def expgamma(double[:] v, redshifts, double gamma, double H0):
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
    cdef double t = v[0]
    cdef double a = v[1]
    cdef double ombar_m = v[2]
    cdef double ombar_de = v[3]
    cdef double z = v[4]
    cdef double dl = v[5]
        
    cdef double Hz = H0 * (ombar_m + ombar_de)**(0.5)
        
    if math.isnan(Hz):
        print('expgamma')
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))
        
    cdef double irate = (1.0-math.exp(gamma))*(1.0-ombar_de/(ombar_de+ombar_m)) /(1.0+z)/Hz
    
    cdef double dtdz = -1.0/((1.0+z) * Hz)
    cdef double dadz = -(1.0+z)**(-2.0)
    cdef double domdz = 3.0*ombar_m /(1.0+z) - irate
    cdef double ddldz = 1.0/Hz
    
    # first derivatives of functions I want to find:
    f = [dtdz,# dt/dz (= f.d wrt z of time)
         dadz,# d(a)/dz (= f.d wrt z of scale factor)
         domdz,# d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)
    
    return f

def txgamma(double[:] v, redshifts, double gamma, double H0):
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
    cdef double t = v[0]
    cdef double a = v[1]
    cdef double ombar_m = v[2]
    cdef double ombar_de = v[3]
    cdef double z = v[4]
    cdef double dl = v[5]
        
    cdef double Hz = H0 * (ombar_m + ombar_de)**(0.5)
        
    if math.isnan(Hz):
        print('txgamma')
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))
        
    cdef double irate = (gamma/(-t+0.0001))*(1.0-ombar_de/(ombar_de+ombar_m)) /(1.0+z)/Hz
    
    cdef double dtdz = -1.0/((1.0+z) * Hz)
    cdef double dadz = -(1.0+z)**(-2.0)
    cdef double domdz = 3.0*ombar_m /(1.0+z) - irate
    cdef double ddldz = 1.0/Hz
    
    # first derivatives of functions I want to find:
    f = [dtdz,# dt/dz (= f.d wrt z of time)
         dadz,# d(a)/dz (= f.d wrt z of scale factor)
         domdz,# d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)
    
    return f

def zxgamma(double[:] v, redshifts, double gamma, double H0):
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
    cdef double t = v[0]
    cdef double a = v[1]
    cdef double ombar_m = v[2]
    cdef double ombar_de = v[3]
    cdef double z = v[4]
    cdef double dl = v[5]
        
    cdef double Hz = H0 * (ombar_m + ombar_de)**(0.5)
        
    if math.isnan(Hz):
        print('zxgamma')
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))
        
    cdef double irate = z*gamma*(1.0-ombar_de/(ombar_de+ombar_m)) /(1.0+z)/Hz
    
    cdef double dtdz = -1.0/((1.0+z) * Hz)
    cdef double dadz = -(1.0+z)**(-2.0)
    cdef double domdz = 3.0*ombar_m /(1.0+z) - irate
    cdef double ddldz = 1.0/Hz
    
    # first derivatives of functions I want to find:
    f = [dtdz,# dt/dz (= f.d wrt z of time)
         dadz,# d(a)/dz (= f.d wrt z of scale factor)
         domdz,# d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)
    
    return f

def gamma_over_z(double[:] v, redshifts, double gamma, double H0):
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
    cdef double t = v[0]
    cdef double a = v[1]
    cdef double ombar_m = v[2]
    cdef double ombar_de = v[3]
    cdef double z = v[4]
    cdef double dl = v[5]
        
    cdef double Hz = H0 * (ombar_m + ombar_de)**(0.5)
        
    if math.isnan(Hz):
        print('gamma_over_z')
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))
        
    cdef double irate = gamma/(z + 0.01)*(1.0-ombar_de/(ombar_de+ombar_m)) /(1.0+z)/Hz
    
    cdef double dtdz = -1.0/((1.0+z) * Hz)
    cdef double dadz = -(1.0+z)**(-2.0)
    cdef double domdz = 3.0*ombar_m /(1.0+z) - irate
    cdef double ddldz = 1.0/Hz
    
    # first derivatives of functions I want to find:
    f = [dtdz,# dt/dz (= f.d wrt z of time)
         dadz,# d(a)/dz (= f.d wrt z of scale factor)
         domdz,# d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)
    
    return f

def zxxgamma(double[:] v, redshifts, double gamma, double H0):
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
    cdef double t = v[0]
    cdef double a = v[1]
    cdef double ombar_m = v[2]
    cdef double ombar_de = v[3]
    cdef double z = v[4]
    cdef double dl = v[5]
        
    cdef double Hz = H0 * (ombar_m + ombar_de)**(0.5)
    
    gamma = abs(gamma)
    
    if math.isnan(Hz):
        print('zxxgamma')
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))
        
    cdef double irate = (z**gamma)*(1.0-ombar_de/(ombar_de+ombar_m)) /(1.0+z)/Hz
    
    cdef double dtdz = -1.0/((1.0+z) * Hz)
    cdef double dadz = -(1.0+z)**(-2.0)
    cdef double domdz = 3.0*ombar_m /(1.0+z) - irate
    cdef double ddldz = 1.0/Hz
    
    # first derivatives of functions I want to find:
    f = [dtdz,# dt/dz (= f.d wrt z of time)
         dadz,# d(a)/dz (= f.d wrt z of scale factor)
         domdz,# d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)
    
    return f

def gammaxxz(double[:] v, redshifts, double gamma, double H0):
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
    cdef double t = v[0]
    cdef double a = v[1]
    cdef double ombar_m = v[2]
    cdef double ombar_de = v[3]
    cdef double z = v[4]
    cdef double dl = v[5]
        
    cdef double Hz = H0 * (ombar_m + ombar_de)**(0.5)
    
    gamma = abs(gamma)
        
    if math.isnan(Hz):
        print('gammaxxz')
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))
        
    cdef double irate = (gamma**z)*(1.0-ombar_de/(ombar_de+ombar_m)) /(1.0+z)/Hz
    
    cdef double dtdz = -1.0/((1.0+z) * Hz)
    cdef double dadz = -(1.0+z)**(-2.0)
    cdef double domdz = 3.0*ombar_m /(1.0+z) - irate
    cdef double ddldz = 1.0/Hz
    
    # first derivatives of functions I want to find:
    f = [dtdz,# dt/dz (= f.d wrt z of time)
         dadz,# d(a)/dz (= f.d wrt z of scale factor)
         domdz,# d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)
    
    return f

def rdecay_m(double[:] v, redshifts, double gamma, double H0):
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
    cdef double t = v[0]
    cdef double a = v[1]
    cdef double ombar_m = v[2]
    cdef double ombar_de = v[3]
    cdef double z = v[4]
    cdef double dl = v[5]
        
    cdef double Hz = H0 * (ombar_m + ombar_de)**(0.5)
        
    if math.isnan(Hz):
        print('rdecay_m')
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))
        
    cdef double irate = gamma*(1.0-ombar_m/(ombar_de+ombar_m)) /(1.0+z)/Hz
    
    cdef double dtdz = -1.0/((1.0+z) * Hz)
    cdef double dadz = -(1.0+z)**(-2.0)
    cdef double domdz = 3.0*ombar_m /(1.0+z) - irate
    cdef double ddldz = 1.0/Hz
    
    # first derivatives of functions I want to find:
    f = [dtdz,# dt/dz (= f.d wrt z of time)
         dadz,# d(a)/dz (= f.d wrt z of scale factor)
         domdz,# d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)
    
    return f

def rdecay_de(double[:] v, redshifts, double gamma, double H0):
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
    cdef double t = v[0]
    cdef double a = v[1]
    cdef double ombar_m = v[2]
    cdef double ombar_de = v[3]
    cdef double z = v[4]
    cdef double dl = v[5]
        
    cdef double Hz = H0 * (ombar_m + ombar_de)**(0.5)
        
    if math.isnan(Hz):
        print('rdecay_de')
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))
        
    cdef double irate = gamma*(1.0-ombar_de/(ombar_de+ombar_m)) /(1.0+z)/Hz
    
    cdef double dtdz = -1.0/((1.0+z) * Hz)
    cdef double dadz = -(1.0+z)**(-2.0)
    cdef double domdz = 3.0*ombar_m /(1.0+z) - irate
    cdef double ddldz = 1.0/Hz
    
    # first derivatives of functions I want to find:
    f = [dtdz,# dt/dz (= f.d wrt z of time)
         dadz,# d(a)/dz (= f.d wrt z of scale factor)
         domdz,# d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)
    
    return f

def rdecay_mxde(double[:] v, redshifts, double gamma, double H0):
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
    cdef double t = v[0]
    cdef double a = v[1]
    cdef double ombar_m = v[2]
    cdef double ombar_de = v[3]
    cdef double z = v[4]
    cdef double dl = v[5]
        
    cdef double Hz = H0 * (ombar_m + ombar_de)**(0.5)
        
    if math.isnan(Hz):
        print('rdecay_mxde')
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))
        
    cdef double irate = gamma*ombar_de*ombar_m /(1.0+z)/Hz
    
    cdef double dtdz = -1.0/((1.0+z) * Hz)
    cdef double dadz = -(1.0+z)**(-2.0)
    cdef double domdz = 3.0*ombar_m /(1.0+z) - irate
    cdef double ddldz = 1.0/Hz
    
    # first derivatives of functions I want to find:
    f = [dtdz,# dt/dz (= f.d wrt z of time)
         dadz,# d(a)/dz (= f.d wrt z of scale factor)
         domdz,# d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)
    
    return f

def rdecay(double[:] v, redshifts, double gamma, double H0):
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
    cdef double t = v[0]
    cdef double a = v[1]
    cdef double ombar_m = v[2]
    cdef double ombar_de = v[3]
    cdef double z = v[4]
    cdef double dl = v[5]
        
    cdef double Hz = H0 * (ombar_m + ombar_de)**(0.5)
        
    if math.isnan(Hz):
        print('rdecay')
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))
        
    cdef double irate = gamma*ombar_de /(1.0+z)/Hz
    
    cdef double dtdz = -1.0/((1.0+z) * Hz)
    cdef double dadz = -(1.0+z)**(-2.0)
    cdef double domdz = 3.0*ombar_m /(1.0+z) - irate
    cdef double ddldz = 1.0/Hz
    
    # first derivatives of functions I want to find:
    f = [dtdz,# dt/dz (= f.d wrt z of time)
         dadz,# d(a)/dz (= f.d wrt z of scale factor)
         domdz,# d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)
    
    return f

def interacting(double[:] v, redshifts, double gamma, double H0):
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
    cdef double t = v[0]
    cdef double a = v[1]
    cdef double ombar_m = v[2]
    cdef double ombar_de = v[3]
    cdef double z = v[4]
    cdef double dl = v[5]
        
    cdef double Hz = H0 * (ombar_m + ombar_de)**(0.5)
        
    if math.isnan(Hz):
        print('interacting')
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))
        
    cdef double irate = gamma/(1.0+z)/Hz
    
    cdef double dtdz = -1.0/((1.0+z) * Hz)
    cdef double dadz = -(1.0+z)**(-2.0)
    cdef double domdz = 3.0*ombar_m /(1.0+z) - irate
    cdef double ddldz = 1.0/Hz
    
    # first derivatives of functions I want to find:
    f = [dtdz,# dt/dz (= f.d wrt z of time)
         dadz,# d(a)/dz (= f.d wrt z of scale factor)
         domdz,# d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)
    
    return f

def LCDM(double[:] v, redshifts, double gamma, double H0):
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
    cdef double t = v[0]
    cdef double a = v[1]
    cdef double ombar_m = v[2]
    cdef double ombar_de = v[3]
    cdef double z = v[4]
    cdef double dl = v[5]
        
    cdef double Hz = H0 * (ombar_m + ombar_de)**(0.5)
        
    if math.isnan(Hz):
        print('LCDM')                                        
        print('z = %s, Hz = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, ombar_m, ombar_de))
        
    cdef double dtdz = -1.0/((1.0+z) * Hz)
    cdef double dadz = -(1.0+z)**(-2.0)
    cdef double domdz = 3.0*ombar_m /(1.0+z)
    cdef double ddldz = 1.0/Hz
    
    # first derivatives of functions I want to find:
    f = [dtdz,# dt/dz (= f.d wrt z of time)
         dadz,# d(a)/dz (= f.d wrt z of scale factor)
         domdz,# d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         0.0,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit desnity(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)
        
    return f