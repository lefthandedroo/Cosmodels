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
cdef double w_de = -1.0   # dark energy


def dlLCDM(double[:] v, redshifts, in_terms, double H0):
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
    
    cdef double da = dl * (1.0+z)**(-2.0)
    print('da = ',da)
    
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
         0.0,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit density(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)

    return f


def waterfall(double[:] v, redshifts, in_terms, double H0):
    """
    Takes in:
        v = list, values at z=0;
        redshifts = list, redshifts to integrate over;
        in_terms = list, rates of interaction;

    Returns a function f = [dt/dz, d(a)/dz, d(ombar_matter)/dz,
                            d(ombar_radiation)/dz, d(ombar_a)/dz,
                            d(ombar_b)/dz, d(ombar_c)/dz,
                            d(ombar_de)/dz, d(z)/dz, d(dl)/dz]
    """

    cdef double t = v[0]
    cdef double a = v[1]
    # ombar_i = rho_i(t) / rho_crit(t0)
    cdef double ombar_m = v[2]
    cdef double ombar_r = v[3]
    cdef double ombar_a = v[4]
    cdef double ombar_b = v[5]
    cdef double ombar_c = v[6]
    cdef double ombar_de = v[7]
    cdef double z = v[8]
    cdef double dl = v[9]

    cdef double in_v = in_terms[0]
    cdef double in_w = in_terms[1]
    cdef double in_x = in_terms[2]
    cdef double in_y = in_terms[3]
    cdef double in_z = in_terms[4]

    cdef double Hz = H0 * (ombar_r +ombar_m +ombar_a
                           +ombar_b +ombar_c +ombar_de)**(0.5)

    if math.isnan(Hz):
        print('waterfall')
        print('z = %s, Hz = %s, in_terms = %s'% (z, Hz, in_terms))
        print('ombar_r = ',ombar_r,'ombar_m = ',ombar_m,
              'ombar_a = ',ombar_a,'ombar_b = ',ombar_b,
              'ombar_c = ',ombar_c,'ombar_de = ',ombar_de)

    cdef double dtdz = -1.0/((1.0+z) * Hz)
    cdef double dadz = -(1.0+z)**(-2.0)
    # d(ombar_i)/dz
    cdef double domdz = (3.0*ombar_m +in_v*ombar_m*ombar_r/Hz)/(1.0+z)                          # w = 0
    cdef double dordz = (4.0*ombar_r -in_v*ombar_m*ombar_r/Hz +in_w*ombar_r*ombar_a/Hz)/(1.0+z) # w = 1/3
    cdef double doadz = (2.7*ombar_a -in_w*ombar_r*ombar_a/Hz +in_x*ombar_a*ombar_b/Hz)/(1.0+z) # w = -0.1
    cdef double dobdz = (1.5*ombar_b -in_x*ombar_a*ombar_b/Hz +in_y*ombar_b*ombar_c/Hz)/(1.0+z) # w = -0.5
    cdef double docdz = (0.6*ombar_c -in_y*ombar_b*ombar_c/Hz +in_z*ombar_c*ombar_de/Hz)/(1.0+z)# w = -0.8
    cdef double dodedz = -in_z*ombar_c*ombar_de/(1.0+z)/Hz                                      # w = -1
    cdef double ddldz = 1.0/Hz

    # first derivatives of functions I want to find:
    f = [dtdz,  # time
         dadz,  # scale factor
         domdz, # w = 0
         dordz, # w = 1/3
         doadz, # w = -0.1
         dobdz, # w = -0.5
         docdz, # w = -0.8
         dodedz,# w = -1
         1.0,   # redshift
         ddldz] # luminosty distance

    return f

def rainbow(double[:] v, redshifts, in_terms, double H0):
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
    cdef double ombar_a = v[4]
    cdef double ombar_b = v[5]
    cdef double ombar_c = v[6]
    cdef double ombar_d = v[7]
    cdef double ombar_e = v[8]
    cdef double ombar_f = v[9]
    cdef double ombar_g = v[10]
    cdef double ombar_h = v[11]
    cdef double ombar_i = v[12]
    cdef double ombar_de = v[13]
    cdef double z = v[14]
    cdef double dl = v[15]

    cdef double in_a = in_terms[0]
    cdef double in_b = in_terms[1]
    cdef double in_c = in_terms[2]
    cdef double in_d = in_terms[3]
    cdef double in_e = in_terms[4]
    cdef double in_f = in_terms[5]
    cdef double in_g = in_terms[6]
    cdef double in_h = in_terms[7]
    cdef double in_i = in_terms[8]
    cdef double in_j = in_terms[9]


    cdef double Hz = H0 * (ombar_r + ombar_m + ombar_a + ombar_b + ombar_c
                           + ombar_d + ombar_e + ombar_f + ombar_g + ombar_h
                           + ombar_i + ombar_de)**(0.5)

    if math.isnan(Hz):
        print('rainbow')
        print('z = %s, Hz = %s, in_terms = %s'% (z, Hz, in_terms))
        print('ombar_m = %s, ombar_r = %s, ombar_de = %s'
              % (ombar_m, ombar_r, ombar_de))

    cdef double dtdz = -1.0/((1.0+z) * Hz)
    cdef double dadz = -(1.0+z)**(-2.0)
    cdef double domdz = (3.0*ombar_m +in_d*ombar_m/Hz)/(1.0+z) # w = 0
    cdef double dordz = (4.0*ombar_r -in_d*ombar_m/Hz +in_c*ombar_r/Hz)/(1.0+z) # w = 1/3
    cdef double doadz = (2.7*ombar_a +in_a*ombar_a/Hz -in_h*ombar_h/Hz)/(1.0+z) # w = -0.1
    cdef double dobdz = (2.4*ombar_b -in_e*ombar_e/Hz +in_i*ombar_m/Hz)/(1.0+z) # w = -0.2
    cdef double docdz = (2.1*ombar_c +in_b*ombar_b/Hz -in_g*ombar_g/Hz)/(1.0+z) -in_f * ombar_f /(1.0+z) /Hz # w = -0.3
    cdef double doddz = (1.8*ombar_d +in_g*ombar_g/Hz +in_h*ombar_h/Hz)/(1.0+z) # w = -0.4
    cdef double doedz = (1.5*ombar_e -in_a*ombar_a/Hz -in_i*ombar_m/Hz)/(1.0+z)# w = -0.5
    cdef double dofdz = (1.2*ombar_f -in_b*ombar_b/Hz +in_e*ombar_e/Hz)/(1.0+z) # w = -0.6
    cdef double dogdz = (0.9*ombar_g +in_f*ombar_f/Hz)/(1.0+z) # w = -0.7
    cdef double dohdz = (0.6*ombar_h +in_h*ombar_h/Hz +in_j*ombar_de/Hz)/(1.0+z) # w = -0.8
    cdef double doidz = (0.3*ombar_i -in_h*ombar_h/Hz -in_j*ombar_de/Hz)/(1.0+z) # w = -0.9
    cdef double dodedz = -in_c*ombar_r/(1.0+z)/Hz # w = -1
    cdef double ddldz = 1.0/Hz

    # first derivatives of functions I want to find:
    f = [dtdz,# dt/dz (= f.d wrt z of time)
         dadz,# d(a)/dz (= f.d wrt z of scale factor)
         domdz,# dw = 0, (ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         dordz,# w = 1/3, d(ombar_r)/dz   (= f.d wrt z of density_r(t) / crit density(t0))
         doadz,# w = -0.1
         dobdz,# w = -0.2
         docdz,# w = -0.3
         doddz,# w = -0.4
         doedz,# w = -0.5
         dofdz,# w = -0.6
         dogdz,# w = -0.7
         dohdz,# w = -0.8
         doidz,# w = -0.9
         dodedz,# w = -1, d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit density(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)

    return f

def exotic(double[:] v, redshifts, in_terms, double H0):
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
    cdef double gamma = in_terms[0]
    cdef double zeta = in_terms[1]

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
         dodedz,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit density(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)

    return f

def late_intxde(double[:] v, redshifts, in_terms, double H0):
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
    cdef double gamma = in_terms[0]

    cdef double Hz = H0 * (ombar_m + ombar_de)**(0.5)

    if math.isnan(Hz):
        print('late_intxde')
        print('z = %s, Hz = %s, gamma = %s, ombar_m = %s, ombar_de = %s'
              %(z, Hz, gamma, ombar_m, ombar_de))

    irate = 0.0

    if z < 0.9:
        irate = ombar_de * ombar_m * gamma/(1.0+z)/Hz

    cdef double dtdz = -1.0/((1.0+z) * Hz)
    cdef double dadz = -(1.0+z)**(-2.0)
    cdef double domdz = 3.0*ombar_m /(1.0+z) - irate
    cdef double ddldz = 1.0/Hz

    # first derivatives of functions I want to find:
    f = [dtdz,# dt/dz (= f.d wrt z of time)
         dadz,# d(a)/dz (= f.d wrt z of scale factor)
         domdz,# d(ombar_m)/dz   (= f.d wrt z of density_m(t) / crit density(t0))
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit density(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)

    return f

def heaviside_late_int(double[:] v, redshifts, in_terms, double H0):
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
    cdef double gamma = in_terms[0]

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
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit density(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)

    return f

def late_int(double[:] v, redshifts, in_terms, double H0):
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
    cdef double gamma = in_terms[0]

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
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit density(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)

    return f

def expgamma(double[:] v, redshifts, in_terms, double H0):
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
    cdef double gamma = in_terms[0]

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
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit density(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)

    return f

def txgamma(double[:] v, redshifts, in_terms, double H0):
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
    cdef double gamma = in_terms[0]

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
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit density(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)

    return f

def zxgamma(double[:] v, redshifts, in_terms, double H0):
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
    cdef double gamma = in_terms[0]

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
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit density(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)

    return f

def gamma_over_z(double[:] v, redshifts, in_terms, double H0):
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
    cdef double gamma = in_terms[0]

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
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit density(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)

    return f

def zxxgamma(double[:] v, redshifts, in_terms, double H0):
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
    cdef double gamma = in_terms[0]

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
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit density(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)

    return f

def gammaxxz(double[:] v, redshifts, in_terms, double H0):
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
    cdef double gamma = in_terms[0]

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
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit density(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)

    return f

def rdecay_m(double[:] v, redshifts, in_terms, double H0):
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
    cdef double gamma = in_terms[0]

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
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit density(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)

    return f

def rdecay_de(double[:] v, redshifts, in_terms, double H0):
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
    cdef double gamma = in_terms[0]

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
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit density(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)

    return f

def rdecay_mxde(double[:] v, redshifts, in_terms, double H0):
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
    cdef double gamma = in_terms[0]

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
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit density(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)

    return f

def rdecay(double[:] v, redshifts, in_terms, double H0):
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
    cdef double gamma = in_terms[0]

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
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit density(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)

    return f

def interacting(double[:] v, redshifts, in_terms, double H0):
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
    cdef double gamma = in_terms[0]

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
         irate,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit density(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)

    return f

def LCDM(double[:] v, redshifts, in_terms, double H0):
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
         0.0,# d(ombar_de)/dz (= f.d wrt z of density_de(t) / crit density(t0))
         1.0,# d(z)/dz (= f.d wrt z of redshift)
         ddldz]# d(dl)/dz (= f.d wrt z of luminosty distance) # H + Hdz*(1+z)

    return f