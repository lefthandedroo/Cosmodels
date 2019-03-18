#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:52:36 2018

@author: BallBlueMeercat
"""

import numpy as np
from datasim import magn

def lnlike(theta, data, sigma, firstderivs_key, names):
    '''
    Takes in:
        theta = numpy.ndarray, guessed values for parameters;
        data = dictionary, {'mag':mag, 'zpicks':zpicks};
        sigma = int/float, error on data;
        firstderivs_key = string, model being tested;
        params = list of dictionaries {string:value} of names and
        current guessed values of parameters being emcee fitted:
            [{'matter':float} = e_m(t)/ec(t0) at t=t0;
            {'Mcorr':float} = corrected absolute mag M;
            {'gamma':float} = interaction term;
            {'zeta':float} = interaction term;
            ... {'parameter':value})].
    Returns:
        float, likelihood for firstderivs_key model with theta parameters.
    '''
    mag = data['mag']
    model_mag, model_da = magn(names, theta, data, firstderivs_key) # mag, but with theta params
    var = sigma**2
    likelihood = -0.5*np.sum((mag-model_mag)**2.0 /var +np.log(2.0*np.pi*var))
    return likelihood

def lnprior(th, key):
    '''
    Takes in:
        th = numpy.ndarray, guessed values for parameters;

    Returns:
        0.0 if all conditions on theta values are met;
        -np.inf if theta values are outside of prior.
    '''
    Mcorr_min, Mcorr_max = -20, -18
    assert len(th) > 0, f'len(th) = 0, theta = {th}'

    # Corrected magnitude of SN Ia and Omega_bar matter
    if Mcorr_min < th[0] < Mcorr_max and (0 < th[1] < 1 or th[1] == 1):

        if key == 'rainbow':
            # ombar_radiation, a_ombar, b_ombar, c_ombar
            if 0 < th[2] < 1 and 0 < th[3] < 1 and 0 < th[4] < 1 and 0 < th[5] < 1:
                # , d_ombar, e_ombar, f_ombar, g_ombar
                if 0 < th[6] < 1 and 0 < th[7] < 1 and 0 < th[8] < 1 and 0 < th[9] < 1:
                    # h_ombar, i_ombar
                    if 0 < th[10] < 1 and 0 < th[11] < 1:
                        l = 0.01
                        # a_in, b_in, c_in
                        if -l < th[12] < l and -l < th[13] < l and -l < th[14] < l:
                            # d_in, e_in, f_in
                            if  -l < th[15] < l and -l < th[16] < l and -l < th[17] < l:
                                # g_in, h_in, i_in
                                if -l < th[18] < l and -l < th[19] < l and -l < th[20] < l:
                                     #j_in, k_in
                                     if -l < th[21] < l and -l < th[22] < l:
                                         return 0.0
        elif key == 'kanangra':
            # ombar_radiation, a_ombar, b_ombar, c_ombar
            if 0 < th[2] < 1 and 0 < th[3] < 1 and 0 < th[4] < 1 and 0 < th[5] < 1:
                # d_ombar, e_ombar
                if 0 < th[6] < 1 and 0 < th[7] < 1:
                    l = 1
                    # a_in, b_in, c_in, d_in, e_in, f_in, g_in
                    if -l<th[8]<l and -l<th[9]<l and -l<th[10]<l:
                        if -l<th[11]<l and -l<th[12]<l and -l<th[13]<l and -l<th[14]<l:
                            return 0.0

        elif key == 'waterfall':
            # ombar_radiation, a_ombar, b_ombar, c_ombar
            if 0 < th[2] < 1 and 0 < th[3] < 1 and 0 < th[4] < 1 and 0 < th[5] < 1:
                l = 1
                # v_in, w_in, x_in, y_in, z_in
                if -0.01 < th[6] < l and -0.01 < th[7] < l and -0.01 < th[8] < l:
                    if -0.01 < th[9] < l and -0.01 < th[10] < l:
#                    if abs(th[6]) < l and abs(th[7]) < l and -0.5< th[8] < l:
#                        if abs(th[9]) < l and abs(th[10]) < l:
                        return 0.0

        elif key == 'stepfall':
            assert len(th) > 3, f'len(thet) = {len(th)}, theta = {th}'
            # ombar_radiation, a_ombar
            if 0 < th[2] < 1 and 0 < th[3] < 1:
                # v_in, w_in, x_in
                l = 1
                if -0.01 < th[4] < l and -0.01 <th[5] < l and -0.01 < th[6] < l:
                    return 0.0

        elif key == 'exotic':
            # ombar_radiation
            if 0 < th[2] < 1:
                # gamma, zeta
#                if -100 < th[3] < 100 and -100 < th[4] < 100:
                if -0.01 < th[3] < 1 and -0.01 < th[4] < 1:

                    return 0.0

        elif key == 'late_intxde':
            if -2 < th[2] < 0.1: # interaction term
                return 0.0

        elif key == 'heaviside_late_int':
            if -1.45 < th[2] < 0.1: # interaction term
                return 0.0

        elif key == 'late_int':
            if -15 < th[2] < 0.1: # interaction term
                return 0.0

        elif key == 'expgamma':
            if -0.1 < th[2] < 1.5: # interaction term
                return 0.0

        elif key == 'txgamma':
            if -0.5 < th[2] < 0.1: # interaction term
                return 0.0

        elif key == 'zxgamma':
            if -10 < th[2] < 0.1: # interaction term
                return 0.0

        elif key == 'zxxgamma':
            if -0.1 < th[2] < 12: # interaction term
                return 0.0

        elif key == 'gammaxxz':
            if -1 < th[2] < 1: # interaction term
                return 0.0

        elif key == 'rdecay':
#            if -2 < th[2] < 0.1: # interaction term
            if -10 < th[2] < 1: # interaction term
                return 0.0

        elif key == 'interacting':
            if -1.5 < th[2] < 0.1: # interaction term
                return 0.0

        elif key == 'LCDM':
            return 0.0

        else:
            if abs(th[2]) < 10: # interaction term
                return 0.0

    return -np.inf

def lnprob(theta, data, sigma, firstderivs_key, names):
    lp = lnprior(theta, firstderivs_key)
#    print(theta)
#    str_theta = str(theta[1])
#    last_digits = int(str_theta[-3:])
#    if last_digits == 11:
#        print(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, data, sigma, firstderivs_key, names)