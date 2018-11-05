3#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:47:44 2018

@author: BallBlueMeercat
"""
import random
import numpy as np
import zodesolve
import tools
import matplotlib.pyplot as plt

def redshift_picks(zmin, zmax, n):
    """
    Takes in:
        zmin = integer lowest redshift;
        zmax = integer highest redshift;
        n = integer number of redshifts to be generated.
    Returns:
        zpicks = list of randomly selected redshifts between zmin and zmax.
    """
    zinterval = (zmax - zmin) / (n*2)
    z_opts = tools.flist(zmin, zmax, zinterval)
    zpicks = random.sample(z_opts, n)
    zpicks = sorted(zpicks)
    return zpicks


def gnoise(mag, mu, sigma):
    """
   Returns:
       mag = mag, each point offset by unique Gaussian noise;
       noise = Gaussian noise.
    """
    n = len(mag)
    noise = np.random.normal(mu,sigma,n)
    mag = mag + noise

#    import matplotlib.pyplot as pl
#    from pylab import figure
#    figure()
#    pl.title('Noise distribution')
#    pl.hist(noise, 100)
#    pl.show()
    return mag


def magn(params, data, firstderivs_key, plot_key=False):
    """
    Finding matter density m, corrected absolute mag M, interaction gamma.

    Takes in:
        params = list of dictionaries {string:value} of names and
        starting values of parameters to be emcee fitted:
            [{'matter':int/float} = e_m(t)/ec(t0) at t=t0;
            {'Mcorr':int/float} = corrected absolute mag M;
            {'gamma':int/float} = interaction term;
            {'zeta':int/float}] = interaction term;
            ... (more)
        data = dictionary w/
            'colour': numpy.ndarray = SN colour;
            'x1': numpy.ndarray = SN stretch correction as;
            'zpicks':list of redshifts sorted in accending order;
            'mag':list of apparent magnitudes;

        firstderivs_key = string, indicates which firstderivs to integrate;
        plot_key = Boolean, to plot or not to plot model figures;
    Returns:
        mag = np.ndarray of apparent mag corresponding to input redshits.
    """
    zpicks = data['zpicks']

    # Corrected absolute magnitude M of SN.
    M = params[1]['Mcorr']

    dlpc, plot_var = zodesolve.zodesolve(params, zpicks, firstderivs_key)

    # Calculating apparent magnitudes of supernovae at the simulated
    # luminosity distances using the distance modulus formula.
    mag = 5 * np.log10(dlpc/10) + M

    if plot_key:
        # Plotting evolution of parameters in the model.
        import plots
        plots.modelcheck(mag, zpicks, plot_var, firstderivs_key)

    return mag


def noisy_mag(mu, sigma, params, data, firstderivs_key):

    model = magn(params, data, firstderivs_key)
    model = np.asarray(model)
    mag = gnoise(model, mu, sigma)

    return mag


def model_comparison(params, zpicks, firstderivs_key, label):
    """
    Takes in:
            params = list of 3 lists of dictionaries with model parameters;
            zpicks = list of redshifts to integrate over, in accending order;
            firstderivs_key = list of 3 strings, firstderivs for each params;
    Action:
        Plots parameter evolution for different params/models specified.

    """
    import plots

    plot_var_list = []

    for i in range(3):
        dlpc, plot_var = zodesolve.zodesolve(params[i], zpicks, firstderivs_key[i])

        # Corrected absolute magnitude M of SN.
        M = params[i][1]['Mcorr']

        # Apparent mags of SN at the luminosity
        # distances using the distance modulus formula.
        mag = 5 * np.log10(dlpc/10) + M

        plot_var['mag'] = mag
        plot_var_list.append(plot_var)

    plots.multi_model_plot(zpicks, firstderivs_key, plot_var_list, label)
    return


def makensavemagnz(params, data, firstderivs_key, mu, sigma, filename):
    '''
    Takes in:

        Parameters used to simulate magnitude:
    m_true                  = e_m(t)/e_crit(t0) at t=t0;
    de_true = 1 - m_true    = e_de(t)/e_crit(t0) at t=t0;
    g_true = interaction term, rate at which DE decays into matter.

        Statistical parameteres of gaussian noise added to data:
    mu =  mean;
    sigma = standard deviation;

    npoints = how many mag and z to generate.

        Model type:
    data_key = string, key for dictionary of interaction modes in firstderivs
    Options: 'Hdecay', 'rdecay', 'rdecay_de', 'rdecay_m', 'interacting', 'LCDM'
    Length of parameters has to correspond to the model being tested.

    filename = string, name of file data is saved to.

    Returns:
        Nothing. Generates redshifts and corresponding magnitudes (according
        to the model specified by data_key) offset by Gaussian noise,
        saves them into a binary file called filename in the working directory.
    '''
    zpicks = data['zpicks']

    mag = noisy_mag(mu, sigma, params, data, firstderivs_key)

    output = mag, zpicks

    # Relative path of output folder.
    save_path = './data/'+filename

    import pickle
    pickle.dump(output, open(save_path, 'wb'))

    return

def magn_plot(params, data, firstderivs_key, plot_key=False):
    """
    Finding matter density m, corrected absolute mag M, interaction gamma.

    Takes in:
        params = list of dictionaries {string:value} of names and
        starting values of parameters to be emcee fitted:
            [{'matter':int/float} = e_m(t)/ec(t0) at t=t0;
            {'Mcorr':int/float} = corrected absolute mag M;
            {'gamma':int/float} = interaction term;
            {'zeta':int/float}] = interaction term;
            ... (more)
        data = dictionary w/
            'colour': numpy.ndarray = SN colour;
            'x1': numpy.ndarray = SN stretch correction as;
            'zpicks':list of redshifts sorted in accending order;
            'mag':list of apparent magnitudes;

        firstderivs_key = string, indicates which firstderivs to integrate;
        plot_key = Boolean, to plot or not to plot model figures;
    Returns:
        mag = np.ndarray of apparent mag corresponding to input redshits.
    """
    zpicks = data['zpicks']

    # Corrected absolute magnitude M of SN.
    M = params[1]['Mcorr']

    dlpc, plot_var = zodesolve.zodesolve(params, zpicks, firstderivs_key)

    # Calculating apparent magnitudes of supernovae at the simulated
    # luminosity distances using the distance modulus formula.
    mag = 5 * np.log10(dlpc/10) + M
    z = plot_var['z']
    if firstderivs_key == 'waterfall':
        plt.figure()
        plt.title(r'$\bar \Omega$ evolution in waterfall')
        plt.xlabel('redshift')
        plt.ylabel(r'$\bar \Omega$')
        plt.plot(z, plot_var['ombar_m'], label='ombar_m vs redshift')
        plt.plot(z, plot_var['ombar_r'], label='ombar_r vs redshift')
        plt.plot(z, plot_var['a_ombar'], label='a_ombar vs redshift')
        plt.plot(z, plot_var['b_ombar'], label='b_ombar vs redshift')
        plt.plot(z, plot_var['c_ombar'], label='c_ombar vs redshift')
        plt.plot(z, plot_var['ombar_de'], label='ombar_de vs redshift')
        plt.legend()

        sum_om = plot_var['ombar_m'] + plot_var['ombar_r'] + plot_var['a_ombar']+ plot_var['b_ombar'] + plot_var['c_ombar'] + plot_var['c_ombar'] +plot_var['ombar_de']
        om_m = plot_var['ombar_m']/sum_om
        om_r = plot_var['ombar_r']/sum_om
        om_a = plot_var['a_ombar']/sum_om
        om_b = plot_var['b_ombar']/sum_om
        om_c = plot_var['c_ombar']/sum_om
        om_de = plot_var['ombar_de']/sum_om

        plt.figure()
        plt.title(r'$\Omega$ evolution in waterfall')
        plt.xlabel('redshift')
        plt.ylabel(r'$\Omega$')
        plt.plot(z, om_m, label = 'om_m')
        plt.plot(z, om_r, label = 'om_r')
        plt.plot(z, om_a, label = 'om_a')
        plt.plot(z, om_b, label = 'om_b')
        plt.plot(z, om_c, label = 'om_c')
        plt.plot(z, om_de, label = 'om_de')
        plt.legend()
        plt.show()

    elif firstderivs_key == 'LCDM':
        plt.figure()
        plt.title(r'$\bar \Omega$ evolution in LCDM')
        plt.xlabel('redshift')
        plt.ylabel(r'$\bar \Omega$')
        plt.plot(z, plot_var['ombar_m'], label='ombar_m vs redshift')
        plt.plot(z, plot_var['ombar_de'], label='ombar_de vs redshift')
        plt.legend()

        sum_om = plot_var['ombar_m'] + plot_var['ombar_de']
        om_m = plot_var['ombar_m']/sum_om
        om_de = plot_var['ombar_de']/sum_om

        plt.figure()
        plt.title(r'$\Omega$ evolution in LCDM')
        plt.xlabel('redshift')
        plt.ylabel(r'$\Omega$')
        plt.plot(z, om_m, label = 'om_m')
        plt.plot(z, om_de, label = 'om_de')
        plt.legend()
        plt.show()

    if plot_key:
        # Plotting evolution of parameters in the model.
        import plots
        plots.modelcheck(mag, zpicks, plot_var, firstderivs_key)

    return mag