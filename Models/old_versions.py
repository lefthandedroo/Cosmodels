#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 19:05:26 2018

@author: BallBlueMeercat
"""
import numpy as np 
import matplotlib.pyplot as plt
import time
import os
import zodesolve
import tools
import datasim
import results
#datasim.py

def magn(params, data, firstderivs_key, plot_key=False):
    """
    Finding matter density m, interaction gamma.
    
    Takes in:
            params = dictionary with true parameters;
            zpicks = list of redshifts to integrate over, in accending order;
            firstderivs_key = string, indicates which firstderivs to integrate;
            plot_key = Boolean, to plot or not to plot model figures;
    Returns:
        mag = np.ndarray of apparent mag corresponding to input redshits.
    """
#    print('@@@ magn has been called')
    if firstderivs_key == 'LCDM':
        params['gamma'] = 0
        del params['gamma']
    
    zpicks = data['zpicks']
    
    # Absolute brightness of supernovae.
    M = -19
    
    dlpc, plot_var = zodesolve.zodesolve(params, zpicks, firstderivs_key)
    
    # Calculating apparent magnitudes of supernovae at the simulated
    # luminosity distances using the distance modulus formula.
    mag = 5 * np.log10(dlpc/10) + M

    if plot_key:
        # Checking evolution of the model.
        import plots
        plots.modelcheck(mag, zpicks, plot_var, firstderivs_key)
        
    return mag

#def magn(params, data, firstderivs_key, plot_key=False):
#    """
#    Finding matter density m, alpha, beta, interaction gamma.
#    Takes in:
#            params = dictionary with true parameters;
#            zpicks = list of redshifts to integrate over, in accending order;
#            firstderivs_key = string, indicates which firstderivs to integrate;
#            plot_key = Boolean, to plot or not to plot model figures;
#    Returns:
#        mag = np.ndarray of apparent mag corresponding to input redshits.
#    """
##    print('@@@ magn has been called')
#    if firstderivs_key == 'LCDM':
#        params['gamma'] = 0
#        del params['gamma']
#    
#    zpicks = data['zpicks']
#    x1 = data['x1']
#    colour = data['colour']
#    
#    # Absolute brightness of supernovae.
#    M = params['M']
#    alpha = params['alpha']
#    beta = params['beta']
#    
#    dlpc, plot_var = zodesolve.zodesolve(params, zpicks, firstderivs_key)
#    
#    # Calculating apparent magnitudes of supernovae at the simulated
#    # luminosity distances using the distance modulus formula.
#    mag = 5 * np.log10(dlpc/10) + M - alpha*x1 +beta*colour
#
#    if plot_key:
#        # Checking evolution of the model.
#        import plots
#        plots.modelcheck(mag, zpicks, plot_var, firstderivs_key)
#        
#    return mag
    
# Slow mag calculation
#    # Calculating apparent magnitudes of supernovae at the simulated
#    # luminosity distances using the distance modulus formula.
#    mag = []   
#    for i in range(len(dlpc)):
#        if dlpc[i] == 0:
#            magnitude = M
#        else:
#            # magnitude from the distance modulus formula
#            magnitude = 5 * math.log10(dlpc[i]/10) + M
#        mag.append(magnitude)
    
#evaluator.py
    
dataname = 'mag_z_LCDM_1000_sigma_'+str(sigma)

# Number of datapoints to be simulated
npoints = 1000
zmax = 2

#    Making data (mag and z)
#dataname = 'mag_z_'+data_key+'_'+str(npoints)+'_sigma_'+str(sigma)
#datasim.makensavemagnz(m_true, g_true, mu, sigma, zpicks, data_key, dataname)

#    Making redshifts to use in this script.
#zpicks = datasim.redshift_picks(0.005, zmax, npoints)

def errorvsdatasize():
    
    data_key = 'LCDM'
    test_key = 'late_int'
    
    data_params = {'m':0.3, 'gamma':0}
    test_params = {'m':0.3, 'gamma':0}
    
    # Script timer.
    timet0 = time.time()
    
    sigma = 0.02
    sigma_max = 0.03
    sigma_step = 0.05
    npoints_min = 1000
    npoints_max = 1100
    npoints_step = 3000
    
    # How many iterations have I signed up for?
    tools.runcount(sigma, sigma_max, sigma_step,
              npoints_min, npoints_max, npoints_step)
    
    decision = input('Happy with the number of iterations? (enter=yes) ')
    if  decision:
        return
    
    # Folder for saving output.
    directory = str(int(time.time()))
    # Relative path of output folder.
    save_path = './emcee_results/'+directory
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    
    run = 0
    
    m_sd_l = []
    m_mean_l = []
    m_vc_l = []    

    g_sd_l = []
    g_mean_l = []
    g_vc_l = [] 
    
    sigma_l = []
    npoints_l = []
    sampler_l = []
    
    while sigma < sigma_max:

        npoints = npoints_min 
        
        # Data to be used:
        mag = datasim.noisy_mag(mu, sigma, npoints, data_params, data_key)
        
        while npoints < npoints_max:
            print('_____________________ run number',run)
          
            propert, sampler = stats.stats(test_params, zpicks, mag, 
                                           sigma, nsteps, save_path, test_key)
            
            m_sd = propert.get('m_sd',0)
            m_mean = propert.get('m_mean', 0)
            m_vc = m_sd/m_mean * 100
            m_vc_l.append(m_vc)
            m_sd_l.append(m_sd)
            m_mean_l.append(m_mean)
            
            g_sd = propert.get('gamma_sd', 0)
            g_mean = propert.get('gamma_mean', 0)
            g_vc = g_sd/g_mean * 100
            g_vc_l.append(g_vc)
            g_sd_l.append(g_sd)
            g_mean_l.append(g_mean)                        
            
            sigma_l.append(sigma)
            npoints_l.append(npoints)
            sampler_l.append(sampler)
            
            npoints += npoints_step
            run += 1
        
        sigma += sigma_step
        
    # Saving plots to run directory.
    # m
    plt.figure()
    plt.xlabel('size of dataset')
    plt.ylabel('standard deviation of marginalised m distribution')
    plt.title('sd of m vs size of dataset, sd of noise = %s'%(sigma))
    plt.scatter(npoints_l, m_sd_l, c='m')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_sd_of_m_.pdf'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    plt.show()
    
    plt.figure()
    plt.xlabel('size of dataset')
    plt.ylabel('mean of marginalised m distribution')
    plt.title('mean of m vs size of dataset, sd of noise = %s'%(sigma))
    plt.scatter(npoints_l, m_mean_l, c='c')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_mean_of_m_.pdf'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    plt.show()
    
    plt.figure()
    plt.xlabel('size of dataset')
    plt.ylabel('variance coefficient in %')
    plt.title('sd/mean of m vs size of dataset, sd of noise = %s'%(sigma))
    plt.scatter(npoints_l, m_vc_l, c='coral')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_cv_of_m_.pdf'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    plt.show()

    # gamma
    plt.figure()
    plt.xlabel('size of dataset')
    plt.ylabel('standard deviation of marginalised gamma distribution')
    plt.title('sd of gamma vs size of dataset, sd of noise = %s'%(sigma))
    plt.scatter(npoints_l, g_sd_l, c='m')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_sd_of_g_.pdf'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    plt.show()
    
    plt.figure()
    plt.xlabel('size of dataset')
    plt.ylabel('mean of marginalised gamma distribution')
    plt.title('mean of gamma vs size of dataset, sd of noise = %s'%(sigma))
    plt.scatter(npoints_l, g_mean_l, c='c')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_mean_of_g_.pdf'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    plt.show()
    
#    plt.figure()
#    plt.xlabel('size of dataset')
#    plt.ylabel('variance coefficient in %')
#    plt.title('sd/mean of gamma vs size of dataset, sd of noise = %s'%(sigma))
#    plt.scatter(npoints_l, g_vc_l, c='coral')        
#    plt.stamp = str(int(time.time()))
#    plt.filename = str(stamp)+'_cv_of_g_.pdf'
#    plt.filename = os.path.join(save_path, filename)
#    plt.savefig(filename)
#    plt.show()
        
    # Saving results to directory.
    results.save(save_path, 'm_vc', m_vc_l)
    results.save(save_path, 'm_sd', m_sd_l)
    results.save(save_path, 'm_mean', m_mean_l)

    results.save(save_path, 'g_vc', g_vc_l)
    results.save(save_path, 'g_sd', g_sd_l)
    results.save(save_path, 'g_mean', g_mean_l)
    
    results.save(save_path, 'sigma', sigma_l)
    results.save(save_path, 'npoints', npoints_l)
    results.save(save_path, 'sampler', sampler_l)
    
    print('directory:',directory)
    
    # Time taken by evaluator. 
    timet1=time.time()
    tools.timer('evaluator', timet0, timet1)
    
    return

#errorvsdatasize()
    

#ln.py
    
#def lnlike(theta, data, sigma, firstderivs_key, ndim):
#    '''
#    Finding matter density m, interaction gamma.
#    '''
#    mag = data['mag']
#    
#    params = {}
#    if ndim == 1:
#        params = {'m':theta}
#    elif ndim == 2:
#        params = {'m':theta[0],'gamma':theta[1]}
#    
#    model = magn(params, data, firstderivs_key)
#    var = sigma**2
#    return -0.5*np.sum((mag-model)**2 /var +0.5*np.log(2*np.pi*var))


#def lnlike(theta, data, sigma, firstderivs_key, ndim):
#    '''
#    Finding matter density m, corrected absolute mag M, interaction gamma.
#    '''    
#    mag = data['mag']
#    
#    params = {}
#    if ndim == 2:
#        params = {'m':theta[0], 'M':theta[1]}
#    elif ndim == 3:
#        params = {'m':theta[0],'M':theta[1], 'gamma':theta[2]}
#    
#    model = magn(params, data, firstderivs_key)
#    var = sigma**2
#    return -0.5*np.sum((mag-model)**2 /var +0.5*np.log(2*np.pi*var))
    
#def lnprior(theta, key):
#    '''
#    Finding matter density m, interaction gamma.
#    '''
#    
#    if key == 'LCDM':
#        m = theta
#        if 0 < m < 1 or m == 1:
#            return 0.0
#    elif key == 'late_int' or 'heaviside_late_int' or 'late_intxde':
#        m, gamma = theta
#        if (0 < m < 1 or m == 1) and -1.45 < gamma < 0.2:
#            return 0.0       
#    elif key == 'rdecay':
#        m, gamma = theta
#        if (0 < m < 1 or m == 1) and -10 < gamma < 0:
#            return 0.0
#    elif key == 'interacting':
#        m, gamma = theta
#        if (0 < m < 1 or m == 1) and abs(gamma) < 1.45:
#            return 0.0
#    elif key == 'expgamma':
#        m, gamma = theta
#        if (0 < m < 1 or m == 1) and abs(gamma) < 25:
#            return 0.0
#    elif key == 'zxxgamma' or 'gammaxxz':
#        m, gamma = theta
#        if (0 < m < 1 or m == 1) and 0 < gamma < 10:
#            return 0.0        
#    else:
#        m, gamma = theta
#        if (0 < m < 1 or m == 1) and abs(gamma) < 10:
#            return 0.0
#        
#    return -np.inf
#    
#
#def lnprior(theta, key):
#    '''
#    Finding matter density m, corrected absolute mag M, interaction gamma.
#    '''  
#    
#    Mmin = -20
#    
#    Mmax = -18
#    
#    if key == 'LCDM':
#        m, M = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax:
#            return 0.0
#    elif key == 'late_int' or 'heaviside_late_int' or 'late_intxde':
#        m, M, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and -1.45 < gamma < 0.2:
#            return 0.0       
#    elif key == 'rdecay':
#        m, M, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and -10 < gamma < 0 :
#            return 0.0
#    elif key == 'interacting':
#        m, M, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(gamma) < 1.45:
#            return 0.0
#    elif key == 'expgamma':
#        m, M, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(gamma) < 25 :
#            return 0.0
#    elif key == 'zxxgamma' or 'gammaxxz':
#        m, M, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and 0 < gamma < 10:
#            return 0.0        
#    else:
#        m, M, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(gamma) < 10:
#            return 0.0
#        
#    return -np.inf
#    
#
#def lnprior(theta, key):
#    '''
#    Finding matter density m, absolute M, alpha, beta, interaction gamma.
#    '''  
#    
#    Mmin, Mmax = -20, -18
#    amax = 5
#    bmax = 5
#    
#    print('key ln prior gets is = ',key)
#    
#    if key == 'LCDM':
#        m, M, alpha, beta = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(alpha) < amax and abs(beta) < bmax:
#            return 0.0
#    elif key == 'late_int' or key == 'heaviside_late_int' or key == 'late_intxde':
#        m, M, alpha, beta, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(alpha) < amax and abs(beta) < bmax and -1.45 < gamma < 0.2:
#            return 0.0       
#    elif key == 'rdecay':
#        m, M, alpha, beta, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(alpha) < amax and abs(beta) < bmax and -10 < gamma < 0 :
#            return 0.0
#    elif key == 'interacting':
#        m, M, alpha, beta, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(alpha) < amax and abs(beta) < bmax and abs(gamma) < 1.45:
#            return 0.0
#    elif key == 'expgamma':
#        m, M, alpha, beta, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(alpha) < amax and abs(beta) < bmax and abs(gamma) < 25 :
#            return 0.0
#    elif key == 'zxxgamma' or key == 'gammaxxz':
#        m, M, alpha, beta, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(alpha) < amax and abs(beta) < bmax and 0 < gamma < 10:
#            return 0.0
#    elif key == 'exotic':
#        m, M, alpha, beta, gamma, zeta = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(alpha) < amax and abs(beta) < bmax and 0 < gamma < 10 and 0 < zeta < 10:
#            return 0.0
#    else:
#        m, M, alpha, beta, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(alpha) < amax and abs(beta) < bmax and abs(gamma) < 10:
#            return 0.0
#        
#    return -np.inf
    
# zodesolve.py

#def odesolve(params, zpicks, firstderivs_key):
#    """
#    Takes in:
#        gamma = interaction constant;
#        m = e_m(t)/ec(t0) at t=t0;
##        de = e_de(t)/ec(t0) at t=t0.
#    Returns: 
#        z = numpoints number of redshifts zmin<z<zmax;
#        dlpc = luminosity distance in pc.
#    
#    """
##    print('@@ zodesolve has been called')
#
#    # Inserting 0 at the front of redshifts to allow initial conditions.
#    zpicks = [0.0] + zpicks
#    
#    # Standard cosmological parameters.
#    H0 = 1
#    c_over_H0 = 4167 * 10**6    # c/H0 in parsecs
#    
#    # Initial conditions at z = 0 (now).
#    dl0 = 0             # luminosity distance
#    rho_c0 = H0**2      # critical density
#    ombar_m0 = params.get('m', 0)                        # e_m(z)/ec(z=0)
#    gamma = params.get('gamma',0)
#    ombar_de0 = params.get('de', rho_c0/rho_c0 -ombar_m0)# e_de(z)/ec(z=0)
#    
#    # ODE solver parameters:
#    abserr = 1.0e-8
#    relerr = 1.0e-6
#    
#    # Pack up the initial conditions and eq of state parameters.
#    v0 = [ombar_m0, ombar_de0, dl0]
#    
#    # Extracting the parsed mode of interaction.
#    firstderivs_function = firstderivs_functions.get(firstderivs_key,0)
#    if firstderivs_function == 0:
#        print("firstderivs_functions dict didn't have the key zodeosolve asked for")
#    
#    # Call the ODE solver.
#    vsol = odeint(firstderivs_function, v0, zpicks, args=(gamma,H0), 
#                  atol=abserr, rtol=relerr)
#            
#    # Separate results into their own arrays:
#    z = np.asarray(zpicks)
#    z = z[1:]
#    ombar_m = vsol[1:,0]
#    ombar_de = vsol[1:,1]
#    dl = vsol[1:,2] * (1+z)  # in units of dl*(H0/c)
#    dlpc = dl * c_over_H0    # dl in parsecs (= vsol[dl] * c/H0)
#    
#    plot_var = dlpc, dl, ombar_m, gamma, ombar_de, ombar_m0, ombar_de0
#    
#    return dlpc, plot_var
