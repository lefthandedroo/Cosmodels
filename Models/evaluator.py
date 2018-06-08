#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 21:45:40 2018

@author: BallBlueMeercat
"""
from pylab import figure, xlabel, ylabel, title, scatter, show, savefig
import time
import os.path
import numpy as np
<<<<<<< HEAD

from results import save
import paramfinder
from tools import runcount, timer
from datasim import data
=======
from scipy import integrate
from ln import lnprob
from results import save
import paramfinder
from tools import runcount, timer, path
from datasim import data


>>>>>>> 43dbdc3f8d07c0144cad572808e36f870b5e87aa

# Parameters used to simulate data:  
m_true = 0.3           # (= e_m(t)/e_crit(t0) at t=t0).
de_true = 1 - m_true   # (de = e_de(t)/e_crit(t0) at t=t0).
g_true = 0             # Interaction term, rate at which DE decays into matter.

# Number of datapoints to be simulated and number of emcee steps.
npoints, nsteps = 1000, 1000

# Statistical parameteres of noise:
mu = 0            # mean
sigma = 0.01      # standard deviation

# Key for the dictionary of interaction modes in firstderivs
# 'Hdecay', 'rdecay', 'rdecay_de', 'rdecay_m', 'interacting', 'LCDM':LCDM
<<<<<<< HEAD
data_firstderivs_key = 'LCDM'
data_params = {'m':m_true}
test_firstderivs_key = data_firstderivs_key#'rdecay'

=======
>>>>>>> 43dbdc3f8d07c0144cad572808e36f870b5e87aa
# Length of parameters has to correspond to the model being tested.
data_key = 'LCDM'
if data_key == 'LCDM':
    data_params = {'m':m_true}
else:
    data_params = {'m':m_true, 'gamma':g_true}

test_key = 'rdecay'




def repeatrun():
    # Changing directory to dedicated folder for saving output.
    save_path, directory = path()

    if test_key == 'LCDM':
        params = {'m':m_true}
    else:
        params = {'m':m_true, 'gamma':g_true}
    i = 0
    while i < 1:
        print('_____________________ run number',i)
        
<<<<<<< HEAD
        mag = data(mu, sigma, npoints, data_params, data_firstderivs_key)
        
        propert, sampler = paramfinder.paramfinder(
                npoints, nsteps, sigma, mu, params, 
                save_path, mag, test_firstderivs_key)
=======
        mag, zpicks = data(mu, sigma, npoints, data_params, data_key)
        
        propert, sampler = paramfinder.paramfinder(
                npoints, nsteps, sigma, mu, params, zpicks, 
                mag, test_key, save_path)
>>>>>>> 43dbdc3f8d07c0144cad572808e36f870b5e87aa
        i += 1
    
    # Saving sampler to directory.
    save(save_path, 'sampler', sampler)
<<<<<<< HEAD
    print('Type of interaction in the model being fitted to data:', 
          test_firstderivs_key)
    print('Data is simulated using',data_firstderivs_key)
    print()
=======

    print('Type of interaction in the model being tested:', 
          test_key)
    print('Data is simulated using',data_key)
>>>>>>> 43dbdc3f8d07c0144cad572808e36f870b5e87aa
    print('directory:',directory)

    return

<<<<<<< HEAD
repeatrun()

def modeltest(mag, test_firstderivs_key, save_path):
    
    propert, sampler = paramfinder.paramfinder(
            npoints, nsteps, sigma, mu, params,
            save_path, mag, test_firstderivs_key)
    
    trace = propert.get('trace', 'no_trace')
    if trace == 'no_trace':
        print('modeltest found no trace in propert from paramfinder')
    
    return trace

def Bfactor():
    
    import zpicks # DO NOT MOVE THIS LINE UP
    zpicks = zpicks.zpicks(0.005, 2, npoints)
    
    import datasim
    model = datasim.mag(params, zpicks, 'LCDM')
    model = np.asarray(model)
    mag = datasim.gnoise(model, mu, sigma)

    trace_1D = modeltest(mag, 'LCDM')
    trace_2D = modeltest(mag, 'rdecay')
    
    return
=======
#repeatrun()


def integrate_posterior_1D(lnprob, xlim, theta, zpicks, 
                           mag, sigma, test_key):
    func = lambda theta: np.exp(lnprob(theta, zpicks, mag, sigma, 
                                       test_key))
    return integrate.dblquad(func, xlim[0], xlim[1])


def integrate_posterior_2D(lnprob, xlim, ylim, theta, zpicks, 
                           mag, sigma, test_key):
    func = lambda theta1, theta0: np.exp(lnprob([theta0, theta1], 
                                                zpicks, mag, sigma, test_key))
    return integrate.dblquad(func, xlim[0], xlim[1],
                             lambda x: ylim[0], lambda x: ylim[1])


def modeltest(npoints, nsteps, sigma, mu, 
              zpicks, mag, test_key, save_path):
    
    print('__________ testing', test_key)
    
    if test_key == 'LCDM':
        params = {'m':m_true}
    else:
        params = {'m':m_true, 'gamma':g_true}
    
    propert, sampler = paramfinder.paramfinder(
            npoints, nsteps, sigma, mu, params, zpicks, 
            mag, test_key, save_path)
    
    trace = propert.get('trace',)    
#    if not trace:
#        print ('modeltest got no trace for %s from paramfinder'%(test_key))
    
    theta = []
    for key in params:
        theta.append(propert.get(key,))
        
    return trace, theta


def Bfactor(npoints, nsteps, sigma, mu, data_params, data_key, M1_key, M2_key):
    # Changing directory to dedicated folder for saving output.
    save_path, directory = path()
    
    print('Generating magnitudes...')
    mag, zpicks = data(mu, sigma, npoints, data_params, data_key)

    
    trace_1D, theta_1D = modeltest(npoints, nsteps, sigma, mu,
                         zpicks, mag, M1_key, save_path)
    print('theta_1D:',theta_1D)
    trace_2D, theta_2D = modeltest(npoints, nsteps, sigma, mu,
                         zpicks, mag, M2_key, save_path)
    print('theta_2D:',theta_2D)
    
#    xlim = zip(trace_1D.min(0), trace_1D.max(0))
#    Z1, err_Z1 = integrate_posterior_1D(lnprob, xlim, theta_1D, 
#                                        zpicks, mag, sigma, M1_key)
#    print("Z1 =", Z1, "+/-", err_Z1)
    
    xlim, ylim = zip(trace_2D.min(0), trace_2D.max(0))
    Z2, err_Z2 = integrate_posterior_2D(lnprob, xlim, ylim, theta_2D, 
                                        zpicks, mag, sigma, M2_key)
    print("Z2 =", Z2, "+/-", err_Z2)
                  
#    print("Bayes factor:", Z2 / Z1)
    
    print('Data is simulated using',data_key)
    print()
    print('directory:',directory)
    
    return trace_1D, trace_2D

Bfactor(npoints, nsteps, sigma, mu, data_params, data_key, data_key, test_key)
>>>>>>> 43dbdc3f8d07c0144cad572808e36f870b5e87aa


def errorvsdatasize(params):
    # Script timer.
    timet0 = time.time()
    
    sigma = 0.02
    sigma_max = 0.03
    sigma_step = 0.05
    npoints_min = 1000
    npoints_max = 1100
    npoints_step = 3000
    
    # How many iterations have I signed up for?
    runcount(sigma, sigma_max, sigma_step,
              npoints_min, npoints_max, npoints_step)
    
    decision = input('Happy with the number of iterations? (enter=yes) ')
    if not decision:
        pass
    else:
        return
    
    # Folder for saving output.
    directory = str(int(time.time()))
    # Relative path of output folder.
    save_path = './results/'+directory
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
        mag, zpicks = data(mu, sigma, npoints, data_params, data_key)
        
        while npoints < npoints_max:
            print('_____________________ run number',run)
            propert, sampler = paramfinder.paramfinder(       
            npoints, nsteps, sigma, mu, params, zpicks, 
            mag, test_key, save_path)
            
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
    figure()
    xlabel('size of dataset')
    ylabel('standard deviation of marginalised m distribution')
    title('sd of m vs size of dataset, sd of noise = %s'%(sigma))
    scatter(npoints_l, m_sd_l, c='m')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_sd_of_m_.png'
    filename = os.path.join(save_path, filename)
    savefig(filename)
    show()
    
    figure()
    xlabel('size of dataset')
    ylabel('mean of marginalised m distribution')
    title('mean of m vs size of dataset, sd of noise = %s'%(sigma))
    scatter(npoints_l, m_mean_l, c='c')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_mean_of_m_.png'
    filename = os.path.join(save_path, filename)
    savefig(filename)
    show()
    
    figure()
    xlabel('size of dataset')
    ylabel('variance coefficient in %')
    title('sd/mean of m vs size of dataset, sd of noise = %s'%(sigma))
    scatter(npoints_l, m_vc_l, c='coral')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_cv_of_m_.png'
    filename = os.path.join(save_path, filename)
    savefig(filename)
    show()

    # gamma
    figure()
    xlabel('size of dataset')
    ylabel('standard deviation of marginalised gamma distribution')
    title('sd of gamma vs size of dataset, sd of noise = %s'%(sigma))
    scatter(npoints_l, g_sd_l, c='m')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_sd_of_g_.png'
    filename = os.path.join(save_path, filename)
    savefig(filename)
    show()
    
    figure()
    xlabel('size of dataset')
    ylabel('mean of marginalised gamma distribution')
    title('mean of gamma vs size of dataset, sd of noise = %s'%(sigma))
    scatter(npoints_l, g_mean_l, c='c')        
    stamp = str(int(time.time()))
    filename = str(stamp)+'_mean_of_g_.png'
    filename = os.path.join(save_path, filename)
    savefig(filename)
    show()
    
#    figure()
#    xlabel('size of dataset')
#    ylabel('variance coefficient in %')
#    title('sd/mean of gamma vs size of dataset, sd of noise = %s'%(sigma))
#    scatter(npoints_l, g_vc_l, c='coral')        
#    stamp = str(int(time.time()))
#    filename = str(stamp)+'_cv_of_g_.png'
#    filename = os.path.join(save_path, filename)
#    savefig(filename)
#    show()
        
    # Saving results to directory.
    save(save_path, 'm_vc', m_vc_l)
    save(save_path, 'm_sd', m_sd_l)
    save(save_path, 'm_mean', m_mean_l)

    save(save_path, 'g_vc', g_vc_l)
    save(save_path, 'g_sd', g_sd_l)
    save(save_path, 'g_mean', g_mean_l)
    
    save(save_path, 'sigma', sigma_l)
    save(save_path, 'npoints', npoints_l)
    save(save_path, 'sampler', sampler_l)
    
    print('directory:',directory)
    
    # Time taken by evaluator. 
    timet1=time.time()
    timer('evaluator', timet0, timet1)
    
    return

#errorvsdatasize()