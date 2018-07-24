#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 21:45:40 2018

@author: BallBlueMeercat
"""
import matplotlib.pyplot as plt
import time
import os.path
import results
import paramfinder
import tools
import datasim

# Parameters used to simulate data:  
m_true = 0.3           # (= e_m(t)/e_crit(t0) at t=t0).
de_true = 1.0 - m_true   # (de = e_de(t)/e_crit(t0) at t=t0).
g_true = 0#-0.50           # Interaction, rate at which DE decays into matter.

# Number of datapoints to be simulated and number of emcee steps.
npoints, nsteps = 1000, 10000
zmax = 2

# Statistical parameteres of noise:
mu = 0.0            # mean
sigma = 0.07      # standard deviation

# Key for the dictionary of interaction modes in firstderivs
# 'expgamma','txgamma','zxgamma','gamma_over_z','zxxgamma','gammaxxz',
#'rdecay_m','rdecay_de','rdecay_mxde','rdecay','interacting','LCDM'
# Length of parameters has to correspond to the model being tested.

data_key = 'LCDM'

test_key = 'expgamma'

if data_key == 'LCDM':
    data_params = {'m':m_true}
else:
    data_params = {'m':m_true, 'gamma':g_true}

if test_key == 'LCDM':
    test_params = {'m':m_true}
else:
    test_params = {'m':m_true, 'gamma':g_true}

#dataname = 'mag_z_'+data_key+'_'+str(npoints)+'_sigma_'+str(sigma)
#datasim.makensavemagnz(m_true, g_true, mu, sigma, zpicks, data_key, dataname)
               
zpicks = datasim.redshift_picks(0.005, zmax, npoints)

def modelcheck():

    datasim.magn(test_params, zpicks, test_key, plot_key=True)

    return

#modelcheck()

#def modelcheck():
#    g = [-10, -5, 0, 5, 10]
#    datasim.magn(test_params, zpicks, test_key, plot_key=True, gamma_list = g)
#    return
#
#modelcheck()

def all_modelcheck():
    
    import firstderivs as f
    
    firstderivs_functions = {'expgamma':f.expgamma,
                         'txgamma':f.txgamma,
                         'zxgamma':f.zxgamma,
                         'gamma_over_z':f.gamma_over_z,
                         'zxxgamma':f.zxxgamma,
                         'gammaxxz':f.gammaxxz,
                         'rdecay_m':f.rdecay_m,
                         'rdecay_de':f.rdecay_de,
                         'rdecay_mxde':f.rdecay_mxde,
                         'rdecay':f.rdecay,                         
                         'interacting':f.interacting,
                         'LCDM':f.LCDM}
    
    for key in firstderivs_functions:
        test_key = key

    if key == 'LCDM':
        test_params = {'m':m_true}
    else:
        test_params = {'m':m_true, 'gamma':g_true}
    
        datasim.noisy_mag(zpicks, mu, sigma, test_params, 
                               test_key, plot_key=True)
    return

#all_modelcheck()


def quickemcee():

    
    i = 0
    while i < 1:
        if i > 0:
            print('_____________________ run number',i)
        
#        mag = datasim.noisy_mag(zpicks, mu, sigma, data_params, data_key)
        mag, zpicks = results.load('./data', 'mag_z_LCDM_1000_sigma_0.05')
        
        firstderivs_functions = ['expgamma']
#        firstderivs_functions = ['expgamma','txgamma','zxgamma','gamma_over_z',
#                                 'zxxgamma','gammaxxz','rdecay_m','rdecay_de',
#                                 'rdecay_mxde','rdecay','interacting','LCDM']
        
        for key in firstderivs_functions:
            
            # Creating a folder for saving output.
            save_path = './quick_emcee/'+str(int(time.time()))+'_'+key
            if not os.path.exists(save_path):
                os.makedirs(save_path)    
    
            if key == 'LCDM':
                test_params = {'m':m_true}
            else:
                test_params = {'m':m_true, 'gamma':g_true}
            
            # Profiling
            import cProfile, pstats, io
            pr = cProfile.Profile()
            pr.enable()
            
            propert, sampler = paramfinder.paramfinder(npoints, nsteps, sigma, 
                                                       mu, test_params, zpicks, 
                                                       mag, save_path, key)
            
            pr.disable()
            s = io.StringIO()
            sortby = 'cumulative'
            ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
            ps.print_stats()
#            print (s.getvalue())
            
            f = open('profiler_quick_emcee.txt','w')
            f.write(s.getvalue())
            f.close()
            
        i += 1
    

    
    # Saving sampler to directory.
    results.save(save_path, 'sampler', sampler)

    print('Model being tested:', 
          test_key)
    print('Data simulated with:',data_key)

    return

quickemcee()


def errorvsdatasize():
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
    if not decision:
        pass
    else:
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
            propert, sampler = paramfinder.paramfinder(       
            npoints, nsteps, sigma, mu, test_params, zpicks, 
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