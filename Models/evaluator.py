#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 21:45:40 2018

@author: BallBlueMeercat
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import time
import os.path
import results
import tools
import datasim
import stats
import plots

# Number of emcee steps.
nsteps = 1000

# Statistical parameteres of noise:
mu = 0.0            # mean
sigma = 0.07        # standard deviation

# Pantheon data:
dataname = './data/lcparam_full_long.txt'
pantheon = pd.read_csv(dataname, sep=" ")

# Reading each txt file column of interest as numpy.ndarray
mag = pantheon.mb.values
x1 = pantheon.x1.values
colour = pantheon.color.values
zpicks = pantheon.zhel.values

# Stacking them together and sorting by accending redshift.
data = np.stack((mag,x1,colour,zpicks), axis=0)
data.sort(axis=-1)

mag = data[0]
#x1 = data[1]
#colour = data[2]
zpicks = data[3]
zpicks = zpicks.tolist()
data_dic = {'mag':mag, 'zpicks':zpicks}

#g1, g2, g3 = 0.0, 0.2, 0.3
#z1, z2, z3 = 0.0, -0.2, -0.4
## Compare param evolution for 3 models, plotting on the same axis.
#p1 = [{'matter':0.3},{'Mcorr':-19.3},{'alpha':0},{'beta':0},{'gamma':g1},{'zeta':z1}]
#p2 = [{'matter':0.3},{'Mcorr':-19.3},{'alpha':0},{'beta':0},{'gamma':g2},{'zeta':z2}]
#p3 = [{'matter':0.3},{'Mcorr':-19.3},{'alpha':0},{'beta':0},{'gamma':g3},{'zeta':z3}]
#
#datasim.multi_modelcheck([p1, p2, p3], zpicks, ['LCDM', 'late_int', 'exotic'],
#    ['$\gamma$='+str(g1)+' $\zeta$='+str(z1), 
#     '$\gamma$='+str(g2)+' $\zeta$='+str(z2), 
#     '$\gamma$='+str(g3)+' $\zeta$='+str(z3)])

firstderivs_functions = [None
            ,'waterfall'
#            ,'rainbow'
            ,'exotic'
            ,'late_intxde'
#            ,'heaviside_late_int'
#            ,'late_int'
#            ,'expgamma'
#            ,'txgamma'         # doesn't converge
#            ,'zxgamma'
#            ,'gamma_over_z'    # doesn't converge
#            ,'zxxgamma'        # gamma forced positive in firstderivs
#            ,'gammaxxz'        # gamma forced positive in firstderivs
#            ,'rdecay_m'
#            ,'rdecay_de'
#            ,'rdecay_mxde'
#            ,'rdecay'               
#            ,'interacting'
            ,'LCDM'
             ]

def modelcheck():    
    params_dic = [{'matter':0.3}, {'Mcorr':-19.3}, {'gamma':0.1}, {'zeta':0.2}]
    
    for test_key in firstderivs_functions:
        if test_key:
            datasim.magn(params_dic, data_dic, test_key, plot_key=True)
    return

#modelcheck()

    
def emcee():
    print('@@@@@@@ Mcor_emcee @@@@@@@')

    for test_key in firstderivs_functions:
            
        if test_key:
            print('---',test_key)
            if test_key == 'waterfall':
                params_dic = [{'matter':0.3},{'Mcorr':-19.3},
                              {'radiation':0.025},{'a_ombar':0.1},
                              {'b_ombar':0.1},{'c_ombar':0.1},
                              {'v_in':0.0},{'w_in':0.0},{'x_in':0.0},
                              {'y_in':0.0},{'z_in':0.0}]
            elif test_key == 'rainbow':
                params_dic = [{'matter':0.3},{'Mcorr':-19.3},{'radiation':0.025},
                              {'a_ombar':0.0},{'b_ombar':0.0},{'c_ombar':0.0},
                              {'d_ombar':0.0},{'e_ombar':0.0},{'f_ombar':0.0},
                              {'g_ombar':0.0},{'h_ombar':0.0},{'i_ombar':0.0},
                              {'q_in':0.0},{'r_in':0.0},{'s_in':0.0},{'t_in':0.0},
                              {'u_in':0.0}]
            elif test_key == 'exotic':
                params_dic = [{'matter':0.3},{'Mcorr':-19.3},{'gamma':0},{'zeta':0}]
            else:
                params_dic = [{'matter':0.3},{'Mcorr':-19.3},{'gamma':0}]
                
            
            # Creating a folder for saving output.
            save_path = './results_emcee/'+str(int(time.time()))+'_'+test_key
            if not os.path.exists(save_path):
                os.makedirs(save_path)    
          
            # Script timer.
            timet0 = time.time()            
            
            # emcee parameter search.
            propert, sampler = stats.stats(params_dic, data_dic, sigma, 
                                           nsteps, save_path, test_key, plot=1)        
            # Time taken by script. 
            timet1=time.time()
            tools.timer('script', timet0, timet1)
            
            # Saving sampler to directory.
            results.save(save_path, 'sampler', sampler)
    
            print('Data:',dataname)

    return

emcee()

def errorvsdatasize():
    
    test_key = 'exotic'
    params_dic = [{'matter':0.3},{'Mcorr':-19.3}, {'gamma':0.0}, {'zeta':0.0}]
    
    # Script timer.
    timet0 = time.time()
    
    sigma = 0.3
    sigma_max = 1
    sigma_step = 0.2
    npoints_min = 2750
    npoints_max = 4000
    npoints_step = 600
    
    # How many iterations have I signed up for?
    N = tools.runcount(sigma, sigma_max, sigma_step,
                   npoints_min, npoints_max, npoints_step)
    
    if input('Happy with '+str(N)+' iterations? (enter=yes) '):
        return
    
    # Folder for saving output.
    directory = str(int(time.time()))+'_'+test_key
    # Relative path of output folder.
    save_path = './results_error_vs_data/'+test_key+'/'+directory
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    
    run = 0
    
    sd_list = []
    mean_list = []
    vc_list = []
    
    sigma_list = []
    npoints_list = []
    sampler_list = []
    
    while sigma < sigma_max:

        npoints = npoints_min
        
        while npoints < npoints_max:
            print('_____________________ run number',run)
            
            zpicks = datasim.redshift_picks(0.005, 2.0, npoints)
            data_dic = {'zpicks':zpicks}
            
            # Magnitudes corresponding to simulated redshifts:
            mag = datasim.noisy_mag(mu, sigma, params_dic, data_dic, 'LCDM')            
            
            data_dic['mag']=mag
            
            propert, sampler = stats.stats(params_dic, data_dic, 
                                           sigma, nsteps, save_path, 
                                           test_key, plot=False)
            
            for key in propert:
                if 'sd' in key:
                    sd = propert.get(key,0)
                    sd_list.append([key, sd])
                elif 'mean' in key:
                    mean = propert.get(key,0)
                    mean_list.append([key, mean])
                    if mean != 0:
                        vc = sd/mean * 100
                        vc_list.append([key[0]+'_vc', vc])
            
            sigma_list.append(sigma)
            npoints_list.append(npoints)
            sampler_list.append(sampler)
            
            npoints += npoints_step
            run += 1
        
        sigma += sigma_step
        sigma = round(sigma, 2)

    for j in range(len(params_dic)):
        sd = []
        mean = []
        vc = []
        for i in range(N):
            index = i*len(params_dic)+j
            
            sd_name = sd_list[index][0]
            sd_initial = sd_name[0]
            sd.append(sd_list[index][1])
            
            mean_name = mean_list[index][0]
            mean_initial = mean_name[0]
            mean.append(mean_list[index][1])
            
            vc_name = vc_list[index][0]
            vc_initial = vc_name[0]
            vc.append(vc_list[index][1])            
            
            i+=1

        fig, ax = plt.subplots()
        ax.scatter(npoints_list, sd, c='r')
        
        # Plotting SD vs dataset size.
        for i, txt in enumerate(sigma_list):
            txt = 'sd = '+ str(txt)
            ax.annotate(txt, (npoints_list[i], sd[i]))
            
        plt.xlabel('Dataset size')
        plt.ylabel('s.d. of a marginalised distribution')
        plt.title(sd_name+' vs dataset size'+
                  '\n s.d. of noise labeled, model '+test_key)
        stamp = str(int(time.time()))
        filename = str(stamp)+'_sd_of_'+sd_initial+'_.png'
        filename = os.path.join(save_path, filename)
        plt.savefig(filename)
        
        # Plotting mean vs dataset size.
        fig, ax = plt.subplots()
        ax.scatter(npoints_list, mean, c='c')
        for i, txt in enumerate(sigma_list):
            txt = 'sd = '+ str(txt)
            ax.annotate(txt, (npoints_list[i], mean[i]))
            
        plt.xlabel('Dataset size')
        plt.ylabel('Mean of a marginalised distribution')
        plt.title(mean_name+' vs dataset size'+
                  '\n s.d. of noise labeled, model '+test_key)
        stamp = str(int(time.time()))
        filename = str(stamp)+'_mean_of_'+mean_initial+'_.png'
        filename = os.path.join(save_path, filename)
        plt.savefig(filename)
        
        # Plotting variance coefficient vs dataset size.
        if len(vc) == N:
            fig, ax = plt.subplots()
            ax.scatter(npoints_list, vc, c='g')
            for i, txt in enumerate(sigma_list):
                txt = 'sd = '+ str(txt)
                ax.annotate(txt, (npoints_list[i], vc[i]))
            
            plt.xlabel('Dataset size')
            plt.ylabel('s.d. /mean x100 of a marginalised distribution')
            plt.title(vc_name+' vs dataset size'+
                      '\n s.d. of noise labeled, model '+test_key)
            stamp = str(int(time.time()))
            filename = str(stamp)+'_vc_of_'+vc_initial+'_.png'
            filename = os.path.join(save_path, filename)
            plt.savefig(filename)
            
            j+=1
    
    plt.show()
        
#     Saving results to directory.
    results.save(save_path, 'vc_list', vc_list)
    results.save(save_path, 'sd_list', sd_list)
    results.save(save_path, 'mean_list', mean_list)
    
    results.save(save_path, 'sigma_list', sigma_list)
    results.save(save_path, 'npoints_list', npoints_list)
    results.save(save_path, 'sampler_list', sampler_list)
    
    print('directory:',directory)
    
    # Time taken by script. 
    timet1=time.time()
    tools.timer('errorvsdatasize', timet0, timet1)
    
    plots.precise_runs(test_key, params_dic, 9, 1.2)
    
    return vc_list, sd_list, mean_list, sigma_list, npoints_list, sampler_list

#vc, sd, mean, sigma, npoints, sampler = errorvsdatasize()

