#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 12:40:48 2018

@author: BallBlueMeercat


"""
import matplotlib.pyplot as plt
import numpy as np
import os
import time
from results import load


def stat(hue, var, var_true, var_name, slnprob, zpicks, 
          mag, sigma, nsteps, nwalkers, save_path, firstderivs_key):
    
    initial = var_name.lower()[:1]
    name_true = var_name[:1] + '_true'
    hue_name = hue
    hue = 'xkcd:'+hue
    
    # Marginalised distribution histogram.
    plt.figure()
#    plt.xlabel(r'$\{}$'.format(name_l))
    plt.xlabel(var_name)
    plt.title('model: '+firstderivs_key+'\n Marginalised distribution of '
              +var_name+' \n nsteps: '+str(nsteps)+', noise: '+str(sigma)
              +', npoints: '+str(len(zpicks))+' '+firstderivs_key)
    plt.hist(var, 50, facecolor=hue)
    stamp = str(int(time.time()))
    filename = str(stamp)+'_'+initial+'_mhist__nsteps_'+str(nsteps) \
    +'_nwalkers_'+str(nwalkers)+'_noise_'+str(sigma) \
    +'_numpoints_'+str(len(zpicks))+'.png'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    
    # Walker steps.
    plt.figure()
    plt.xlabel(var_name)
    plt.title('model: '+firstderivs_key+'\n lnprobability of '+var_name
              +' \n nsteps: '+str(nsteps)+', noise: '+str(sigma)
              +', npoints: '+str(len(zpicks)))
    plt.plot(var, slnprob, '.', color=hue)
    stamp = str(int(time.time()))
    filename = str(stamp)+'_'+initial+'_steps__nsteps_'+str(nsteps) \
    +'_nwalkers_'+str(nwalkers)+'_noise_'+str(sigma) \
    +'_numpoints_'+str(len(zpicks))+'.png'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    
    # Chains.
    plt.figure()
    plt.xlabel('step number')
#    plt.ylabel(r'$\{}$'.format(name_l))
    plt.ylabel(var_name)
    plt.title('model: '+firstderivs_key+'\n flatchains, '+name_true+
              ' in '+hue_name+' \n nsteps: '+str(nsteps)+', noise: '
              +str(sigma)+', npoints: '+str(len(zpicks)))
    plt.plot(var.T, '-', color='k', alpha=0.3)
    plt.axhline(var_true, color=hue)
    stamp = str(int(time.time()))
    filename = str(stamp)+'_'+initial+'_chain__nsteps_'+str(nsteps) \
    +'_nwalkers_'+str(nwalkers)+'_noise_'+str(sigma) \
    +'_numpoints_'+str(len(zpicks))+'.png'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    
    plt.show(block=False)
    
    return


def precise_runs(firstderivs_key, params_dic, p, x):
    '''
    Takes in:
        firstderivs_key = string, model tested;
        params_dic = list of dicts, parameters used to generate 'LCDM' data;
        p = flot/int, cutoff precision in % for vc (params with mean !=0);
        x = flot/int, cutoff precision for s.d.
    '''
    
    
    # Results folder to search through.
    directory = os.path.join('./results_error_vs_data/'+firstderivs_key) 
    
    # Lists to be populated with contents of results folders. 
    sd_list = []
    mean_list = []
    vc_list = []
    sigma_list = []
    npoints_list = []
    
    # Creating a list of folders for each run that the model was tested. 
    folders = []
    for d in os.walk(directory):
        folders.append(d[0])
    folders.pop(0)
    
    # Colecting sd, mean, vc of marginalised distributions for all parameters, and 
    # the dataset properties such as sigma of noise added and number of points used.
    for folder in folders:
        sd_list += load(folder, 'sd_list.p')
        mean_list += load(folder, 'mean_list.p')
        vc_list += load(folder, 'vc_list.p')
        sigma_list += load(folder, 'sigma_list.p')
        npoints_list += load(folder, 'npoints_list.p')
    
    n_param = len(params_dic)   # How many parameters were fitted.
    n_run = int(len(vc_list) / n_param)   # How many runs were recorded.
    
    # For each parameter:
    for j in range(n_param):
        initial = vc_list[j][0][0] # First letter of the parameter fitted.
        
        sd = []
        mean = []
        vc = []
        # Results recorded for each run:
        for i in range(n_run):
            index = i*n_param+j # index of current parameter's results, 
                                # given number of parameters fitted
            
            sd.append(sd_list[index][1])
            mean.append(mean_list[index][1])
            vc.append(vc_list[index][1])            
            
            i+=1 # Onto the next run.
        
        # Converting to np.ndarrays to find & remove rows with repeating npoints.
        sd = np.array(sd)
        vc = np.asarray(vc)
        sigma = np.asarray(sigma_list)
        npoints = np.asarray(npoints_list)
    
        # Since input m and M are !=0, fitted m and M are unlikely to have mean=0, 
        # so precision can be based on the variance coefficient vc=sd/mean*100.
        if initial == 'm' or initial == 'M':
            # Indicies of rows with vc > p%.
            vc_above_p_index = np.where(abs(vc) > p)
            # Eliminating parameters found with precision > p%.
            p_stack = np.stack((npoints, sigma, vc), axis=1)
            p_stack = np.delete(p_stack, vc_above_p_index, 0)
            
            # Removing all but the noisiest run for each dataset size.
            for l in range(len(p_stack)-1):
                k = l+1 # next line
                npoints_l = p_stack[l][0] # Size of dataset on run l.
                npints_k = p_stack[k][0]  # Size of dataset on run k.
                if npoints_l == npints_k:   # If dataset sizes are the same, then
                    sigma_l = p_stack[l][1]   # compare sd of noise added to data.
                    sigma_k = p_stack[k][1]   # and leave the noisier run results.
                    if sigma_l > sigma_k:
                        p_stack = np.delete(p_stack, k, 0)
                    elif sigma_l < sigma_k:
                        p_stack = np.delete(p_stack, l, 0)
                l+=1
    
            # Splitting npoints, sigma added to data and variance 
            # coefficient into one dimentional np.arrays.        
            p_npoints, p_sigma, p_vc= np.hsplit(p_stack, 3)
            p_npoints = p_npoints.flatten()
            p_sigma = p_sigma.flatten()
            p_vc = p_vc.flatten()
            
            # Plotting dataset size vs noise added to data for all runs, and runs 
            # where parameters were found withing precision p, with vc annotated.
            fig, ax = plt.subplots()
            ax.scatter(npoints, sigma, c='c', label=('all runs')) 
            ax.scatter(p_npoints, p_sigma, c='m', label=('vc < %s%%'%(p)))
            # Annotating vc.
            for i, txt in enumerate(p_vc):
                txt = str(round(txt,2))
                ax.annotate(txt, (p_npoints[i], p_sigma[i]))
              
            plt.xlabel('Dataset size')
            plt.ylabel('Sigma of noise added to data')
            plt.title('Runs where '+initial+' was recovered within '+
                      str(p)+' percent. \n Variance coefficients annotated.')  
            plt.legend()      
            
        # Since input interaction terms are 0, recovered terms often have mean=0, 
        # hence precision is based on stadard deviation to avoid division by 0 
        # when calculating a variance coefficient.
        else:
            # Indicies of rows with sd > x.
            sd_above_x_index = np.where(abs(sd) > x)
            
            # Eliminating parameters found outside of sd < x.
            x_stack = np.stack((npoints, sigma, sd), axis=1)
            x_stack = np.delete(x_stack, sd_above_x_index, 0)
    
            # Removing all but the noisiest run for each dataset size.
            for l in range(len(x_stack)-1):
                k = l+1 # next line
                npoints_l = x_stack[l][0] # Size of dataset on run l.
                npints_k = x_stack[k][0]  # Size of dataset on run k.
                if npoints_l == npints_k:   # If dataset sizes are the same, then
                    sigma_l = x_stack[l][1]   # compare sd of noise added to data.
                    sigma_k = x_stack[k][1]   # and leave the noisier run results.
                    if sigma_l > sigma_k:
                        x_stack = np.delete(x_stack, k, 0)
                    elif sigma_l < sigma_k:
                        x_stack = np.delete(x_stack, l, 0)
                l+=1
    
            # Splitting npoints, sigma added to data and standard 
            # deviation into one dimentional np.arrays.        
            x_npoints, x_sigma, x_sd= np.hsplit(x_stack, 3)
            x_npoints = x_npoints.flatten()
            x_sigma = x_sigma.flatten()
            x_sd = x_sd.flatten()        
            
            # Plotting dataset size vs noise added to data for all runs, and runs 
            # where parameters were found with s.d. < x, with s.d. annotated.
            fig, ax = plt.subplots()
            ax.scatter(npoints, sigma, c='c', label=('all runs')) 
            ax.scatter(x_npoints, x_sigma, c='m', label=('s.d. < %s'%(x)))
            # Annotating vc.
            for i, txt in enumerate(x_sd):
                txt = str(round(txt,2))
                ax.annotate(txt, (x_npoints[i], x_sigma[i]))
              
            plt.xlabel('Dataset size')
            plt.ylabel('Sigma of noise added to data')
            plt.title('Runs where '+initial+' was recovered with s.d. < '+
                      str(x)+'. \n Standard deviations annotated.')  
            plt.legend()      
            
        j+=1    # Onto the next parameter.
    
    plt.show()
    
    return
      

def modelcheck(mag, zpicks, plot_var, firstderivs_key):
    
    t = plot_var.get('t')
    dl = plot_var.get('dl')
    a = plot_var.get('a')
    ombar_m = plot_var.get('ombar_m')
    ombar_r = plot_var.get('ombar_r')
    gamma = plot_var.get('gamma')
    zeta = plot_var.get('zeta')
    ombar_de = plot_var.get('ombar_de')
    
    t = -t
    
    if min(ombar_m) < 0:
        print('unphysical ombar_m', str(min(ombar_m)))
    elif min(ombar_de) < 0:
        print('unphysical ombar_de', str(min(ombar_de)))    
    
    # Scale factor vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel('a')
    plt.grid(True)
    plt.plot(zpicks, a, 'xkcd:crimson', lw=1)
    plt.title('Scale factor evolution, model = %s, $\gamma$ = %s, $\zeta$ = %s'
              %(firstderivs_key, gamma, zeta))
    
    Hz = []
    for i in range(len(ombar_m)):
        H = (ombar_m[i] + ombar_de[i])**(1/2)
        if np.isnan(H):
            print('plots.modelcheck got NaN value for H')
        Hz.append(H)

    # H vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel('H')
    plt.grid(True)
    plt.plot(zpicks, Hz, color='xkcd:blue', lw=1)
    plt.title(r'H evolution, model = %s, $\gamma$ = %s, $\zeta$ = %s'
              %(firstderivs_key, gamma, zeta))
    
    if ombar_r is None:
        # ombar_m, ombar_de vs redshift.
        plt.figure()
        plt.xlabel('redshift $z$')
        plt.ylabel(r'$\bar \Omega $')
        plt.grid(True)
        plt.plot(zpicks, ombar_m, label=r'$\bar \Omega_{m}$', 
                 color='xkcd:coral', linestyle=':')
        plt.plot(zpicks, ombar_de, label=r'$\bar \Omega_{DE}$', 
                 color='xkcd:aquamarine')
        plt.legend()
        plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, model = %s, $\gamma$ = %s'
              %(firstderivs_key, gamma))
        
        # ombar_m, ombar_de vs redshift log x axis.
        plt.figure()
        plt.xlabel('redshift $z$')
        plt.ylabel(r'$\bar \Omega $')
        plt.grid(True)
        plt.semilogx(zpicks, ombar_m, label=r'$\bar \Omega_{m}$', 
                     color='xkcd:coral', linestyle=':')
        plt.semilogx(zpicks, ombar_de, label=r'$\bar \Omega_{DE}$', 
                     color='xkcd:aquamarine')
        plt.legend()
        plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, model = %s, $\gamma$ = %s'
              %(firstderivs_key, gamma))
        
    else:
        # ombar_m, ombar_de, ombar_r evolution
        plt.figure()
        plt.xlabel('redshift $z$')
        plt.ylabel(r'$\bar \Omega $')
        plt.grid(True)
        plt.plot(zpicks, ombar_m, label=r'$\bar \Omega_{m}$', 
                 color='xkcd:coral', linestyle=':')
        plt.plot(zpicks, ombar_de, label=r'$\bar \Omega_{DE}$', 
                 color='xkcd:aquamarine')
        plt.plot(zpicks, ombar_r, label=r'$\bar \Omega_{r}$', 
                 color='xkcd:black')
        plt.legend()
        plt.title(r'$\bar \Omega_{r}$ evolution, model = %s, $\gamma$ = %s, $\zeta$ = %s'
              %(firstderivs_key, gamma, zeta)) 

    # Luminosity distance vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel('$d_L$*($H_0$/c)')
    plt.grid(True)     
    plt.plot(zpicks, dl, 'xkcd:lightgreen', lw=1)
    plt.title('$d_L$ evolution, model = %s, $\gamma$ = %s, $\zeta$ = %s'
              %(firstderivs_key, gamma, zeta))

    # Redshift vs -time.
    plt.figure()
    plt.xlabel('age')
    plt.ylabel('redshift $z$')
    plt.grid(True)
    plt.plot(t, zpicks, 'm', lw=1)
    plt.title('Redshift evolution, model = %s, $\gamma$ = %s, $\zeta$ = %s'
          %(firstderivs_key, gamma, zeta))

    # Scale factor vs -time.
    plt.figure()
    plt.xlabel('age')
    plt.ylabel('a')
    plt.grid(True)
    plt.plot(t, a, color='xkcd:crimson', lw=1)
    plt.title('Scale factor evolution, model = %s, $\gamma$ = %s, $\zeta$ = %s'
          %(firstderivs_key, gamma, zeta))
    
    # Magnitude vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel('magnitude')
    plt.title('Magnitude evolution, model = %s, $\gamma$ = %s, $\zeta$ = %s'
          %(firstderivs_key, gamma, zeta))
    plt.scatter(zpicks, mag, marker='.', lw='1', c='xkcd:tomato')

    plt.show()
    return


def multi_modelcheck(zpicks, firstderivs_key, plot_var_list, label):
    
    mag = []
    t = []
    dl = []
    a = []
    ombar_m = []
    ombar_r = []
    ombar_de = []
    ombar_m0 = []
    ombar_r0 = []
    ombar_de0 = []
    gamma = []
    zeta = []
    
    for i in range(len(plot_var_list)):
        for plot_var in plot_var_list:
            mag.append(plot_var.get('mag'))
            t.append(plot_var.get('t'))
            dl.append(plot_var.get('dl'))
            a.append(plot_var.get('a'))
            ombar_m.append(plot_var.get('ombar_m'))
            ombar_r.append(plot_var.get('ombar_r'))
            ombar_de.append(plot_var.get('ombar_de'))
            ombar_m0.append(plot_var.get('ombar_m0'))
            ombar_r0.append(plot_var.get('ombar_r0'))
            ombar_de0.append(plot_var.get('ombar_de0'))
            gamma.append(plot_var.get('gamma'))
            zeta.append(plot_var.get('zeta'))
    
    # Changing time into age
    t[0], t[1], t[2] = -t[0], -t[1], -t[2]
    
    # Scale factor vs redshift.
    fig, ax = plt.subplots()
    plt.xlabel('redshift $z$')
    plt.ylabel('a')
    plt.grid(True)
    ax.plot(zpicks, a[0], 'r-', label=label[0])
    ax.plot(zpicks, a[1], 'g-.', label=label[1])
    ax.plot(zpicks, a[2], 'b:', label=label[2])
    plt.title('Scale factor evolution, model = %s'%(firstderivs_key))
    ax.legend()
    
    # H vs redshift.
    Hz_1 = (ombar_m[0] + ombar_de[0])**(1/2)
    Hz_2 = (ombar_m[1] + ombar_de[1])**(1/2)
    Hz_3 = (ombar_m[2] + ombar_de[2])**(1/2)
    
    fig, ax = plt.subplots()
    plt.xlabel('redshift $z$')
    plt.ylabel('H')
    plt.grid(True)
    ax.plot(zpicks, Hz_1, 'r-', label=label[0])
    ax.plot(zpicks, Hz_2, 'g-.', label=label[1])
    ax.plot(zpicks, Hz_3, 'b:', label=label[2])
    plt.title('H evolution, model = %s'%(firstderivs_key))
    ax.legend()

    if ombar_r[0] is None:
        
        # ombar_m, ombar_de vs redshift.
        fig, ax = plt.subplots()
        plt.xlabel('redshift $z$')
        plt.ylabel(r'$\bar \Omega $')
        plt.grid(True)
        ax.plot(zpicks, ombar_m[0], 'r:', label='m')
        ax.plot(zpicks, ombar_m[1], 'g:', label='m')
        ax.plot(zpicks, ombar_m[2], 'b:', label='m')
        ax.plot(zpicks, ombar_de[0], 'm-.', label='de')
        ax.plot(zpicks, ombar_de[1], 'k-.', label='de')
        ax.plot(zpicks, ombar_de[2], 'c-.', label='de')
        plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, model = %s'
                  %(firstderivs_key))
        ax.legend()
    
        # ombar_m, ombar_de vs redshift log x axis.
        fig, ax = plt.subplots()
        plt.xlabel('redshift $z$')
        plt.ylabel(r'$\bar \Omega $')
        plt.grid(True)
        ax.semilogx(zpicks, ombar_m[0], 'r:', label='m')
        ax.semilogx(zpicks, ombar_m[1], 'g:', label='m')
        ax.semilogx(zpicks, ombar_m[2], 'b:', label='m')
        ax.semilogx(zpicks, ombar_de[0], 'm-.', label='de')
        ax.semilogx(zpicks, ombar_de[1], 'k-.', label='de')
        ax.semilogx(zpicks, ombar_de[2], 'c-.', label='de')
        plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, model = %s'
                  %(firstderivs_key))
        ax.legend()

    else:        
        # ombar_r, ombar_m, ombar_de vs redshift.
        fig, ax = plt.subplots()
        plt.xlabel('redshift $z$')
        plt.ylabel(r'$\bar \Omega $')
        plt.grid(True)
        ax.plot(zpicks, ombar_r[0], 'r:', label='r')
        ax.plot(zpicks, ombar_r[1], 'g:', label='r')
        ax.plot(zpicks, ombar_r[2], 'b:', label='r')
        ax.plot(zpicks, ombar_m[0], 'r-.', label='m')
        ax.plot(zpicks, ombar_m[1], 'g-.', label='m')
        ax.plot(zpicks, ombar_m[2], 'b-.', label='m')
        ax.plot(zpicks, ombar_de[0], 'm-', label='de')
        ax.plot(zpicks, ombar_de[1], 'k-', label='de')
        ax.plot(zpicks, ombar_de[2], 'c-', label='de')
        plt.title(r'$\bar \Omega_{r}$, $\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, model = %s'
                  %(firstderivs_key))
        ax.legend()
    
        # ombar_r, ombar_m, ombar_de vs redshift log x axis.
        fig, ax = plt.subplots()
        plt.xlabel('redshift $z$')
        plt.ylabel(r'$\bar \Omega $')
        plt.grid(True)
        ax.semilogx(zpicks, ombar_r[0], 'r:', label='r')
        ax.semilogx(zpicks, ombar_r[1], 'g:', label='r')
        ax.semilogx(zpicks, ombar_r[2], 'b:', label='r')        
        ax.semilogx(zpicks, ombar_m[0], 'r-.', label='m')
        ax.semilogx(zpicks, ombar_m[1], 'g-.', label='m')
        ax.semilogx(zpicks, ombar_m[2], 'b-.', label='m')
        ax.semilogx(zpicks, ombar_de[0], 'm-', label='de')
        ax.semilogx(zpicks, ombar_de[1], 'k-', label='de')
        ax.semilogx(zpicks, ombar_de[2], 'c-', label='de')
        plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, model = %s'%(firstderivs_key))
        ax.legend()
    
        # ombar_r vs redshift.
        fig, ax = plt.subplots()
        plt.xlabel('redshift $z$')
        plt.ylabel('$\Omega_{m}$')
        plt.grid(True)
        ax.plot(zpicks, ombar_r[0], 'r-', label=label[0])
        ax.plot(zpicks, ombar_r[1], 'g-.', label=label[1])
        ax.plot(zpicks, ombar_r[2], 'b:', label=label[2])
        plt.title('$\Omega_{m}$ evolution, model = %s'%(firstderivs_key))
        ax.legend()
	
    # Luminosity distance vs redshift.
    fig, ax = plt.subplots()
    plt.xlabel('redshift $z$')
    plt.ylabel('$d_L$*($H_0$/c)')
    plt.grid(True)
    ax.plot(zpicks, dl[0], 'r-', label=label[0])
    ax.plot(zpicks, dl[1], 'g-.', label=label[1])
    ax.plot(zpicks, dl[2], 'b:', label=label[2])
    plt.title('$d_L$ evolution, model = %s'%(firstderivs_key))
    ax.legend()

    # Redshift vs age.
    fig, ax = plt.subplots()
    plt.xlabel('age')
    plt.ylabel('redshift $z$')
    plt.grid(True)
    ax.plot(t[0], zpicks, 'r-', label=label[0])
    ax.plot(t[1], zpicks, 'g-.', label=label[1])
    ax.plot(t[2], zpicks, 'b:', label=label[2])
    plt.title('Redshift evolution, model = %s'%(firstderivs_key))
    ax.legend()
    
    # Scale factor vs -time.
    fig, ax = plt.subplots()
    plt.xlabel('age')
    plt.ylabel('a')
    plt.grid(True)
    ax.plot(t[0], a[0], 'r-', label=label[0])
    ax.plot(t[1], a[1], 'g-.', label=label[1])
    ax.plot(t[2], a[2], 'b:', label=label[2])
    plt.title('Scale factor evolution, model = %s'%(firstderivs_key))
    ax.legend()

    # Magnitude vs redshift.
    fig, ax = plt.subplots()
    plt.xlabel('redshift $z$')
    plt.ylabel('magnitude')
    plt.scatter(zpicks, mag[0], c='r', marker='.', lw=0.1, label=label[0])
    plt.scatter(zpicks, mag[1], c='g', marker='.', lw=0.1,  label=label[1])
    plt.scatter(zpicks, mag[2], c='b', marker='.', lw=0.2, label=label[2])
    plt.title('Magnitude evolution, model = %s'%(firstderivs_key))
    ax.legend()
    
    plt.show()
    return