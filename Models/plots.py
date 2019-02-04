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

import matplotlib as mpl
mpl.style.use('default') # has to be switched on to set figure size
mpl.style.use('fivethirtyeight')
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['grid.color'] = 'white'

def stat_emcee(hue, var, var_true, var_name, slnprob, zpicks,
          mag, sigma, nsteps, nwalkers, save_path, firstderivs_key):
    mpl.style.use('fivethirtyeight')
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
    plt.plot(var.T, '-', color='k', alpha=0.3, lw=2)
    plt.axhline(var_true, color=hue, lw=2)
    stamp = str(int(time.time()))
    filename = str(stamp)+'_'+initial+'_chain__nsteps_'+str(nsteps) \
    +'_nwalkers_'+str(nwalkers)+'_noise_'+str(sigma) \
    +'_numpoints_'+str(len(zpicks))+'.png'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)

    plt.show(block=False)

    return

def stat_DNest(hue, var, var_true, var_name, slnprob, zpicks,
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

def precise_runs(firstderivs_key, names, values, p, x):
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

    n_param = len(values)   # How many parameters were fitted.
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
    Hz = plot_var.get('Hz')
    fluid_names = plot_var.get('fluid_names')
    fluid_arr = plot_var.get('fluid_arr')
    int_terms = plot_var.get('int_terms')
    da = plot_var.get('da')
    dV = plot_var.get('dV')

    t = -t

    if np.isnan(Hz).any():
        print('plots.modelcheck got NaN value for H')

    # Fluid evolution.
    plt.figure()
    plt.xlabel('$z$')
    plt.ylabel(r'$\bar \Omega $')
    plt.grid(True)
    for i in range(len(fluid_arr)):
        plt.plot(zpicks, fluid_arr[i], label=fluid_names[i])
    plt.legend()
    plt.title(r'$\bar \Omega_{X}$ evolution'
              f'\n Model: {firstderivs_key}, int_terms = {int_terms}')

    # Evolution of the angular diameter distance.
    plt.figure()
    plt.title('Angular diameter distance evolution'
              f'\n Model: {firstderivs_key}, int_terms = {int_terms}')
    plt.xlabel('z')
    plt.ylabel(r'$ d_A (H_0 / c)$')
    plt.grid(True)
    da_index = np.argmin(da)
    if da[da_index] < 0:
        plt.plot(zpicks[:,da_index], da[:, da_index])
        plt.plot(zpicks[da_index, :], da[da_index, :])
    else:
        plt.plot(zpicks, da)

    # Evolution of dV.
    plt.figure()
    plt.title(r'$d_V$ evolution'
              f'\n Model: {firstderivs_key}, int_terms = {int_terms}')
    plt.xlabel('z')
    plt.ylabel(r'$ d_V (z)$ [Mpc]')
    plt.grid(True)
    plt.plot(zpicks, dV)

    # Luminosity distance vs redshift.
    plt.figure()
    plt.xlabel('$z$')
    plt.ylabel('$d_L$*($H_0$/c)')
    plt.grid(True)
    plt.plot(zpicks, dl)
    plt.title('$d_L$ evolution'+'\n Model: %s, int_terms = %s'
              %(firstderivs_key, int_terms))

    # H vs redshift.
    plt.figure()
    plt.xlabel('$z$')
    plt.ylabel('H')
    plt.grid(True)
    plt.plot(zpicks, Hz)
    plt.title(r'$H(z)$ evolution'+'\n Model: %s, int_terms = %s'
              %(firstderivs_key, int_terms))

    # Scale factor vs redshift.
    plt.figure()
    plt.xlabel('$z$')
    plt.ylabel('a')
    plt.grid(True)
    plt.plot(zpicks, a)
    plt.title('Scale factor evolution with redshift '
              +'\n Model: %s, int_terms = %s'
              %(firstderivs_key, int_terms))

    # Scale factor vs age.
    plt.figure()
    plt.xlabel('Age')
    plt.ylabel('a')
    plt.grid(True)
    plt.plot(t, a)
    plt.title('Scale factor evolution with time'+'\n Model: %s, int_terms = %s'
          %(firstderivs_key, int_terms))

    # Redshift vs age.
    plt.figure()
    plt.xlabel('Age')
    plt.ylabel('$z$')
    plt.grid(True)
    plt.plot(t, zpicks)
    plt.title('Redshift evolution'+'\n Model: %s, int_terms = %s'
          %(firstderivs_key, int_terms))

    # Magnitude vs redshift.
    plt.figure()
    plt.xlabel('$z$')
    plt.ylabel('Magnitude')
    plt.title('Magnitude evolution'+'\n Model: %s, int_terms = %s'
          %(firstderivs_key, int_terms))
    plt.scatter(zpicks, mag, marker='.')

    plt.show()
    return


def multi_modelcheck(data, keys, plot_var_list):
    print('multi_modelcheck')

#    plt.rcParams['axes.facecolor'] = 'white'
#    plt.rcParams['figure.facecolor'] = 'white'
#    plt.rcParams['grid.color'] = 'white'
#    plt.rcParams['axes.edgecolor'] = 'black' # frame

    zpicks = data['zpicks']
    data_mag = data['mag']

    model_names = ''
    for model_name in keys:
        model_names += model_name+', '
    model_names = model_names[:-2]

    mag = []
    t = []
    dl = []
    a = []
    Hz = []
    da = []
    dV = []
    fluid_names = []
    fluid_arr = []
    int_terms = []

    for plot_var in plot_var_list:
        mag.append(plot_var.get('mag'))
        t.append(plot_var.get('t'))
        dl.append(plot_var.get('dl'))
        a.append(plot_var.get('a'))
        Hz.append(plot_var.get('Hz'))
        da.append(plot_var.get('da'))
        dV.append(plot_var.get('dV'))
        fluid_names.append(plot_var.get('fluid_names'))
        fluid_arr.append(plot_var.get('fluid_arr'))
        int_terms.append(plot_var.get('int_terms'))

    # Changing time t into age
    age=[]
    for item in t:
        item = -item
        age.append(item)

    # Fluid evolution.
    plt.figure()
    plt.xlabel('$z$')
    plt.ylabel(r'$\bar \Omega $')
    for i in range(len(fluid_arr)):
        for j in range(len(fluid_arr[i])):
            if j % 2 == 0:
                linestyle=(0, ())
                plt.plot(zpicks, fluid_arr[i][j], lw=2, color="C{}".format(i), ls=linestyle, label='$\gamma$ = %s'%(int_terms[i]))
            else:
                linestyle=(0, (3, 1, 1, 1))
                plt.plot(zpicks, fluid_arr[i][j], lw=2, color="C{}".format(i), ls=linestyle)
#    plt.title(r'$\bar \Omega_{X}$, models: %s.'%(model_names))
    plt.legend()

    # log x axis fluid evolution vs redshift.
    plt.figure()
    plt.xlabel('$z$')
    plt.ylabel(r'$\bar \Omega $')
    for i in range(len(fluid_arr)):
        for j in range(len(fluid_arr[i])):
            if j % 2 == 0:
                linestyle=(0, ())
                plt.semilogx(zpicks, fluid_arr[i][j], lw=2, color="C{}".format(i), ls=linestyle, label='$\gamma$ = %s'%(int_terms[i]))
            else:
                linestyle=(0, (3, 1, 1, 1))
                plt.semilogx(zpicks, fluid_arr[i][j], lw=2, color="C{}".format(i), ls=linestyle)
#    plt.title(r'$\bar \Omega_{X}$, models: %s.'%(model_names))
    plt.legend()

#    # Evolution of the angular diameter distance.
#    plt.figure()
##    plt.title('Angular diameter distance evolution'+f'\n Models: {model_names}')
#    plt.xlabel('z')
#    plt.ylabel(r'$ (H_0/c) \ D_A $')
#    for i in range(len(da)):
#        plt.plot(zpicks, da[i], label='$\gamma$ = %d, %s'%(int_terms[i], keys[i]))
#    plt.legend()

#     # Evolution of dV.
#    plt.figure()
##    plt.title(r'$d_V$ evolution, models: %s.'%(model_names))
#    plt.xlabel('z')
#    plt.ylabel(r'$ d_V (z)$ [Mpc]')
#    for i in range(len(dV)):
#        plt.plot(zpicks, dV[i], label='$\gamma$ = %d, %s'%(int_terms[i], keys[i]))
#    plt.legend()

#    # Luminosity distance vs redshift.
#    plt.figure()
#    plt.xlabel('$z$')
#    plt.ylabel('$D_L$*($H_0$/c)')
#    for i in range(len(dl)):
#        plt.plot(zpicks, dl[0], label='$\gamma$ = %d, %s'%(int_terms[i], keys[i]))
##    plt.title('$d_L$ evolution, models: %s.'%(model_names))
#    plt.legend()
#
#    # H vs redshift.
#    plt.figure()
#    plt.xlabel('$z$')
#    plt.ylabel('H')
#    for i in range(len(Hz)):
#        plt.plot(zpicks, Hz[i], label='$\gamma$ = %d, %s'%(int_terms[i], keys[i]))
##    plt.title('H evolution, models: %s.'%(model_names))
#    plt.legend()
#
#    # Scale factor vs redshift.
#    plt.figure()
#    plt.xlabel('$z$')
#    plt.ylabel('a')
#    for i in range(len(a)):
#        plt.plot(zpicks, a[0], label='$\gamma$ = %d, %s'%(int_terms[i], keys[i]))
##    plt.title('Scale factor evolution, models: %s.'%(model_names))
#    plt.legend()
#
#    # Scale factor vs age.
#    plt.figure()
#    plt.xlabel('Age')
#    plt.ylabel('a')
#    for i in range(len(age)):
#        plt.plot(age[i], a[i], label='$\gamma$ = %d, %s'%(int_terms[i], keys[i]))
##    plt.title('Scale factor evolution, models: %s.'%(model_names))
#    plt.legend()
#
#    # Redshift vs age.
#    plt.figure()
#    plt.xlabel('Age')
#    plt.ylabel('$z$')
#    for i in range(len(age)):
#        plt.plot(age[i], zpicks, label='$\gamma$ = %d'%(int_terms[i]))
##    plt.title('Redshift evolution, models: %s.'%(model_names))
#    plt.legend()

    # Magnitude vs redshift.
    plt.figure()
    plt.xlabel('$z$')
    plt.ylabel('Magnitude')

    for i in range(len(mag)):
        plt.plot(zpicks, mag[i], lw=2, label='$\gamma$ = %s'%(int_terms[i]))
#    plt.title('Magnitude evolution, models: %s.'%(model_names))
    plt.scatter(zpicks, data_mag, s=50, marker='o', c='darkslategrey', alpha=0.2, label='Pantheon')
    plt.legend()

    plt.show()
    return