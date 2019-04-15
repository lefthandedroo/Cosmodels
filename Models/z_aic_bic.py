#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 19:43:00 2019

@author: BallBlueMeercat
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import results
import matplotlib as mpl
import tools
#mpl.style.use('default') # has to be switched on to set figure size
mpl.style.use('fivethirtyeight')
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['grid.color'] = 'white'

test_keys_order_dic = {'rainbow':22,'niagara':21,'kanangra':20,'waterfall':19,'stepfall':18,'exotic':17,'late_intxde':16,'heaviside_late_int':15,'heaviside_sudden':14,'late_int':13,'expgamma':12,'txgamma':11,'zxgamma':10,'gamma_over_z':9,'zxxgamma':8,'gammaxxz':7,'rdecay_m':6,'rdecay_de':5,'rdecay_mxde':4,'rdecay':3,'interacting':2,'LCDM':0,'rLCDM':1}

# DNest4 output options sigma_synth = , sigma = 'pantheon', etc.
#error_options = 'pantheon', 0.9, 0.3, 0.11, 0.07, 0.05, 0.001
error_options = None, 0.001
txt_height = 25
label_pos = 0 # 9 centre top, 0 top left
print_marker = False

for error in error_options:

    if error == 'pantheon':
        rootdir = './results_Bfactor/recent/'+str(error)
    elif not error:
        continue
    else:
        rootdir = './results_Bfactor/recent/'+str(error)+'_synth'
    bfactor_max = 0
    best_aic = 1e100
    aic_dic = {}
    best_bic = 1e100
    bic_dic = {}
    aic_array = np.zeros(len(test_keys_order_dic))
    bic_array = np.zeros(len(test_keys_order_dic))
    bfactor_array = np.ones(len(test_keys_order_dic))
    name_list = ['Model Name'] * len(test_keys_order_dic)

    # Getting LCDM log(Z)
    for subdir, dirs, files in os.walk(rootdir):
        for directory in dirs:
            model_name = directory[19:]
            if model_name == 'LCDM':
                path = os.path.join(subdir, directory)
                LCDM_log_Z = results.load(path, 'evidence.p')
                break

    for subdir, dirs, files in os.walk(rootdir):
        for directory in dirs:
            path = os.path.join(subdir, directory)
            model_name = directory[19:]
            order = test_keys_order_dic[model_name]

            names, values = tools.names_values(model_name)
            ndim = len(names)


            # Getting aic and bic
            sample_info = pd.read_csv(path+'/sample_info.txt', sep=" ")
            log_like = sample_info.log_likelihood.values
            max_log_like = max(log_like)

            aic = 2*ndim -2*max_log_like
            aic_dic[model_name] = aic
            aic_array[order] = aic
            if aic < best_aic:
                best_aic = aic
                best_aic_model = model_name
            bic = np.log(1048)*ndim -2*max_log_like
            bic_dic[model_name] = bic
            bic_array[order] = bic
            if bic < best_bic:
                best_bic = bic
                best_bic_model = model_name

            # Getting Z and bayes factor
            model_log_Z = results.load(path, 'evidence.p')
            if print_marker:
                print(model_name, model_log_Z)
            log_Bfactor = LCDM_log_Z - model_log_Z
            Bfactor = np.exp(log_Bfactor)
            bfactor_array[order] = Bfactor
            name_list[order] = model_name
            ndim = str(len(names))
            if Bfactor > bfactor_max:
                bfactor_max = Bfactor
            if np.isinf(Bfactor):
                name_list[order] = f'{model_name} {ndim} infinite'
                bfactor_array[order] = 10e250
            else:
                name_list[order] = f'{model_name} {ndim}'
                bfactor_array[order] = Bfactor
            if print_marker:
                print('log(Bfactor) =',log_Bfactor)
                print('(LCDM Z) / ('+model_name, 'Z) = ', str(Bfactor))
                print()

    aic_array -= best_aic
    bic_array -= best_bic

    x_fill = np.arange(-1, 30, 1)


    barWidth = 0.25
    r1 = np.arange(len(test_keys_order_dic))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]

    plt.figure()
    threshold = 1
    fig, ax = plt.subplots()
    ax.bar(r1[1:], bfactor_array[1:], alpha=0.5, width=barWidth, edgecolor='white', label=r'Bfactor')# $\sigma$ = '+str(error))
    ax.bar(r2[1:], aic_array[1:], alpha=0.5, width=barWidth, edgecolor='white', label=r'AIC')# $\sigma$ = '+str(error))
    ax.bar(r3[1:], bic_array[1:], alpha=0.5, width=barWidth, edgecolor='white', label=r'BIC')# $\sigma$ = '+str(error))
    for i, key in enumerate(name_list):
        if i > 0:
            plt.text(i-0.35, txt_height, key, {'ha': 'left', 'va': 'bottom','fontname':'monospace'}, rotation=90)
    y_max = bfactor_max+bfactor_max/2
    # strong evidence for alternative
#    ax.fill_between(x_fill, 0, 0.1, facecolor='crimson', alpha=0.7)
    # moderate evidence for alternative
#    ax.fill_between(x_fill, 0.1, 0.3, facecolor='crimson', alpha=0.5)
    # weak evidence for alternative
#    ax.fill_between(x_fill, 0.3, 1, facecolor='crimson', alpha=0.2)
    plt.axhline(y=1,lw=1, color='k')
    # weak evidence for LCDM
#    ax.fill_between(x_fill, 1, 3, facecolor='mediumseagreen', alpha=0.1)
    # moderate evidence for LCDM
#    ax.fill_between(x_fill, 3, 10, facecolor='mediumseagreen', alpha=0.3)
    # strong evidence for LCDM
#    ax.fill_between(x_fill, 10, y_max, facecolor='mediumseagreen', alpha=0.5)
    plt.xticks([])
    ax.set_yscale('log')
    plt.ylim(0.15, y_max)
    plt.xlim(0.5, 23)
    plt.legend(loc=label_pos) # 9 centre top, 0 top left
    plt.show()