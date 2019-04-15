#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 12:25:25 2018

@author: usyd
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import matplotlib as mpl
import tools
#mpl.style.use('default') # has to be switched on to set figure size
mpl.style.use('fivethirtyeight')
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['grid.color'] = 'white'

test_keys_order_dic = {'rainbow':22,'niagara':21,'kanangra':20,'waterfall':19,'stepfall':18,'exotic':17,'late_intxde':16,'heaviside_late_int':15,'heaviside_sudden':14,'late_int':13,'expgamma':12,'txgamma':11,'zxgamma':10,'gamma_over_z':9,'zxxgamma':8,'gammaxxz':7,'rdecay_m':6,'rdecay_de':5,'rdecay_mxde':4,'rdecay':3,'interacting':2,'LCDM':0,'rLCDM':1}

error_list = 'pantheon', 0.9, 0.3, 0.11, 0.07, 0.05, 0.001
#error_list = 0.05, 0.001
for error in error_list:
    if error == 'pantheon':
        rootdir = './results_Bfactor/recent/'+str(error)
    else:
        rootdir = './results_Bfactor/recent/'+str(error)+'_synth'
    best_aic = 1e100
    aic_dic = {}
    best_bic = 1e100
    bic_dic = {}
    aic_list = [None] * len(test_keys_order_dic)
    bic_list = [None] * len(test_keys_order_dic)
    name_list = ['Model Name'] * len(test_keys_order_dic)
#    print('\n', error)
    for subdir, dirs, files in os.walk(rootdir):
        for directory in dirs:
            path = os.path.join(subdir, directory)
            model_name = directory[19:]
            name_list[test_keys_order_dic[model_name]] = model_name

            sample_info = pd.read_csv(path+'/sample_info.txt', sep=" ")
            log_like = sample_info.log_likelihood.values
            max_log_like = max(log_like)

            names, values = tools.names_values(model_name)
            ndim = len(names)

            aic = 2*ndim -2*max_log_like
            aic_dic[model_name] = aic
            aic_list[test_keys_order_dic[model_name]] = aic
            if aic < best_aic:
                best_aic = aic
                best_aic_model = model_name
            bic = np.log(1048)*ndim -2*max_log_like
            bic_dic[model_name] = bic
            bic_list[test_keys_order_dic[model_name]] = bic
            if bic < best_bic:
                best_bic = bic
                best_bic_model = model_name
#            print(model_name, aic)
#    print('best aic',best_aic_model, best_aic)
#    print('best bic',best_bic_model, best_bic)

    x_fill = np.arange(-1, 30, 1)
    y_max = 0

    plt.figure()
    number = range(len(aic_list))
    threshold = 1
    fig, ax = plt.subplots()
    ax.bar(number, aic_list, alpha=0.5, label=r'$\sigma$ = '+str(error))
    for i, key in enumerate(name_list):
        plt.text(i-0.35, 10, key, {'ha': 'left', 'va': 'bottom','fontname':'monospace'}, rotation=90)
        if aic_list[i] > y_max:
            y_max = aic_list[i]
    y_max = y_max+y_max/2
    # strong evidence for alternative
    ax.fill_between(x_fill, 0, 0.1, facecolor='crimson', alpha=0.7)
    # moderate evidence for alternative
    ax.fill_between(x_fill, 0.1, 0.3, facecolor='crimson', alpha=0.5)
    # weak evidence for alternative
    ax.fill_between(x_fill, 0.3, 1, facecolor='crimson', alpha=0.2)
    # weak evidence for LCDM
    ax.fill_between(x_fill, 1, 3, facecolor='mediumseagreen', alpha=0.1)
    # moderate evidence for LCDM
    ax.fill_between(x_fill, 3, 10, facecolor='mediumseagreen', alpha=0.3)
    # strong evidence for LCDM
    ax.fill_between(x_fill, 10, y_max, facecolor='mediumseagreen', alpha=0.5)
    plt.xticks([])
#    ax.set_yscale('log')
    plt.ylim(0, y_max)
    plt.xlim(-0.6, 22.6)
    plt.legend(loc=0) # 9 centre top, 0 top left
    plt.show()