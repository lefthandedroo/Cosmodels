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
import matplotlib as mpl
import tools
mpl.style.use('default') # has to be switched on to set figure size
mpl.style.use('fivethirtyeight')
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['grid.color'] = 'white'

test_keys_order_dic = {'rainbow':22,'niagara':21,'kanangra':20,'waterfall':19,'stepfall':18,'exotic':17,'late_intxde':16,'heaviside_late_int':15,'heaviside_sudden':14,'late_int':13,'expgamma':12,'txgamma':11,'zxgamma':10,'gamma_over_z':9,'zxxgamma':8,'gammaxxz':7,'rdecay_m':6,'rdecay_de':5,'rdecay_mxde':4,'rdecay':3,'interacting':2,'LCDM':0,'rLCDM':1}

# DNest4 output options sigma_synth = , sigma = 'pantheon', etc.
error_options = 'pantheon', 0.9, 0.3, 0.11, 0.07, 0.05, 0.001
#error_options = None, 0.001
dic_of_aic = {}
dic_of_bic = {}
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
    name_list = ['Model Name'] * len(test_keys_order_dic)


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

            name_list[order] = model_name
            ndim = str(len(names))
            name_list[order] = f'{model_name} {ndim}'


#    aic_array -= best_aic
#    bic_array -= best_bic
    dic_of_aic[error] = aic_array
    dic_of_bic[error] = bic_array


#    print(aic_array)
#    print(bic_array)

    txt_height = max(aic_array)
    label_pos = 0 # 9 centre top, 0 top left
    y_max = max(bic_array)+10


    barWidth = 0.25
    r1 = np.arange(len(test_keys_order_dic))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]

#    plt.figure()
#    fig, ax = plt.subplots()
#    ax.bar(r2, aic_array, alpha=0.5, width=barWidth, facecolor="C{}".format(0), edgecolor='white', label=r'AIC')# $\sigma$ = '+str(error))
#    ax.bar(r3, bic_array, alpha=0.5, width=barWidth, facecolor="C{}".format(1), edgecolor='white', label=r'BIC')# $\sigma$ = '+str(error))
#    for i, key in enumerate(name_list):
#        plt.text(i+barWidth, txt_height, key, {'ha': 'left', 'va': 'bottom','fontname':'monospace'}, rotation=90)
#    plt.xticks([])
##    ax.set_yscale('log')
#    plt.ylim(0, y_max)
#    plt.xlim(0, 23)
#    plt.legend(frameon=False, loc=label_pos, title='$\sigma =$'+str(error)) # 9 centre top, 0 top left
#    plt.show()

max_change = 0
plt.figure()
for key in dic_of_aic:
#    if key == 'pantheon':
#        continue
    residual = dic_of_aic['pantheon']-dic_of_aic[key]
    percentage = residual/dic_of_aic['pantheon'] * 100
    if max(abs(percentage)) > max_change:
        max_key = key
        max_change = max(abs(percentage))
        print(max_key, max_change)
    plt.plot(residual, label=key)
plt.legend(loc='lower right', ncol=3)
plt.show()

print('max change in aic = ',max_change, 'at sigma = ',key)









