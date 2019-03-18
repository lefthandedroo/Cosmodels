#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 12:25:25 2018

@author: usyd
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import results
import matplotlib as mpl
#mpl.style.use('default') # has to be switched on to set figure size
mpl.style.use('fivethirtyeight')
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['grid.color'] = 'white'

# Data_name options 'synth, 'amanullah', 'pantheon', 'recent'.
rootdir = './results_Bfactor/recent/0.07'

print_marker = False

ratio_dic ={'LCDM':1}
LCDM_log_Z = None

def compare(LCDM_log_Z, model_name, model_log_Z):

    log_Bfactor = LCDM_log_Z - model_log_Z


    Bfactor = np.exp(log_Bfactor)

    ratio_dic[model_name] = Bfactor

    if print_marker:
        print('log(Bfactor) =',log_Bfactor)
        print('(LCDM Z) / ('+model_name, 'Z) = ', str(Bfactor))
        print()

    return

for subdir, dirs, files in os.walk(rootdir):
    for directory in dirs:
        model_name = directory[19:]
        if model_name == 'LCDM':
            path = os.path.join(subdir, directory)
            LCDM_log_Z = results.load(path, 'evidence.p')
            if print_marker:
                print(model_name, LCDM_log_Z)
                print()
            break

for subdir, dirs, files in os.walk(rootdir):
    for directory in dirs:
        path = os.path.join(subdir, directory)
        model_name = directory[19:]
        if model_name == 'LCDM':
            pass
        else:
            model_log_Z = results.load(path, 'evidence.p')
            if print_marker:
                print(model_name, model_log_Z)
            compare(LCDM_log_Z, model_name, model_log_Z)



objects = ratio_dic.keys()
y_pos = np.arange(len(objects))
performance = ratio_dic.values()
threshold = 1
fig, ax = plt.subplots()
ax.bar(y_pos, performance, alpha=0.5)
plt.axhline(y=0.1,lw=1, color="C{}".format(0), label='strong evidence for alternative')
plt.axhline(y=0.3,lw=1, color="C{}".format(1), label='moderate evidence for alternative')
plt.axhline(y=1,lw=1, color='k', label='weak evidence for alternative')
plt.axhline(y=3,lw=1, color="C{}".format(2), label='weak evidence for LCDM')
plt.axhline(y=10,lw=1,color="C{}".format(3), label='moderate evidence for LCDM')
plt.axhline(y=10,lw=1, color="C{}".format(4), label='strong evidence for LCDM')

#plt.xticks(y_pos, objects)
ax.set_yscale('log')
plt.legend()
plt.show()
