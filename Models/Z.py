#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 12:25:25 2018

@author: usyd
"""

import numpy as np
import os
import results

# Data_name options 'man, 'sn', 'compare'.
data_name = 'man'     # fitted to aritifical LCDM data
#data_name = 'compare' # fitted to Pantheon

rootdir = os.path.join('./results_Bfactor', data_name)
LCDM_log_Z = None

def compare(LCDM_log_Z, model_name, model_log_Z):

    log_Bfactor = LCDM_log_Z - model_log_Z
    print('log(Bfactor) =',log_Bfactor)

    Bfactor = np.exp(log_Bfactor)

    print('(LCDM Z) / ('+model_name, 'Z) = ', str(Bfactor))
    print()

    return

for subdir, dirs, files in os.walk(rootdir):
    for directory in dirs:
        model_name = directory[19:]
        if model_name == 'LCDM':
            path = os.path.join(subdir, directory)
            LCDM_log_Z = results.load(path, 'evidence.p')
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
            print(model_name, model_log_Z)
            compare(LCDM_log_Z, model_name, model_log_Z)
