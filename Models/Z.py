#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 12:25:25 2018

@author: usyd
"""

import numpy as np
import os
import results

rootdir = './results_Bfactor'
LCDM_log_Z = 353.97527367893576

def compare(LCDM_log_Z, model_name, model_log_Z):
    
    log_Bfactor = LCDM_log_Z - model_log_Z
    print()
    print('log(Bfactor) =',log_Bfactor)
    
    Bfactor = np.exp(log_Bfactor)
    
    print('(LCDM Z) / ('+model_name, 'Z) = ', str(Bfactor))
    print()
   
    return
    

for subdir, dirs, files in os.walk(rootdir):
    for directory in dirs:

        path = os.path.join(subdir, directory)
        model_log_Z = results.load(path, 'evidence.p') 
        model_name = directory[19:]
        if model_name == 'LCDM':
            pass
        else:
            compare(LCDM_log_Z, model_name, model_log_Z)