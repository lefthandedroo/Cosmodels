#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 12:25:25 2018

@author: usyd
"""

import numpy as np

def compare(LCDM_log_Z, sigma, model_name, model_log_Z):
    
    Bfactor = LCDM_log_Z - model_log_Z
    print()
    print(model_name)
    print('log(Bfactor) =',Bfactor)
    
    Bfactor = np.exp(Bfactor)
    
#    print('sigma =', sigma)
    print('(LCDM Z) / ('+model_name, 'Z) = ', str(Bfactor))
    print()
    
#    Bfactor =  model_log_Z - LCDM_log_Z
#    
##    print('log(Bfactor) =',Bfactor)
#    
#    Bfactor = np.exp(Bfactor)
#    
#    print(model_name,'Z / LCDM Z = ', str(Bfactor))
    
    return Bfactor
    
#compare(-32844.59417770218, 0.07, 'late_intxde', -32844.55643113441)
#compare(-32844.59417770218, 0.07, 'heaviside_late_int', -31040.911807978944)
#compare(-32844.59417770218, 0.07, 'late_int', -31041.24720478364)
#compare(-32844.59417770218, 0.07, 'expgamma', -24232.177248562988)
#compare(-32844.59417770218, 0.07, 'txgamma', 325.3981263878408)
#compare(-32844.59417770218, 0.07, 'zxgamma', -20298.745374150178)
#compare(-32844.59417770218, 0.07, 'gamma_over_z', 229.16420577388232)
#compare(-32844.59417770218, 0.07, 'zxxgamma', -24243.51620038459)
#compare(-32844.59417770218, 0.07, 'gammaxxz', -19860.093792940486)
#compare(-32844.59417770218, 0.07, 'rdecay_de', -621.290963339889)
#compare(-32844.59417770218, 0.07, 'rdecay', -32844.30784323075)
#compare(-32844.59417770218, 0.07, 'rdecay_de', -621.290963339889)
compare(-32844.50826284908, 0.07, 'not real model', -32847.50826284908)