#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 12:25:25 2018

@author: usyd
"""

import numpy as np

def compare(LCDM_log_Z, sigma, model_name, model_log_Z):
    
    Bfactor = LCDM_log_Z - model_log_Z
    
    print('log(Bfactor) =',Bfactor)
    
    Bfactor = np.exp(Bfactor)
    
#    print('sigma =', sigma)
    print('LCDM Z / ', model_name, 'Z = ', str(Bfactor))
    
#    Bfactor =  model_log_Z - LCDM_log_Z
#    
##    print('log(Bfactor) =',Bfactor)
#    
#    Bfactor = np.exp(Bfactor)
#    
#    print(model_name,'Z / LCDM Z = ', str(Bfactor))
    
    return Bfactor
    
compare(-32844.50826284908, 0.07, 'expgamma', -24232.425212751157)
compare(-32844.50826284908, 0.07, 'rdecay', -32844.57535456529)
compare(-32844.50826284908, 0.07, 'txgamma', -23339.420419805938)
compare(-32844.50826284908, 0.07, 'zxgamma', -20298.791888243453)
#compare(-32844.50826284908, 0.07, 'not real model', -32847.50826284908)