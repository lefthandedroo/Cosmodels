#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 16:17:52 2018

@author: BallBlueMeercat
"""
import numpy as np
import lnprior

# prior:  if -0.1 < gamma < 0.1 and 0.299 < m < 0.301 and 0.699 < de < 0.701:

a = [0.05, 0.01, 0, -10]
b = [0.3, 0.3009, -10, 0.2999]
c = [0.7, -10, 0.6999, 0.7009]
prob = [0.0, -np.inf, -np.inf, -np.inf]

def check(gamma, m, de, slnprob):
    theta = np.column_stack((gamma, m, de))
    
    badtheta_finite_slnprob = 0
    goodtheta_inf_slnprob = 0
    
    badtheta_inf_slnprob = 0
    goodtheta_finite_slnprob = 0
    
    for i in range(len(slnprob)):
        lp = lnprior.lnprior(theta[i])
        
        if not np.isfinite(lp):             # if theta is outside prior
            if np.isfinite(slnprob[i]):     # and has finite slnprob
                print('theta = %s, slnprob = %s'%(theta,slnprob))
                badtheta_finite_slnprob +=1
                
            if not np.isfinite(slnprob[i]): # and has inf slnprob as it should
                badtheta_inf_slnprob +=1
                
        if np.isfinite(lp):                 # if theta is within prior
            if not np.isfinite(slnprob[i]): # but has inf slnprob
                print('theta = %s, slnprob = %s'%(theta,slnprob))
                goodtheta_inf_slnprob +=1
            
            if np.isfinite(slnprob[i]):     # & has finite slnprob as it should
                goodtheta_finite_slnprob +=1

    wrong = badtheta_finite_slnprob + goodtheta_inf_slnprob
    right = badtheta_inf_slnprob + goodtheta_finite_slnprob
#    print('Bad theta with finite slnprob found: ',badtheta_finite_slnprob)
#    print('Good theta with inf slnprob found: ',goodtheta_inf_slnprob)
#    print('Bad theta with inf slnprob found: ',badtheta_inf_slnprob)
    print('Good theta with finite slnprob found: ',goodtheta_finite_slnprob)
#    print('wrong = %s, right = %s'%(wrong, right))
    
#check(a, b, c, prob)