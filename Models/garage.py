#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""
import numpy as np

params = [{'matter':0.3},{'Mcorr':-19.3},{'alpha':0},
                   {'beta':0},{'gamma':10},{'zeta':20}]

# Packing up interaction terms:
in_terms = np.zeros((1, (len(params)-4))).flatten()

k = 0
for i in range(4,len(params)):
    for key in params[i]:
        in_terms[k] = params[i][key]
        k+=1
print(in_terms)       