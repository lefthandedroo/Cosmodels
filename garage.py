#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""

import numpy as np
from scipy import integrate

def df(f,t):
   return f-2.

def df_stop(f,t):
   return f < 0.0

f0 = 1.
t0 = 0.
t_max = 5.
nout = 100
ts = np.linspace(t0,t_max,nout)


fs = [f0,]
df_continue = True
i = 0
while df_continue:
    f = integrate.odeint(df,fs[i],[ts[i],ts[i+1]])
    i+=1
    if i==nout-1:
        df_continue = False
    elif df_stop(f[1][0],ts[i+1]):
        df_continue = False
    else:
        fs.append( f[1][0] )

fs = np.array( fs )
