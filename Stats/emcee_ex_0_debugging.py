#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 17:12:08 2017

@author: BallBlueMeercat
"""
import numpy as np

def lnlike(theta, x, y, yerr):
    m, b = theta
    model = m * x + b
    inv_sigma2 = 1.0/(yerr**2 + model**2)
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))


def lnprior(theta):
    m, b = theta
    if -5.0 < m < 0.5 and 0.0 < b < 10.0:
        return 0.0
    return -np.inf   

def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)    

m_true = 0.1
b_true = 3
theta = [m_true, b_true]

N = 50
sigma = 0.05
mu = 0 

x = np.random.rand(N)
yerr = np.random.normal(mu,sigma,N)
y = m_true*x+b_true
y += yerr   

lnlike(theta, x, y, yerr)
