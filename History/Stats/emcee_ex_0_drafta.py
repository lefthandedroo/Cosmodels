#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 14:28:29 2017

Geraint's exercise

chose m, c (parameters) for a straight line
from the line pick N points (N=3, 5, 50,10, 1000)
pick sigma (size of the noise)
randomly deviate(offset) points in y direction by using 
sigma*random number from normal distribution
sigma the same for all points
then define the likelihood use likelihood for dataset 
with gaussian error
Lookup how to write the eqution for a likelihood
(and then use log likelihood)
plug into emcee

draw a plot of c vs m displaying the walkers' walk
produce marginalised distribution - historgram 
+
for every m and c found - plot them together 

"""
import corner
import emcee
import matplotlib.pyplot as plt
import numpy as np


# Parameters

m_true = 1  # intercept
c_true = 3  # slope
N = 1000       # number of datapoints
sigma = 0.05 # standard deviation
mu = 0      # mean

# Predction of range possible for parameters
m_guess = 10
b_guess = 3

m = np.linspace(-m_guess, m_guess, 1)
b = np.linspace(-b_guess, b_guess, 1)
# Datapoints

x = np.random.rand(N)
y = m_true*x + c_true

#print('x is '+str(x))
#print('y is '+str(y))
plt.plot(x,y, 'gx')


# Noise

noise = np.random.normal(mu,sigma,N)
#print('noise is '+str(noise))


# Signal, offset by Gaussian noise

signal = y + noise
#print('signal is '+str(signal))
plt.plot(x,signal, 'rx')


# Log likelihood

def logLike(sigma, mu, x):
#    print('sigma is = '+str(sigma))
#    print('mu is = '+str(mu))
    logLike =np.array([0])
    for i in range(m):
        for j in range(b):
            for k in range (signal):
                logLike =logLike + (m*x + b - y)/sigma**2
    return logLike


# Emcee

ndim, nwalkers = 2, 4
p0 = [np.random.rand(ndim) for i in range(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers, ndim, logLike, args=[mu, x])
sampler.run_mcmc(p0, 1000)


# Corner plot
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
fig = corner.corner(samples, labels=["$m$", "$c$"],
                      truths=[m_true, c_true])
fig.savefig("triangle.png")












