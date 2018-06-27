#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 20:20:03 2018

@author: BallBlueMeercat
"""
import time
import dnest4
import numpy as np
import numpy.random as rng
from datasim import magn
from results import load
from tools import timer
#from scipy.special import erf

firstderivs_key = 'LCDM'
sigma = 0.01

# Load the data
mag, zpicks = load('./data', 'mag_z_LCDM_1000_sigma_0.01')

g_max = 0.1



class Model(object):
    """
    Specify the model in Python.
    """
#    def __init__(self):
#        """
#        Parameter values *are not* stored inside the class
#        """
#        pass

    def from_prior(self):
        """
        Unlike in C++, this must *return* a numpy array of parameters.
        """
        m = rng.rand()
        g = rng.uniform(-g_max, g_max)
        
        return np.array([m, g])

    def perturb(self, params):
        """
        Unlike in C++, this takes a numpy array of parameters as input,
        and modifies it in-place. The return value is still logH.
        """
        logH = 0.0
        which = rng.randint(2)
        
#        print(params[which])
        if which == 0:
            log_m = np.log(params[which])
            log_m += dnest4.randh()
            # Note the difference between dnest4.wrap in Python and
            # DNest4::wrap in C++. The former *returns* the wrapped value.
            log_m = dnest4.wrap(log_m, 0, 1)
            params[which] = np.exp(log_m)
            
        elif which == 1:
            g = params[which]
            g += dnest4.randh()
            # Note the difference between dnest4.wrap in Python and
            # DNest4::wrap in C++. The former *returns* the wrapped value.
            g = dnest4.wrap(g, -g_max, g_max)
            params[which] = g

        return logH

    def log_likelihood(self, params):
        """
        Gaussian sampling distribution.
        """
        m, g = params

        theta = {'m':m,'gamma':g}
        
        model = magn(theta, zpicks, firstderivs_key)
        
        var = sigma**2
        return -0.5*np.sum((mag-model)**2 /var +0.5*np.log(2*np.pi*var))

# Create a model object and a sampler
model = Model()
sampler = dnest4.DNest4Sampler(model,
                               backend=dnest4.backends.CSVBackend(".",
                                                                  sep=" "))

# Set up the sampler. The first argument is max_num_levels
#gen = sampler.sample(max_num_levels=30, num_steps=1000, new_level_interval=10000,
#                      num_per_step=10000, thread_steps=100,
#                      num_particles=5, lam=10, beta=100, seed=1234)

# num_per_step can be down to a few thousand
gen = sampler.sample(max_num_levels=30, num_steps=1000, new_level_interval=1000,
                      num_per_step=1000, thread_steps=100,
                      num_particles=5, lam=10, beta=100, seed=1234)

ti = time.time()
# Do the sampling (one iteration here = one particle save)
for i, sample in enumerate(gen):
    print("# Saved {k} particles.".format(k=(i+1)))
tf = time.time()

timer('Sampling', ti, tf)

# Run the postprocessing
dnest4.postprocess()

# LCDM
#log(Z) = -1622866.8534441872
#Information = 14.078678027261049 nats.
#Effective sample size = 129.22232212112772
#time 297min 50s

#log(Z) = -1622866.790641218
#Information = 13.905435690656304 nats.
#Effective sample size = 167.73507536834273
#time 34 min

#rdecay
#log(Z) = -1622866.8177826053
#Information = 13.970533838961273 nats.
#Effective sample size = 85.54638980461822
#Sampling time:   37min 5s