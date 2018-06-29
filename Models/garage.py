#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""

from __future__ import print_function, division

import os
import sys
import six
import numpy as np

# import DNest4
from dnest4 import randh, wrap
import dnest4

# import model and data
from createdata import *

LN2PI = np.log(2.*np.pi)
LNSIGMA = np.log(sigma)

# define prior range values
cmin = -10.  # lower range on c (the same as the uniform c prior lower bound)
cmax = 10.   # upper range on c (the same as the uniform c prior upper bound)

mmu = 0.     # mean of the Gaussian prior on m
msigma = 10. # standard deviation of the Gaussian prior on m

# Create model class
class DNest4Model(object):
    """
    Specify the model in Python.
    """
    def __init__(self):
        """
        Parameter values *are not* stored inside the class
        """
        pass

    def from_prior(self):
        """
        Draw points from the prior distribution.

        Returns:
            :class:`np.ndarray`: a numpy array containing the parameter values.
        
        """

        m = np.random.normal(mmu, msigma)
        c = np.random.uniform(cmin,  cmax)

        return np.array([m, c])

    def perturb(self, params):
        """
        Perturb the current parameters by proposing a new position. This takes a numpy array of
        parameters as input, and modifies it in-place. In just perturbs one parameter at a time.

        Args:
            params (:class:`np.ndarray`): the current parameters (this is modified by the function)

        Returns:
            float: the natural logarithm of the Metropolis-Hastings proposal ratio.
        """
        logH = 0.0 # log of the Metropolis-Hastings prior x proposal ratio
        
        # randomly choose which parameter to perturb 
        which = np.random.randint(len(params))
        mmu = 0.
        msigma = 10.
        if which == 0:
            # update H for Gaussian prior
            logH -= -0.5*((params[which]-mmu)/msigma)**2
        params[which] += 1.*randh() # scale factor of 1. (randh is a heavy-tailed distribution to occasionally sample distant values, although this isn't really required in this case)

        if which == 0:
            # update H for Gaussian prior
            logH += -0.5*((params[which]-mmu)/msigma)**2
        else:
            # wrap c value so that it stays within the prior range
            params[which] = wrap(params[which], cmin, cmax)

        return logH

    def log_likelihood(self, params):
        """
        Gaussian sampling distribution.
        """
        m, c = params # unpack parameters
        
        norm = -0.5*M*LN2PI - M*LNSIGMA
        chisq = np.sum(((data - straight_line(x, m, c))/sigma)**2)
        return norm - 0.5*chisq

# Create a model object and a sampler
model = DNest4Model()
sampler = dnest4.DNest4Sampler(model, backend=dnest4.backends.CSVBackend(".", sep=" "))

# Set up the sampler. The first argument is max_num_levels
gen = sampler.sample(max_num_levels=30, num_steps=1000, new_level_interval=10000,
                     num_per_step=10000, thread_steps=100, num_particles=5, lam=10, beta=100)

# Do the sampling (one iteration here = one particle save)
for sample in enumerate(gen):
    pass

# Run the postprocessing to get marginal likelihood and generate posterior samples
logZdnest4, infogaindnest4, _ = dnest4.postprocess(plot=False);

postsamples = np.loadtxt('posterior_sample.txt')

print(six.u('Marginalised evidence is {}'.format(logZdnest4)))

print('Number of posterior samples is {}'.format(postsamples.shape[0]))

# plot posterior samples (if corner.py is installed)
try:
    import matplotlib as mpl
    mpl.use("Agg") # force Matplotlib backend to Agg
    import corner # import corner.py
except ImportError:
    sys.exit(1)

fig = corner.corner(postsamples, labels=[r"$m$", r"$c$"], truths=[m, c])
fig.savefig('DNest4.png')

#m_true = 0.3
#g_true = 0.
#
#def integrate_posterior_1D(lnprob, xlim, zpicks, 
#                           mag, sigma, test_key, ndim):
#    
#    func = lambda theta: np.exp(lnprob(theta, zpicks, 
#                                       mag, sigma, test_key, ndim))    
#    return quad(func, xlim[0], xlim[1])
#
#
#def integrate_posterior_2D(lnprob, xlim, ylim, zpicks, 
#                           mag, sigma, test_key, ndim):
#    
#    print('xlim',xlim)
#    print('ylim', ylim)
#    
#    func = lambda theta1, theta0: np.exp(lnprob([theta0, theta1], zpicks, 
#                                                mag, sigma, test_key, ndim))
#    return dblquad(func, xlim[0], xlim[1],
#                             lambda x: ylim[0], lambda x: ylim[1])
#
#
#def modeltest(npoints, nsteps, sigma, mu, 
#              zpicks, mag, test_key, save_path):
#    
#    print('__________ testing', test_key)
#    
#    if test_key == 'LCDM':
#        params = {'m':m_true}
#    else:
#        params = {'m':m_true, 'gamma':g_true}
#    
#    propert, sampler = paramfinder.paramfinder(
#            npoints, nsteps, sigma, mu, params, zpicks, 
#            mag, test_key, save_path)
#    
#    trace = propert.get('trace',)    
##    if not trace:
##        print ('modeltest got no trace for %s from paramfinder'%(test_key))
#        
#    return trace, sampler
#
#
#def Bfactor(npoints, nsteps, sigma, mu, data_params, data_key, M1_key, M2_key):
#    # Changing directory to dedicated folder for saving output.
#    save_path, directory = path()
#    
#    print('Generating magnitudes...')
#    mag, zpicks = data(mu, sigma, npoints, data_params, data_key)
#    
##    trace_1D, sampler_1D = modeltest(npoints, nsteps, sigma, mu,
##                         zpicks, mag, M1_key, save_path)
#    
#    trace_2D, sampler_2D = modeltest(npoints, nsteps, sigma, mu,
#                         zpicks, mag, M2_key, save_path)
#
##    print('__________ integration 1D')
##    # Integration:
##    xlim = trace_1D.min(0), trace_1D.max(0)
##    Z1, err_Z1 = integrate_posterior_1D(lnprob, xlim, 
##                                        zpicks, mag, sigma, M1_key, 1)
##    print("Z1 =", Z1, "+/-", err_Z1)
#
#    print('__________ integration 2D')
#    xlim, ylim = zip(trace_2D.min(0), trace_2D.max(0))
#    Z2, err_Z2 = integrate_posterior_2D(lnprob, xlim, ylim, 
#                                        zpicks, mag, sigma, M2_key, 2)
#    print("Z2 =", Z2, "+/-", err_Z2)
#
##    print("Bayes factor:", Z2 / Z1)  
##    print('Data is simulated using',data_key)
##    print()
##    print('directory:',directory)
#    
#    return sampler_2D
# 
##sampler_2D = Bfactor(npoints, nsteps, sigma, mu, 
##                     data_params, data_key, data_key, test_key)





#import numpy as np
#
## Choose the "true" parameters.
#m_true = -0.9594
#b_true = 4.294
#f_true = 0.534
#
## Generate some synthetic data from the model.
#N = 50
#x = np.sort(10*np.random.rand(N))
#yerr = 0.1+0.5*np.random.rand(N)
#y = m_true*x+b_true
#y += np.abs(f_true*y) * np.random.randn(N)
#y += yerr * np.random.randn(N)
#
#A = np.vstack((np.ones_like(x), x)).T
#C = np.diag(yerr * yerr)
#cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
#b_ls, m_ls = np.dot(cov, np.dot(A.T, np.linalg.solve(C, y)))
#
#
#def lnlike(theta, x, y, yerr):
#    m, b, lnf = theta
#    model = m * x + b
#    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
#    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))
#
#import scipy.optimize as op
#nll = lambda *args: -lnlike(*args)
#result = op.minimize(nll, [m_true, b_true, np.log(f_true)], args=(x, y, yerr))
#m_ml, b_ml, lnf_ml = result["x"]
#
#ndim, nwalkers = 3, 6
#pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
#
#print('pos',pos)



from pylab import figure, plot, xlabel, ylabel, title, show
import matplotlib.pyplot as pl

import paramfinder
# Model parameteres:  
m_true = 0.3           # (= e_m(t)/e_crit(t0) at t=t0).
de_true = 1 - m_true   # (de = e_de(t)/e_crit(t0) at t=t0).
gamma_true = 0.0       # Interaction term, rate at which DE decays into matter.

# Number of datapoints to be simulated and number of emcee steps.
npoints, nsteps = 10000, 10000

# Statistical parameteres:
mu = 0          # mean
sigma = 0.085     # standard deviation

def stepevaluator():
    
    steps = []
    standev = []
    meanlist = []
    cvlist = []
    
    nsteps = 100
    run = 0
    while nsteps < 2000:
        print('_____________________ run number',run)
        steps.append(nsteps)
        propert, sampler = paramfinder.paramfinder(npoints, nsteps, sigma, mu, m_true)
        sd, mean = propert
        standev.append(sd)
        meanlist.append(mean)
        cv = sd/mean
        cvlist.append(cv)
        
        cv = sd / mean     # Coefficient of variation.
#        print('cv:',str(cv))
#        if cv < 0.008:
#            print('nsteps', nsteps)
#            break
        
        nsteps += 50
        run += 1
    
    figure()
    xlabel('emcee steps')
    ylabel('standard deviation')
    plot(steps, standev, '.')
    title('standard deviation of m found vs steps taken')
    show()
    
    return steps, sampler, standev, meanlist

#steps, sampler, standev, meanlist = stepevaluator()


def nevaluator():
    
    numpoints = []
    standev = []
    meanlist = []
    cvlist = []
    
    npoints = 100
    run = 0
    while npoints < 35000:    #35000
        print('_____________________ run number',run)
        numpoints.append(npoints)
        propert, sampler = paramfinder.paramfinder(npoints, nsteps, sigma, mu, m_true)
        sd, mean = propert
        standev.append(sd)
        meanlist.append(mean)
        cv = sd/mean
        cvlist.append(cv)
        npoints += 1000
        run += 1

    
    figure()
    xlabel('number of datapoints used')
    ylabel('standard deviation')
    plot(numpoints, standev, 'b.', label='standard deviation')
    plot(numpoints, cvlist, 'r.', label='coefficient of variation')
    pl.legend()
    title('sd and cv of m found vs dataset size')
    show()
    
    return numpoints, sampler, standev, meanlist

#numpoints, sampler, standev, meanlist = nevaluator()


def errevaluator():
    
    error = []
    standev = []
    meanlist = []
    cvlist = []
    
    sigma = 1e-9
    run = 0
    while sigma < 0.1:
        print('_____________________ run number',run)
        error.append(sigma)
        propert, sampler = paramfinder.paramfinder(npoints, nsteps, sigma, mu, m_true)
        sd, mean = propert
        standev.append(sd)
        meanlist.append(mean)
        cv = sd/mean
        cvlist.append(cv)
        sigma += 0.003
        run += 1
    
    figure()
    xlabel('standard deviation of noise')
    ylabel('standard deviation of parameter distribution')
    plot(error, standev, 'b.', label='standard deviation')
    plot(error, cvlist, 'r.', label='coefficient of variation')
    pl.legend()
    title('sd and cv of m found vs sd of noise added')
    show()
    
    return error, sampler, standev, meanlist

#error, sampler, standev, meanlist = errevaluator()