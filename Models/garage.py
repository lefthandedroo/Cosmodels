#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""
import numpy as np


from scipy import stats
from emcee import PTSampler, EnsembleSampler
try:
    import matplotlib.pyplot as plt
except:
    plt = None


def lnlike_geom(p, y):
    """
    The probability distribution of the number of failures before the first
    success, supported on the set { 0, 1, 2, 3, ... }.
    See http://en.wikipedia.org/wiki/Geometric_distribution
    """
    return (y * np.log(1. - p) + np.log(p)).sum()


# Beta(0.1, 5) prior on geometrical distribution ``p`` parameter
# - VERY informative (but i am trying to replicate ``pymc`` result).
def lnprior_geom(p):
    return stats.beta.logpdf(p, 0.1, 5.)


def lnpost_geom(p, y):
    lnpr = lnprior_geom(p)
    if np.isinf(lnpr):
        return -np.inf
    else:
        return lnlike_geom(p, y) + lnpr


def lnlike_pois(lambda_, y):
    return stats.poisson.logpmf(y, lambda_).sum()


# Uniform(0, 1000) prior on Poisson distribution ``lambda`` parameter
def lnprior_pois(lambda_):
    return stats.uniform.logpdf(lambda_, 0., 1000.)


def lnpost_pois(lambda_, y):
    lnpr = lnprior_pois(lambda_)
    if np.isinf(lnpr):
        return -np.inf
    else:
        return lnlike_pois(lambda_, y) + lnpr


if __name__ == '__main__':
    y = np.array([0, 1, 2, 3, 8])
    ntemps = 20
    nwalkers = 50
    ndim = 1

    print ("Check prob. functions using EnsembleSampler")
    print ("Checking geometrical model...")
    sampler = EnsembleSampler(nwalkers, ndim, lnpost_geom, args=(y,))
    p0 = np.random.uniform(low=0, high=1.0, size=(nwalkers, ndim))
    for p, lnprob, lnlike in sampler.sample(p0, iterations=200):
        pass
    sampler.reset()
    for p, lnprob, lnlike in sampler.sample(p, iterations=2000, thin=10):
        pass

    if plt:
        plt.hist(sampler.flatchain[::10, :], bins=30, normed=True)
        plt.xlabel(r'p for geometrical model')
        plt.ylabel(r'p(p | data)')
        plt.show()

    print ("Checking Poisson model...")
    sampler = EnsembleSampler(nwalkers, ndim, lnpost_pois, args=(y,))
    p0 = stats.gamma.rvs(3., 2., size=nwalkers * ndim).reshape((nwalkers, ndim))
    for p, lnprob, lnlike in sampler.sample(p0, iterations=200):
        pass
    sampler.reset()
    for p, lnprob, lnlike in sampler.sample(p, iterations=2000, thin=10):
        pass

    if plt:
        plt.close()
        plt.hist(sampler.flatchain[::10, :], bins=30, normed=True)
        plt.xlabel(r'$\lambda$ for Poisson model')
        plt.ylabel(r'p($\lambda$ | data)')
        plt.show()

    print ("Estimating evidence for geometrical model using PTSampler")
    sampler = PTSampler(ntemps, nwalkers, ndim, lnlike_geom, lnprior_geom,
                        loglargs=(y,))
    p0 = np.random.uniform(low=0, high=1.0, size=(ntemps, nwalkers, ndim))
    for p, lnprob, lnlike in sampler.sample(p0, iterations=200):
        pass
    sampler.reset()
    for p, lnprob, lnlike in sampler.sample(p, lnprob0=lnprob,
                                            lnlike0=lnlike,
                                            iterations=2000, thin=10):
        pass

    # Check perfomance of Tl  (this code is from @farr)
    print ("Temperature swap acceptance rates are ")
    for b, rate in zip(sampler.betas, sampler.tswap_acceptance_fraction):
        print ('T = ', 1.0/b, ' accept = ', rate)
    print

    if plt is not None:
        plt.close()
        # Print a plot of the TI integrand:
        mean_logls = np.mean(sampler.lnlikelihood.reshape((ntemps, -1)), axis=1)
        betas = sampler.betas
        plt.plot(betas, betas*mean_logls) # \int d\beta <logl> = \int d\ln\beta \beta <logl>
        plt.xscale('log')
        plt.xlabel(r'$\beta$')
        plt.ylabel(r'$\beta \left\langle \ln L \right\rangle_\beta$')
        plt.title('Thermodynamic Integration Integrand')
        plt.show()

    logZ_geom, uncertainty = sampler.thermodynamic_integration_log_evidence()
    print ("Estimated log evidence for geometrical model = ", logZ_geom, "+/-", uncertainty)

    print ("Estimating evidence for Poisson model using PT")
    sampler = PTSampler(ntemps, nwalkers, ndim, lnlike_pois, lnprior_pois,
                        loglargs=(y,))
    p0 = stats.gamma.rvs(3., 2., size=nwalkers * ndim * ntemps).reshape((ntemps,
                                                                         nwalkers,
                                                                         ndim))
    for p, lnprob, lnlike in sampler.sample(p0, iterations=200):
        pass
    sampler.reset()
    for p, lnprob, lnlike in sampler.sample(p, lnprob0=lnprob,
                                            lnlike0=lnlike,
                                            iterations=2000, thin=10):
        pass

    # Check perfomance of Tl (this code is from @farr)
    print ("Temperature swap acceptance rates are ")
    for b, rate in zip(sampler.betas, sampler.tswap_acceptance_fraction):
        print ('T = ', 1.0/b, ' accept = ', rate)
    print

    if plt is not None:
        # Print a plot of the TI integrand:
        plt.close()
        mean_logls = np.mean(sampler.lnlikelihood.reshape((ntemps, -1)), axis=1)
        betas = sampler.betas
        plt.plot(betas, betas*mean_logls) # \int d\beta <logl> = \int d\ln\beta \beta <logl>
        plt.xscale('log')
        plt.xlabel(r'$\beta$')
        plt.ylabel(r'$\beta \left\langle \ln L \right\rangle_\beta$')
        plt.title('Thermodynamic Integration Integrand')
        plt.show()

    logZ_pois, uncertainty = sampler.thermodynamic_integration_log_evidence()
    print ("Estimated log evidence for Poisson model = ", logZ_pois, "+/-", uncertainty)

    print ("Bayes factor of geometrical to Poison model = ", np.exp(logZ_geom -
                                                                   logZ_pois))

###############################################################################

def lnprior(theta):
#    print(' lnprior has been called')   
#    print('lnprior speaking: theta = ',theta)
#    print('lnprior -- theta', theta)
    
    # Checking if theta is an array or a scalar.
    if len(theta) == 1:
        m = theta
        if 0 < m < 1 or m == 1:
            return 0.0
    elif len(theta) == 2:
        m, gamma = theta
        if (0 < m < 1 or m == 1) and abs(gamma) < 0.1:
            return 0.0
    elif len(theta) == 3:
        m, gamma, de = theta
        if (0 < m < 1 or m == 1) and abs(gamma) < 0.1 and (0 < de < 1):
            return 0.0
        
    return -np.inf



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