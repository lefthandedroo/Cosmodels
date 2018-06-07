#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""
import numpy as np
data = np.array([[ 0.42,  0.72,  0.  ,  0.3 ,  0.15,
                   0.09,  0.19,  0.35,  0.4 ,  0.54,
                   0.42,  0.69,  0.2 ,  0.88,  0.03,
                   0.67,  0.42,  0.56,  0.14,  0.2  ],
                 [ 0.33,  0.41, -0.22,  0.01, -0.05,
                  -0.05, -0.12,  0.26,  0.29,  0.39, 
                   0.31,  0.42, -0.01,  0.58, -0.2 ,
                   0.52,  0.15,  0.32, -0.13, -0.09 ],
                 [ 0.1 ,  0.1 ,  0.1 ,  0.1 ,  0.1 ,
                   0.1 ,  0.1 ,  0.1 ,  0.1 ,  0.1 ,
                   0.1 ,  0.1 ,  0.1 ,  0.1 ,  0.1 ,
                   0.1 ,  0.1 ,  0.1 ,  0.1 ,  0.1  ]])


import matplotlib.pyplot as plt
import seaborn as sns; sns.set() # set default plot styles

x, y, sigma_y = data
fig, ax = plt.subplots()
ax.errorbar(x, y, sigma_y, fmt='ok', ecolor='gray')
ax.set(xlabel='x', ylabel='y', title='input data');

def polynomial_fit(theta, x):
    """Polynomial model of degree (len(theta) - 1)"""
    return sum(t * x ** n for (n, t) in enumerate(theta))


from scipy import stats

def logL(theta, model=polynomial_fit, data=data):
    """Gaussian log-likelihood of the model at theta"""
    x, y, sigma_y = data
    y_fit = model(theta, x)
    return sum(stats.norm.logpdf(*args)
               for args in zip(y, y_fit, sigma_y))

from scipy import optimize

def best_theta(degree, model=polynomial_fit, data=data):
    theta_0 = (degree + 1) * [0]
    neg_logL = lambda theta: -logL(theta, model, data)
    return optimize.fmin_bfgs(neg_logL, theta_0, disp=False)

theta1 = best_theta(1)
theta2 = best_theta(2)

xfit = np.linspace(0, 1, 1000)
fig, ax = plt.subplots()
ax.errorbar(x, y, sigma_y, fmt='ok', ecolor='gray')
ax.plot(xfit, polynomial_fit(theta1, xfit), label='best linear model')
ax.plot(xfit, polynomial_fit(theta2, xfit), label='best quadratic model')
ax.legend(loc='best', fontsize=14)
ax.set(xlabel='x', ylabel='y', title='data');

degrees = np.arange(1, 10)
thetas = [best_theta(d) for d in degrees]
logL_max = [logL(theta) for theta in thetas]

fig, ax = plt.subplots(1, 2, figsize=(14, 5))
ax[0].plot(degrees, logL_max)
ax[0].set(xlabel='degree', ylabel='log(Lmax)')
ax[1].errorbar(x, y, sigma_y, fmt='ok', ecolor='gray')
ylim = ax[1].get_ylim()
for (degree, theta) in zip(degrees, thetas):
    if degree not in [1, 2, 9]: continue
    ax[1].plot(xfit, polynomial_fit(theta, xfit),
               label='degree={0}'.format(degree))
ax[1].set(ylim=ylim, xlabel='x', ylabel='y')
ax[1].legend(fontsize=14, loc='best');


def log_prior(theta):
    # size of theta determines the model.
    # flat prior over a large range
    if np.any(abs(theta) > 100):
        return -np.inf  # log(0)
    else:
        return 200 ** -len(theta)

def log_likelihood(theta, data=data):
    x, y, sigma_y = data
    yM = polynomial_fit(theta, x)
    return -0.5 * np.sum(np.log(2 * np.pi * sigma_y ** 2)
                         + (y - yM) ** 2 / sigma_y ** 2)

def log_posterior(theta, data=data):
    theta = np.asarray(theta)
    return log_prior(theta) + log_likelihood(theta, data)

import emcee

def compute_mcmc(degree, data=data,
                   log_posterior=log_posterior,
                   nwalkers=50, nburn=1000, nsteps=2000):
    ndim = degree + 1  # this determines the model
    rng = np.random.RandomState(0)
    starting_guesses = rng.randn(nwalkers, ndim)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[data])
    sampler.run_mcmc(starting_guesses, nsteps)
    trace = sampler.chain[:, nburn:, :].reshape(-1, ndim)
    return trace

trace_2D = compute_mcmc(1)
trace_3D = compute_mcmc(2)

import pandas as pd
columns = [r'$\theta_{0}$'.format(i) for i in range(3)]
df_2D = pd.DataFrame(trace_2D, columns=columns[:2])

with sns.axes_style('ticks'):
    jointplot = sns.jointplot(r'$\theta_0$', r'$\theta_1$',
                              data=df_2D, kind="hex");
                              
df_3D = pd.DataFrame(trace_3D, columns=columns[:3])

# get the colormap from the joint plot above
cmap = jointplot.ax_joint.collections[0].get_cmap()

with sns.axes_style('ticks'):
    grid = sns.PairGrid(df_3D)
    grid.map_diag(plt.hist, bins=30, alpha=0.5)
    grid.map_offdiag(plt.hexbin, gridsize=50, linewidths=0, cmap=cmap)
    
from scipy import integrate

def integrate_posterior_2D(log_posterior, xlim, ylim, data=data):
    func = lambda theta1, theta0: np.exp(log_posterior([theta0, theta1], data))
    return integrate.dblquad(func, xlim[0], xlim[1],
                             lambda x: ylim[0], lambda x: ylim[1])

def integrate_posterior_3D(log_posterior, xlim, ylim, zlim, data=data):
    func = lambda theta2, theta1, theta0: np.exp(log_posterior([theta0, theta1, theta2], data))
    return integrate.tplquad(func, xlim[0], xlim[1],
                             lambda x: ylim[0], lambda x: ylim[1],
                             lambda x, y: zlim[0], lambda x, y: zlim[1])

xlim, ylim = zip(trace_2D.min(0), trace_2D.max(0))
Z1, err_Z1 = integrate_posterior_2D(log_posterior, xlim, ylim)
print("Z1 =", Z1, "+/-", err_Z1)

xlim, ylim, zlim = zip(trace_3D.min(0), trace_3D.max(0))
Z2, err_Z2 = integrate_posterior_3D(log_posterior, xlim, ylim, zlim)
print("Z2 =", Z2, "+/-", err_Z2)
              
print("Bayes factor:", Z2 / Z1)






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