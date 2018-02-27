#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
chose m, c (parameters) for a straight line
from the line pick N points (N=3, 5, 50, 100, 1000)
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
for every m and c found - plot them together 
plot data and error bars, 
plot what the actual model is
find max likelihood
and m and b corresponding to max L
draw the line that they give
try for different sigmas


modify to find the points of the max posterior distribution
use those to plot the best line
but increase number of dimensions, think of curve that requires 
4-5 parameters
(say polynomial)
do a multi dimensional search, need more walkers
and more steps

Another thing - grab hold of cosmological equations
revise general relativity
need to use the 2nd order friedmann equations + the continuity eqns
numerically integrate the cosmological equations for standard 
cosmology where 
omega = 1 for matter universe
also
omega = 1 30% matter, 70% dark energy

Also read paper The Influence of Evolving Dark Energy on Cosmology

Try to UNDERSTAND
"""
import corner
import emcee
import matplotlib.pyplot as pl
import numpy as np
import scipy.optimize as op


# Input
    
# "True" parameters.
m_true = 1  # intercept
b_true = 3  # slope

N = 50          # number of datapoints
sigma = 0.5     # standard deviation
mu = 0          # mean

ndim, nwalkers = 2, 4
nsteps = 10000
burnin = 1000


# Functions

def lnlike(theta, x, y, sigma):
    m, b = theta
    model = m * x + b
    inv_sigma2 = 1.0/(sigma**2)
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))    
 
def lnprior(theta):
    m, b = theta
    if -5.0 < m < 5 and 0.0 < b < 10.0:
        return 0.0
    return -np.inf
        
def lnprob(theta, x, y, sigma):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, sigma)    

   

# Generating noisy data from the model y.
x = np.random.rand(N)               # picking random points on x-axis
yerr = np.random.normal(mu,sigma,N) # Gaussian noise
y = m_true*x+b_true                 # model (straight line)
y += yerr                           # data, offset in y with noise


# Finding a "good" place to start using alternative method to emcee.
nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [m_true, b_true], args=(x, y, yerr))
m_ml, b_ml = result["x"]    

    
# Initializing walkers in a Gaussian ball around the max likelihood. 
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]    
    

# Sampler setup
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, sigma))
sampler.run_mcmc(pos, nsteps)


# Mean acceptance fraction. In general, acceptance fraction has an entry 
# for each walker so, in this case, it is a 50-dimensional vector.
print("Mean acceptance fraction:", np.mean(sampler.acceptance_fraction))


# Corner plot (walkers' walk + histogram).
samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
fig = corner.corner(samples, labels=["$m$", "$b$"], truths=[m_true, b_true])
fig.savefig("triangle.png")


# Marginalised distribution (histogram) plot.
pl.hist(sampler.flatchain[:,0], 100)
pl.show()


## Plotting lines of best fit using m, b parameters found in each step.
#i = 0
#while i < (nwalkers*nsteps):
##    print('flatchain[[i],[0]] is = '+str(sampler.flatchain[[i],[0]]))
##    print('flatchain[[i],[1]] is = '+str(sampler.flatchain[[i],[1]]))
#    y_best_fit = sampler.flatchain[[i],[0]]*x + sampler.flatchain[[i],[1]]
#    i += 1
#    pl.plot(y_best_fit, x)
#pl.show()


# Best line of fit found by emcee.
bi = np.argmax(sampler.lnprobability)   # "best index" - index with highest 
                                        # posterior probability
                                        
mbest = sampler.flatchain[bi,0]         #parameters with the highest 
bbest = sampler.flatchain[bi,1]         # posterior probability

print('best index is = '+str(bi))
print('mbest is = '+str(mbest))
print('bbest is = '+str(bbest))

# plot of data with errorbars + model
pl.errorbar(x, y, yerr=sigma, fmt='o')
xt = np.linspace(0,1,10)
yt = m_true * xt + b_true
model, = pl.plot(xt,yt,lw='3', c='g')
ybest = mbest * xt + bbest
best_fit, = pl.plot(xt,ybest,lw='3', c='r')
pl.legend([model, best_fit], ['Model', 'Best Fit'])

pl.show








