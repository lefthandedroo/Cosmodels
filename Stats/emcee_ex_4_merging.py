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
modify to find parameters with max posterior distribution
use those to plot the best line
increase number of dimensions, think of curve that requires 4-5 parameters
(say polynomial)
do a multi dimensional search, need more walkers
and more steps

look up first two erros -  whether they are to do with python version

Try to UNDERSTAND

Notable results:
    best index is = 3010835
    abest is = 3.9419662932
    bbest is = -3.01946040697
    cbest is = 0.990232737609
    dbest is = 15.0034779775
    ebest is = 1.50005168141
    Mean acceptance fraction: 0.50711475
    Number of steps: 100000
    Number of walkers: 200
    Sampler time: 63min 42s
    Total time:   65min 30s
"""
import corner
import emcee
import logging
import matplotlib.pyplot as pl
import numpy as np
import scipy.optimize as op
import sys
import time


try:
    timet0 = time.time()    # starting script timer
    
    
    # Input
    # "True" parameters.
    a_true = 0.1
    b_true = -3
    c_true = 0.5
    d_true = 0.1
    e_true = 12
    
    N = 20        # number of datapoints
    sigma = 0.75    # standard deviation
    mu = 0          # mean
    
    ndim, nwalkers = 5, 12
    nsteps = 1000
    burnin = 500
    
    
    # Functions
    def lnlike(theta, x, y, sigma):
        a, b, c, d, e = theta
        model = a * x**4 + b * x**2 + c * x + d + e*np.sin(x)
        inv_sigma2 = 1.0/(sigma**2)
        return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))    
     
    def lnprior(theta):
        a, b, c, d, e = theta
        if (-5.0 < a < 5 and -5.0 < b < 5.0 and 0.0 < c < 1.0 and 0.0 < d < 20 
            and -3.0 < e < 30):
            return 0.0
        return -np.inf
            
    def lnprob(theta, x, y, sigma):
        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike(theta, x, y, sigma)    
       
    
    # Generating noisy data from the model y.
    x = np.random.rand(N)*4                 # picking random points on x-axis
    yerr = np.random.normal(mu,sigma,N)     # Gaussian noise
    y = a_true * x**4 + b_true * x**2 + c_true * x + d_true + e_true*np.sin(x)     
    y += yerr                               # data, offset in y with noise
    
    
    # Finding a "good" place to start using alternative method to emcee.
    nll = lambda *args: -lnlike(*args)
    result = op.minimize(nll, [a_true, b_true, c_true, d_true, e_true], 
                         args=(x, y, yerr))
    a_ml, b_ml, c_ml, d_ml, e_ml = result["x"]    
    
        
    # Initializing walkers in a Gaussian ball around the max likelihood. 
    pos = [result["x"] + 1*np.random.randn(ndim) for i in range(nwalkers)]    
        
    
    # Sampler setup
    times0 = time.time()    # starting emcee timer
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, sigma))
    sampler.run_mcmc(pos, nsteps)
    
    times1=time.time()      # stopping emcee timer
    times=times1 - times0   # time to run emcee
    timesmin = round((times / 60),1)    # minutes
    timessec = round((times % 60),1)    # seconds
    
    
    # Corner plot (walkers' walk + histogram).
    samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
    fig = corner.corner(samples, labels=["$a$", "$b$", "$c$", "$d$", "$e$"], 
                        truths=[a_true, b_true, c_true, d_true, e_true])
    fig.savefig('nsteps'+str(nsteps)+str(time.strftime("%c"))+
                'nwalkers'+str(nwalkers)+'.png')
    
    
    # Marginalised distribution (histogram) plot.
    pl.hist(sampler.flatchain[:,0], 100)
    pl.show()
    
    
    # Plotting lines of best fit using a 100-strong sample of parameters.
    xl = np.linspace(0,4,100)
    #for a, b, c, d, e in samples[np.random.randint(len(samples), size=100)]:
     #   pl.plot(xl, a * xl**4 + b * xl**2 + c * xl + d +
      #          e*np.sin(xl), color="k", alpha=0.1)
    pl.plot(xl, a_true * xl**4 + b_true * xl**2 + c_true * xl + d_true + 
            e_true*np.sin(xl),color="r", lw=2, alpha=0.8)
    pl.errorbar(x, y, yerr=yerr, fmt=".k")
    pl.show()
    
    
    # Best line of fit found by emcee.
    bi = np.argmax(sampler.lnprobability)   # index with highest post prob                                       
    abest = sampler.flatchain[bi,0]         # parameters with the highest 
    bbest = sampler.flatchain[bi,1]         # posterior probability
    cbest = sampler.flatchain[bi,2]
    dbest = sampler.flatchain[bi,3]
    ebest = sampler.flatchain[bi,4]
    
    
    # plot of data with errorbars + model
    pl.errorbar(x, y, yerr=sigma, fmt='o', alpha=0.3)
    xt = np.linspace(0,4,100)
    yt = (a_true * xt**4 + b_true * xt**2 + c_true * xt + d_true 
          + e_true * np.sin(xt))
    model, = pl.plot(xt,yt,lw='3', c='g')
    ybest = (abest * xt**4 + bbest * xt**2 + cbest * xt + dbest 
             + ebest * np.sin(xt))
    best_fit, = pl.plot(xt,ybest,lw='3', c='r')
    pl.legend([model, best_fit], ['Model', 'Best Fit'])
    pl.show
    
    
    timet1=time.time()      # stopping script time
    timet=timet1-timet0     # total time to run script
    timetmin = round((timet / 60),1)  # minutes
    timetsec = round((timet % 60),1)  # seconds
    
    
    # Results getting printed:
    print('best index is =',str(bi))
    print('abest is =',str(abest))
    print('bbest is =',str(bbest))
    print('cbest is =',str(cbest))
    print('dbest is =',str(dbest))
    print('ebest is =',str(ebest))
    # Mean acceptance fraction. In general, acceptance fraction has an entry 
    # for each walker so, in this case, it is a 50-dimensional vector.
    print('Mean acceptance fraction:', np.mean(sampler.acceptance_fraction))
    print('Number of steps:', str(nsteps))
    print('Number of walkers:', str(nwalkers))
    print('Sampler time:',str(int(timesmin))+'min'
          ,str(int(timessec))+'s')
    print('Total time:  ',str(int(timetmin))+'min'
          ,str(int(timetsec))+'s')
    
    
except Exception as e:
        logging.error('Caught exception:',str(e))
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno))






