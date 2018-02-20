#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
use 2nd order friedmann equations + the continuity equation
numerically integrate the cosmological equations for standard cosmology where
omega = 1 for matter universe
also
omega = 1 30% matter, 70% dark energy
luminosity distance plots
https://arxiv.org/pdf/astro-ph/9905116.pdf
pl.plot(sampler.flatchain[:,0])
pl.plot(sampler.flatchain[:,1])
pl.plot(sampler.flatchain[:,0],sampler.flatchain[:,1],'.')
interaction term next
>>ODE solver (ODE int)
integrate both scale factor and continuity equation at the same time
change units
use time = 1/h0 and distance c/h0
AND
integrate backwards from today so from 0 ti -t
ok to ask Luke
later : movify to include co-moving distance (Ryden?)
+ Luminosity distance
find out how old a universe with only matter and at 
critical density would be today 
integrate over time, t0 will be now, where a crosses 0 is the beginning
so find time in years
then tell how a changes with time and the age of the universe in universe with 
critical 
density today, single fluid an with equation of state w=-1/3
Dodn't always ask what to do
event catcher - stop integration once a reaches 0
put in more than one fluid in the universe
maybe http://www.ni.gsu.edu/~rclewley/PyDSTool/FrontPage.html
m + DE universe (age)
Create plot of omega lambda vs omaga matter across age of universe 
0.3 and 0.7 today, what were they in the past 
Krysztoff, Luke
[try a_dot as equation and not just a_dot, i.e. ditch a_dotdot, but that is if 
the equations are correct. if solver doens't have a problem when feeding data 
in at each calculation then a_dot as a_dot and not equation shouldn't matter] -
worked out in the end
git test
what interaction terms could there be?
add lambda to Omega calculation 
add luminosity distance to integration and plot
move input to top (z, omega)
find absolute magniteds for supernova (-19)
Classic paper on Gaussian processes from Ben's lecture
1)Generate a sample of supernova from omegam=0.3 omegade=0.7, lamb=0 model:
    a) convert luminosity distances into parsecs using 
        D_L = what I found * c/H0 = what I found * 4167 Mpc
    b) Pick a random distrib of redshifts up to z=2: try 10, 100, 1000.
    c) Plug the corresponding luminosity distances (=r in parcecs) and M (=-19) 
        into distance modulus (mag - M = 5 * log10(r/10) ) 
        to find their apparent magnitudes mag.
    d) Add 10% gaussian noise to apparent magnitudes mag.
Write pseudocode first to get the understanding about what I'm doing.
lookup how to call a function from a function from a function to keep them separate
but be able to call in one line
neaten up the code   
2) Use mcmc to re-fit the data to that sample and work out what parameters 
    (constraints) I get on omegam and omegade
generate 10^6 points for t (up to z=2)
pick random points between 0<z<2
interpolate to find "exact" m for specific z I randomly picked
need to be happy that the itnerpolation gives an accurate answer 
(numerical recipy's  book has interpaoplation)
https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.interp.html
output vsol and check type
10 point sample does not generate enough points
SPLIT INTO MODULES

3) What would happen if error were 1%? or 50%? How does it change my parameters
4) What if I generate a cosmology that has an interaction term? 
5) Choose another cosmology, generate fake data and try the same with omegam,
    omegade and lamb
(try the standard one with no lambda and see if you get lambda = 0 back)
6) Working towards distributions of omega lambda, omega matter and interaction term.
make numintegra8 give D_L in parsecs
put plots and calculation into separate functions


check p use vs sigma, why are errobars on the Model with errorbars plot so small?

Make plot of everything vs everything through saving outcomes of variables and plot them after
change interaction term to gamma instead of labda, mindful of lambda teh cosm constant
Ask Sue Yang re debugging python to check type of error, underflow/overflow
Run msim with slightly different parameters, plot redshift vs dlmpc and see if 
anything looks super strange (some parameters might be unphysical)


NOT AN ASSIGNMENT
"""
# cosmo
import numpy as np
import random
from pylab import figure, scatter

import corner
import emcee
import logging
import matplotlib.pyplot as pl
import scipy.optimize as op
import sys
import time

import flist
import gnoise
import msim
import lnlike
import lnprob

timet0 = time.time()    # starting script timer

# Model specific parameters.  
# Interaction term, rate at which DE decays into matter.
lamb_true = 0
# Fraction of matter and dark energy compared to critical density at t=t0.
m_true = 0.3
de_true = 0.7

# Number of datapoints to be simulated.
n = 100 #10, 1000

# Statistical parameters:
# Percentage error on apparent magnitudes.
p = 10
sigma = p/100       # standard deviation
mu = 0              # mean

# emcee parameters:
ndim, nwalkers = 3, 6
nsteps = 1000
burnin = 200

 



# Code:

# Picking redshifts to investigate.
zmin = 0.001
zmax = 2        # Largest meausured z for a supernovae is 2.
zinterval = (zmax - zmin) / (n*2)
z_opts = flist.flist(zmin, zmax, zinterval)
zpicks = random.sample(z_opts, n)
zpicks = np.asarray(zpicks)

# Generating apparent magnitues mag at redshift z<2 (calculated from
# luminosity distances given by LambdaCMD with parameters stated above.
#theta = lamb, m, de
model = msim.msim(lamb_true, m_true, de_true, n, p, zpicks)
model = np.asarray(model)

mag, noise = gnoise.gnoise(model, mu, sigma, n)
#print('noise in code body is = ', noise)



try:
    # Finding a "good" place to start using alternative method to emcee.
    nll = lambda *args: -lnlike.lnlike(*args)  # type of nll is: <class 'function'>
    result = op.minimize(nll, [lamb_true, m_true, de_true], 
                         args=(n, p, zpicks, mag, noise))
    lamb_ml, m_ml, de_ml = result["x"]    
    
        
    # Initializing walkers in a Gaussian ball around the max likelihood. 
    pos = [result["x"] + 1*np.random.randn(ndim) for i in range(nwalkers)]    
        
    
    # Sampler setup
    times0 = time.time()    # starting emcee timer
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob.lnprob, args=(n, p, zpicks, mag, sigma))
    sampler.run_mcmc(pos, nsteps)
    
    times1=time.time()      # stopping emcee timer
    times=times1 - times0   # time to run emcee
    timesmin = round((times / 60),1)    # minutes
    timessec = round((times % 60),1)    # seconds
    
    
    # Corner plot (walkers' walk + histogram).
    samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
    fig = corner.corner(samples, labels=["$lamb$", "$m$", "$de$"], 
                        truths=[lamb_true, m_true, de_true])
    fig.savefig('nsteps'+str(nsteps)+str(time.strftime("%c"))+
                'nwalkers'+str(nwalkers)+'.png')
    
    
    # Marginalised distribution (histogram) plot.
    pl.hist(sampler.flatchain[:,0], 100)
    pl.show()
    
    # Plotting lines of best fit using a 100-strong sample of parameters.
    figure()
    pl.title('Model with errorbars')
    scatter(zpicks, model, color="r", lw=2, alpha=0.8)
    pl.errorbar(zpicks, model, yerr=sigma, fmt=".k")
    pl.show()
    
    # Best line of fit found by emcee.
    bi = np.argmax(sampler.lnprobability)   # index with highest post prob                                       
    lambbest = sampler.flatchain[bi,0]      # parameters with the highest 
    mbest = sampler.flatchain[bi,1]         # posterior probability
    debest = sampler.flatchain[bi,2]
    
    # plot of data with errorbars + model
    figure()
    pl.title('Model and Best Fit')
    pl.errorbar(zpicks, mag, yerr=sigma, fmt='o', alpha=0.3)
    modelt = msim.msim(lamb_true, m_true, de_true, n, p, zpicks)
    model_fit = scatter(zpicks, modelt, lw='3', c='g')
    magbest = msim.msim(lambbest, mbest, debest, n, p, zpicks)
    best_fit = scatter(zpicks,magbest,lw='3', c='r')
    pl.legend([model_fit, best_fit], ['Model', 'Best Fit'])
    pl.show()
    
    
    # Results getting printed:
    print('best index is =',str(bi))
    print('lambbest is =',str(lambbest))
    print('mbest is =',str(mbest))
    print('debest is =',str(debest))
  
    # Mean acceptance fraction. In general, acceptance fraction has an entry 
    # for each walker so, in this case, it is a 50-dimensional vector.
    print('Mean acceptance fraction:', np.mean(sampler.acceptance_fraction))
    print('Number of steps:', str(nsteps))
    print('Number of walkers:', str(nwalkers))
    print('Sampler time:',str(int(timesmin))+'min'
          ,str(int(timessec))+'s')
    
    
except Exception as e:
        logging.error('Caught exception:',str(e))
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno))


timet1=time.time()      # stopping script time
timet=timet1-timet0     # total time to run script
timetmin = round((timet / 60),1)  # minutes
timetsec = round((timet % 60),1)  # seconds
print('Total time:  ',str(int(timetmin))+'min',str(int(timetsec))+'s')