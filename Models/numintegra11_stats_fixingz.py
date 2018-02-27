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
check p use vs sigma, why are errobars on the Model with errorbars plot so small?
call the interactin term gamma everywhere, are any lambda - cosm constant? no
make numintegra8 give D_L in parsecs



3) What would happen if error were 1%? or 50%? How does it change my parameters
4) What if I generate a cosmology that has an interaction term? 
5) Choose another cosmology, generate fake data and try the same with omegam,
    omegade and lamb
(try the standard one with no lambda and see if you get lambda = 0 back)
6) Working towards distributions of omega lambda, omega matter and interaction term.


Put constraint that m+de can't be more than 1, ask Geraint if that's reasonable 

separate plotting function
Save outcomes of variables and plot everything vs everything after script finishes
Ask Sue Yang re debugging python to check type of error, underflow/overflow
Run odesolve with slightly different parameters, plot redshift vs dlmpc and see if 
anything looks super strange (some parameters might be unphysical)


NOT AN ASSIGNMENT
"""
import numpy as np
import time

import gnoise
import msim
import zpicks
import stats

# Starting script timer.
timet0 = time.time()


# Parameters:

# Model specific parameters.  
gamma_true = 0  # Interaction term, rate at which DE decays into matter.
m_true = 0.3    # (= e_m(t)/e_crit(t0) at t=t0).
de_true = 0.7   # (de = e_de(t)/e_crit(t0) at t=t0).

# Number of datapoints to be simulated.
n = 100 #10, 1000

# Statistical parameters:
mu = 0          # mean
sigma = 1       # standard deviation



# Code:

# Picking redshifts to investigate.
zmin, zmax = 0.001, 2
zpicks = zpicks.zpicks(zmin, zmax, n)

# Generating apparent magnitues mag at redshift z < zmax (calculated from
# luminosity distances given by LambdaCMD with parameters stated above.
model = msim.msim(gamma_true, m_true, de_true, n, zpicks)
model = np.asarray(model)
mag, noise = gnoise.gnoise(model, mu, sigma, n)
#print('noise in code body is = ', noise)


stats.stats(gamma_true, m_true, de_true, n, zpicks, mag, noise, sigma)


# Time taken by the script. 
timet1=time.time()      # stopping script time
timet=timet1-timet0     # total time to run script
timetmin = round((timet / 60),1)  # minutes
timetsec = round((timet % 60),1)  # seconds
print('Total time:  ',str(int(timetmin))+'min',str(int(timetsec))+'s')