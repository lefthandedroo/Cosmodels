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
put in mode than one fluid in the universe

matter + de universe (age)
"""

from scipy.integrate import odeint


def vectorfield(v, t, w):
    a, a_dot, e_dash = v

    # Create f = [a_dot, a_dotdot, e'_dot]:
    f = [a_dot, (-a/2)*e_dash*(1+3*w), -3*(a_dot/a)*e_dash*(1+w)]
    return f


# Parameters
H0 = 1       # Hubble parameter at t=now
#Dh = c/H0   # Hubble distance
tH = 1.0/H0  # Hubble time

# Eq of state parameter
w_m = 0.0     # matter
w_de = -1.0   # cosmological constant (dark energy)
w_r = 1/3     # radiation
w = w_m #+ w_de + w_r

# Initial conditions
# a0 = scale factor, a_dot = speed, e_dash0 = e0/ec0
a0 = 1.0
a_dot0 = 1.0
e_dash0 = 1.0

# ODE solver parameters
abserr = 1.0e-8
relerr = 1.0e-6
numpoints = 250

# Create the time samples for the output of the ODE solver.
stoptime = -0.49
t = [stoptime*tH * float(i) / (numpoints - 1) for i in range(numpoints)]

# Pack up the parameters and initial conditions:
v0 = [a0, a_dot0, e_dash0]

# Call the ODE solver.
vsol = odeint(vectorfield, v0, t, args=(w,), atol=abserr, rtol=relerr)


# Plot the solution
from pylab import figure, plot, xlabel, grid, legend, title
from matplotlib.font_manager import FontProperties


figure()

xlabel('t in 1/H0')
grid(True)
lw = 1

plot(t, vsol[:,0], 'r', linewidth=lw)
plot(t, vsol[:,1], 'b', linewidth=lw)
#plot(t, vsol[:,2], 'g', linewidth=lw)

legend((r'$a$', r'$\.a$', r'$\'\epsilon$'), prop=FontProperties(size=16))
title('Matter only')