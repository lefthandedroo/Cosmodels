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
"""

from scipy.integrate import odeint


def vectorfield(v, t, w):
    a, a_dot, e_dash = v
#    w = p

    # Create f = (a',a_dot',e_dash'):
    f = [a_dot, -a / 2 * e_dash * (1 + 3 * w), - 3 * a_dot / a * e_dash * (1+ w)]
    return f

# Parameters
#O_0 = 1.0     # Omega for flat universe
#O_m = 0.3
#O_r = 0.0
#O_de = 0.7
#k = 0       # curvature of flat universe
#G = 1       # Newton's gravitatinal constant
#c = 1       # speed of light
H0 = 1      # Hubble parameter at t=now
#Dh = c/H0   # Hubble distance
th = 1.0/H0   # Hubble time

# Eq of state parameter
w_m = 0.0     # matter
w_de = -1.0   # cosmological constant (dark energy)
w_r = 1/3   # radiation
w = w_m + w_de + w_r
#print('w is: ',w)

## pressure
#p_m = omega_m * rho_m     # matter
#p_de = omega_de * rho_de  # DE
#p = p_m + p_de

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
# I use a large number of points, only because I want to make
# a plot of the solution that looks nice.
t = [-2*th * float(i) / (numpoints - 1) for i in range(numpoints)]
#t = range((-2*th), 0, 1)

# Pack up the parameters and initial conditions:
#p = [w]
v0 = [a0, a_dot0, e_dash0]

#print('up to line 97')

# Call the ODE solver.
vsol = odeint(vectorfield, v0, t, args=(w,),
              atol=abserr, rtol=relerr)


# Plot the solution
#from numpy import loadtxt
from pylab import figure, plot, xlabel, grid, legend, title
from matplotlib.font_manager import FontProperties


figure()

xlabel('t in 1/H0')
grid(True)
lw = 1

plot(t, vsol[:,0], 'r', linewidth=lw)
plot(t, vsol[:,1], 'b', linewidth=lw)
#plot(t, vsol[:,2], 'r', linewidth=lw)

legend((r'$\.a$', r'$\"a$', r'$edashdot$'), prop=FontProperties(size=16))
title('Matter only')