#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Coupled spring-mass system
"""

def vectorfield(w, t, p):
    """
    Defines the differential equations for the coupled 
    spring-mass system.

    Arguments:
        w :  vector of the state variables:
                  w = [x1,y1,x2,y2]
        t :  time
        p :  vector of the parameters:
                  p = [m1,m2,k1,k2,L1,L2,b1,b2]
    """
    x1, y1, x2, y2 = w
    m1, m2, k1, k2, L1, L2, b1, b2 = p

    # Create f = (x1',y1',x2',y2'):
    f = [y1,(-b1 * y1 - k1 * (x1 - L1) + k2 * (x2 - x1 - L2)) / m1,y2,(-b2 * y2 - k2 * (x2 - x1 - L2)) / m2]
    return f

# Use ODEINT to solve the differential equations defined 
# by the vector field
from scipy.integrate import odeint

# Parameter values
# Masses:
m1 = 1.0
m2 = 1.0
# Spring constants
k1 = 8.0
k2 = 8.0
# Natural lengths
L1 = 1.5
L2 = 1.5
# Friction coefficients
b1 = 0.8
b2 = 0.8

# Initial conditions
# x1 and x2 are the initial displacements; y1 and y2 are 
# the initial velocities
x1 = 0.0
y1 = 0.0
x2 = 3.5
y2 = 0.0

# ODE solver parameters
abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 10.0
numpoints = 250

# Create the time samples for the output of the ODE solver.
# I use a large number of points, only because I want to make
# a plot of the solution that looks nice.
t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]

# Pack up the parameters and initial conditions:
p = [m1, m2, k1, k2, L1, L2, b1, b2]
w0 = [x1, y1, x2, y2]

# Call the ODE solver.
wsol = odeint(vectorfield, w0, t, args=(p,),
              atol=abserr, rtol=relerr)
        
# Plot the solution
#from numpy import loadtxt
from pylab import figure, plot, xlabel, grid, legend, title, savefig
from matplotlib.font_manager import FontProperties


figure(1, figsize=(6, 4.5))

xlabel('t')
grid(True)
lw = 1

plot(t, wsol[:,0], 'y', linewidth=lw)
plot(t, wsol[:,2], 'g', linewidth=lw)
plot(t, wsol[:,1], 'r', linewidth=lw)
plot(t, wsol[:,3], 'b', linewidth=lw)

legend((r'$x_1$', r'$x_2$', r'$v_1$', r'$v_2$'), prop=FontProperties(size=16))
title('Mass Displacements for the\nCoupled Spring-Mass System')
savefig('two_springs.png', dpi=100)