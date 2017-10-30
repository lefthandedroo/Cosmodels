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

what interaction terms could there be?

git test

NOT AN ASSIGNMENT
"""

from scipy.integrate import odeint
import numpy as np
from pylab import figure, plot, xlabel, grid, legend, title
from matplotlib.font_manager import FontProperties

def vectorfield(v, t, w):
    a, a_dot, e_dash_m, e_dash_de, omega0m, omega0de = v
    w_m, w_de = w
    lamb_m = -2e-2 #10
    lamb_de = 2e-2
    # Create f = [a_dot, a_dotdot, e'_dotm, e'_dotde, omegam_dot, omegade_dot]:
    f = [a_dot, 
         (-a/2)*(e_dash_m*(1+3*w_m)+e_dash_de*(1+3*w_de)), 
         -3*(a_dot/a)*(e_dash_m*(1+w_m+lamb_m/3*a/a_dot)),
         -3*(a_dot/a)*(e_dash_de*(1+w_de+lamb_de/3*a/a_dot)),
         (H0**2)*e_dash_m*(2*a/a_dot)*(1-a*(-a/2)*(e_dash_m*(1+3*w_m)+
                           e_dash_de*(1+3*w_de))/a_dot**2)+(a/a_dot)**2*
                           (-3*a_dot/a*e_dash_m*(1+w_m)),
         (H0**2)*e_dash_de*(2*a/a_dot)*(1-a*(-a/2)*(e_dash_m*(1+3*w_m)+
                            e_dash_de*(1+3*w_de))/a_dot**2)+(a/a_dot)**2*
                            (-3*a_dot/a*e_dash_de*(1+w_de))]
        
    return f

# a past which to discard values, value chosen by looking at the plot 
# set arbitrarily - sometimes jumps over the result(?)
a_d = 10e-3 

# time in 1/H0 to integrate until, if this time isn't long enough to reach a_d
# then "time" will be added to integration time until a_d is reached
# 0.665 matter only, 0.96 m+de, -0.49
time = 0.6

# Parameters
H0 = 1       # Hubble parameter at t=now
#Dh = c/H0   # Hubble distance
tH = 1.0/H0  # Hubble time
# G = 1
# c = 1

# Eq of state parameter
w_m = 0.0     # matter
w_de = -1.0   # cosmological constant (dark energy)
w_r = 1/3     # radiation

# Initial conditions
# a0 = scale factor, a_dot = speed, e_dash0 = e0/ec0
a0 = 1.0
a_dot0 = 1.0
e_dash0m = 0.3
e_dash0de = 0.7
omega0m = 0.3
omega0de = 0.7

# ODE solver parameters
abserr = 1.0e-8
relerr = 1.0e-6
numpoints = 250

stoptime = -time

# Plot the solution
figure()
xlabel('t in 1/H0')
grid(True)
lw = 1

while True:
#    print('stoptime is:',str(stoptime),'1/H0')
    # Create the time samples for the output of the ODE solver.
    t = [stoptime*tH * float(i) / (numpoints - 1) for i in range(numpoints)]
    
    # Pack up the parameters and initial conditions:
    v0 = [a0, a_dot0, e_dash0m, e_dash0de, omega0m, omega0de]
    
    w = [w_m, w_de]
    
    # Call the ODE solver.
    vsol = odeint(vectorfield, v0, t, args=(w,), atol=abserr, rtol=relerr)
    
    # Remove unwanted results from the plot
    # Separate results into their own arrays
    a = vsol[:,0]
    a_dot = vsol[:,1]
    e_dashm = vsol[:,2]
    e_dashde = vsol[:,3]
    omegam = vsol[:,4]
    omegade = vsol[:,5]
    
    # Find where results start to get strange
    blowups = np.where(a < a_d)    # tuple with indecies of a
                                   # so small it blows up a_dot
    blowups = np.asarray(blowups)  # converting to np array   

    if blowups.any():                       
        blowup = blowups[0,0]      # first instance of a being too small
    else: 
        stoptime -= time
        continue
    
    # Remove the values after the index when a is too small
    t_cut = np.asarray(t)
    t_cut = t_cut[:blowup]
    a = a[:blowup]
    a_dot = a_dot[:blowup]
    e_dashm = e_dashm[:blowup]
    e_dashde = e_dashde[:blowup]
    omegam = omegam[:blowup]
    omegade = omegade[:blowup]
    
    # Age
    age = t_cut[np.argmin(t_cut)]
    age = round(age, 2)
    
    # plotting selected resutls
    # a and a_dot
    xlabel('t in 1/H0')
    lw = 1
    plot(t_cut, a, 'r', linewidth=lw)
    plot(t_cut, a_dot, 'b', linewidth=lw)
    legend((r'$a$', r'$\.a$'), prop=FontProperties(size=16))
    title('Cut results for $\omega$ = %s, age = %s 1/H0'%(w,age))
    
    while False:
        # e_dashm
        figure()
        xlabel('t in 1/H0')
        lw = 1
        plot(t_cut, e_dashm, 'g', linewidth=lw)
        title('Cut $\epsilon_m \'$ for $\omega$ = %s'%(w))
        
        # e_dashde
        figure()
        xlabel('t in 1/H0')
        lw = 1
        plot(t_cut, e_dashde, 'm', linewidth=lw)
        title('Cut $\epsilon_{DE} \'$ for $\omega$ = %s'%(w))
        
        break
    
    # omegam and omegade
    figure()
    plot(t_cut, omegam, 'c', linewidth=lw)
    plot(t_cut, omegade, 'k', linewidth=lw)
    legend((r'$\Omega_m$', r'$\Omega_{DE}$'), prop=FontProperties(size=16))
    title('Cut results for $\omega$ = %s, age = %s 1/H0'%(w,age))
    
    break
        


while False:  
    figure()
    xlabel('t in 1/H0')
    grid(True)
    lw = 1
    
    # plotting all results
    plot(t, vsol[:,0], 'r', linewidth=lw)
    plot(t, vsol[:,1], 'b', linewidth=lw)
    
    legend((r'$a$', r'$\.a$', r'$\'\epsilon$'), prop=FontProperties(size=16))
    title('Complete results for $\omega$ = %s'%(w))

    break