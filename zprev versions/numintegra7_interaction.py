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





NOT AN ASSIGNMENT
"""

from scipy.integrate import odeint
import numpy as np
from pylab import figure, plot, scatter, xlabel, grid, legend, title, annotate
from matplotlib.font_manager import FontProperties

#def vectorfield(v, t, w):
#    """
#    Takes in    v = values at time now
#                w = omega parameters
#                lamb = interaction constant.
#                
#    Returns a function with     a_dot, a_dotdot, 
#                                e'_dotm, e'_dotde, 
#                                omegam_dot, omegade_dot,
#    ready to be integrated with odeint.
#    Uses different lambdas for different fluids.
#    """
#    (a, a_dot, e_dash_m, e_dash_de, omega0m, omega0de) = v
#    (w_m, w_de) = w
#    
#    lamb_m = -2e-2 #10
#    lamb_de = 2e-2
#    
#    # Create f = [a_dot, a_dotdot, e'_dotm, e'_dotde, omegam_dot, omegade_dot]:
#    f = [# a_dot
#         a_dot,
#         # a_dotdot
#         (-a/2) * (e_dash_m * (1+3*w_m) + e_dash_de * (1+3*w_de)),
#         # e'_dotm
#         -3 * (a_dot/a) * (e_dash_m * (1 + w_m + lamb_m/3 * a/a_dot)),
#         # e'_dotde
#         -3 * (a_dot/a) * (e_dash_de * (1 + w_de + lamb_de/3 * a/a_dot)),
#         # omegam_dot
#         (H0**2) * e_dash_m * (2*a/a_dot) * (1 - a * (-a/2) 
#             * (e_dash_m * (1+3*w_m) + e_dash_de * (1+3*w_de)) 
#             / a_dot**2) + (a/a_dot)**2 * (-3 * a_dot/a * e_dash_m * (1+w_m)),
#        # omegade_dot
#         (H0**2) * e_dash_de * (2*a/a_dot) * (1 - a * (-a/2) 
#             * (e_dash_m * (1+3*w_m) + e_dash_de * (1+3*w_de)) 
#             / a_dot**2) + (a/a_dot)**2 * (-3 * a_dot/a * e_dash_de * (1+w_de))
#         ]
#        
#    return f

def vectorfield(v, t, w, lamb):
    """
    Takes in    v = values at time now
                w = omega parameters
                lamb = interaction constant.
                
    Returns a function with     a_dot, a_dotdot, 
                                e'_dotm, e'_dotde, 
                                omegam_dot, omegade_dot,
    ready to be integrated with odeint.
    Uses same lambdas for all fluids.
    """
    (a, a_dot, e_dashm, e_dashde, omegam, omegade) = v
    (w_m, w_de) = w
    
    # Create f = [a_dot, a_dotdot, e'_dotm, e'_dotde, omegam_dot, omegade_dot]:
    f = [# a_dot
         a_dot,
         # a_dotdot
         (-a/2) * (e_dashm * (1+3*w_m) + e_dashde * (1+3*w_de)), 
         # e'_dotm
         -3 * (a_dot/a) * (e_dashm * (1 + w_m -lamb/3 * a/a_dot 
              * e_dashde/e_dashm )),
         # e'_dotde
         -3 * (a_dot/a) * (e_dashde * (1 + w_de +lamb/3 * a/a_dot)),
         # omegam_dot
         (H0**2) * e_dashm * (2*a/a_dot) * (1 - a * (-a/2) 
             * (e_dashm * (1+3*w_m) + e_dashde * (1+3*w_de)) 
             / a_dot**2) + (a/a_dot)**2 * (-3 * a_dot/a * e_dashm * (1+w_m)),
         # omegade_dot        
         (H0**2) * e_dashde * (2*a/a_dot) * (1 - a * (-a/2) 
             * (e_dashm * (1+3*w_m) + e_dashde * (1+3*w_de)) 
             / a_dot**2) + (a/a_dot)**2 * (-3 * a_dot/a * e_dashde * (1+w_de))
         ]
        
    return f


# Parameters (script specific)
    
# Interaction term, rate at which DE decays into matter.
lamb = 0.25

# Last value for a before results are to be considered 
# invalid due to close proximity to big bang, a_d is set
# somewhat arbitrarily - sometimes jumps over the result(?).
a_d = 10e-6 

# Time (in 1/H0) to integrate until.  If this time isn't long enough for a to 
# decrease to a_d then stoptime will be extended by time until a_d is reached.
# 0.665 matter only, 0.96 standard m+de
time = 0.6

# Standard cosmological parameters.
H0 = 1       # Hubble parameter at t=now
tH = 1.0/H0  # Hubble time
# Eq of state for known fluids:
w_r = 1/3     # radiation
w_m = 0.0     # matter
w_de = -1.0   # cosmological constant (dark energy?)

# Initial conditions at time = t0.
a0 = 1.0        # scale factor
a_dot0 = 1.0    # speed of expansion
e_dash0m = 0.3  # e_m(t)/ec(t0)
e_dash0de = 0.7 # e_de(t)/ec(t0)
omega0m = 0.3   # e_m(t)/ec(t)
omega0de = 0.7  # e_de(t)/ec(t)

# ODE solver parameters:
abserr = 1.0e-8
relerr = 1.0e-6
numpoints = 250

stoptime = -time # Integrating back in time as time now is t0.

while True:
    # Create time samples for the ODE solver.
    t = [stoptime * tH * float(i) / (numpoints - 1) for i in range(numpoints)]
    
    # Pack up the initial conditions and eq of state parameters.
    v0 = [a0, a_dot0, e_dash0m, e_dash0de, omega0m, omega0de]
    w = [w_m, w_de]
    
    # Call the ODE solver.
    vsol = odeint(vectorfield, v0, t, args=(w,lamb,), atol=abserr, rtol=relerr)
    
    # Remove unwanted results which are too close to big bang from the plot.
    # Separate results into their own arrays:
    a = vsol[:,0]
    a_dot = vsol[:,1]
    e_dashm = vsol[:,2]
    e_dashde = vsol[:,3]
    omegam = vsol[:,4]
    omegade = vsol[:,5]

    
    # Find where results start to get strange (smaller than a_d):
    blowups = np.where(a < a_d)    # Tuple with indecies of a so
                                   # small that other results blow up.                             
    blowups = np.asarray(blowups)  # Converting to np array.

    if blowups.any():              # Check if instances of a < a_d exist.   
        blowup = blowups[0,0]
    else:                          # If no instance of a < a_d was found
        stoptime -= time           # then integrate further back in time.
        continue
    
    
    # Remove the values after the index of first instance of a < a_d.
    t_cut = np.asarray(t)
    t_cut = t_cut[:blowup]
    a = a[:blowup]
    a_dot = a_dot[:blowup]
    e_dashm = e_dashm[:blowup]
    e_dashde = e_dashde[:blowup]
    omegam = omegam[:blowup]
    omegade = omegade[:blowup]
    omegatot = np.add(omegam,omegade)
    
    # Find when density of DE was equal to that of matter.  Rounding the omegas
    # to avoid not finding an instance of equality.
    equiv = np.argmin(abs(omegam - omegade))    # Index of time of equivlence.
    equiv = np.asarray(equiv)                   # Converting to nparray.
    
    # Age of the universe.
    age = t_cut[np.argmin(t_cut)]
    age = -round(age, 2)
    
    # Plotting selected results:
    # a and a_dot
    figure()
    xlabel('time in $H_0^{-1}$')
    grid(True)
    plot(t_cut, a, 'r', t_cut, a_dot, 'b', lw=1)
    legend((r'$a$', r'$\.a$'), prop=FontProperties(size=16))
    title('Cut results for $\omega$ = %s, $\lambda$ = %s, age = %s $H_0^{-1}$'
          %(w,lamb,age))
    
    while False:    # Looped to make it faster to switch plots on and off.
        # e_dashm
        figure()
        xlabel('time in $H_0^{-1}$')
        plot(t_cut, e_dashm, 'g', lw=1)
        title('Cut $\epsilon_m \'$ for $\omega$ = %s'%(w))
        
        # e_dashde
        figure()
        xlabel('time in $H_0^{-1}$')
        plot(t_cut, e_dashde, 'm', lw=1)
        title('Cut $\epsilon_{DE} \'$ for $\omega$ = %s'%(w))
        
        break
    
    # omegam and omegade calculated using the integrator:
    figure()
    xlabel('time in $H_0^{-1}$')
    grid(True)
    plot(t_cut, omegam, 'c', t_cut, omegade, 'k', lw=1)
    plot(t_cut, omegatot, 'tab:orange', ls= 'dashed', lw=1)
    scatter(t_cut[equiv], omegade[equiv], s=80, 
            facecolors='none', edgecolors='r')
    legend((r'$\Omega_m$', r'$\Omega_{DE}$', r'$\Omega_{(m + DE)}$', 
            r'$t_{m = DE}$'), prop=FontProperties(size=16))
    annotate('$t_{m = DE}$ = %s $H_0^{-1}$'%(round(t_cut[equiv], 4)), 
            xy=(t_cut[equiv], 0.45), 
            xytext=(-0.35, np.min(omegade)+0.03), 
            arrowprops = dict(facecolor='white', shrink=0.05))
    title('Cut results for $\omega$ = %s, $\lambda$ = %s, age = %s $H_0^{-1}$'
          %(w,lamb,age))
    
    # Verified omegam and omegade using e_dashm and e_dashde.
    vomegam = e_dashm / (e_dashm + e_dashde)
    vomegade = e_dashde / (e_dashm + e_dashde)
    
    figure()
    xlabel('time in $H_0^{-1}$')
    grid(True)
    plot(t_cut,vomegam,t_cut,vomegade)
    legend((r'$\Omega_m$', r'$\Omega_{DE}$'), prop=FontProperties(size=16))
    title('Verified $\Omega_m$ and $\Omega_{DE}$ for $\lambda$ = %s, age = %s $H_0^{-1}$'
          %(lamb,age))
    break
        

# Complete results with blow up resulting from a approaching big bang point.
while False:  
    figure()
    xlabel('time in $H_0^{-1}$')
    grid(True)
    
    # Plotting complete results.
    plot(t, vsol[:,0], 'r', lw=1)
    plot(t, vsol[:,1], 'b', lw=1)
    
    legend((r'$a$', r'$\.a$', r'$\'\epsilon$'), prop=FontProperties(size=16))
    title('Complete results for $\omega$ = %s'%(w))

    break