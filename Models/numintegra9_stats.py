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
   




2) Use mcmc to re-fit the data to that sample and work out what parameters 
    (constraints) I get on omegam and omegade

3) What would happen if error were 1%? or 50%? How does it change my parameters

4) What if I generate a cosmology that has an interaction term? 

5) Choose another cosmology, generate fake data and try the same with omegam,
    omegade and lamb
(try the standard one with no lambda and see if you get lambda = 0 back)

6) Working towards distributions of omega lambda, omega matter and interaction term.


make numintegra8 give D_L in parsecs
put plots and calculation into separate functions


lookup how to call a function from a function froa a function to keep them separate
but be able to call in one line
neaten up the code


NOT AN ASSIGNMENT
"""
# cosmo
import numpy as np
import random
from scipy.integrate import odeint
from math import log10
from pylab import figure, plot, scatter, xlabel, grid, legend, title, annotate
from matplotlib.font_manager import FontProperties
import corner
import emcee
import logging
import matplotlib.pyplot as pl
import scipy.optimize as op
import sys
import time

def firstderivs(v, t, w, lamb):
    """
    Takes in:
        v = values at t=10;
        w = omega parameters;
        lamb = interaction constant.
                
    Returns a function with:    a_dot, a_dotdot, 
                                e'_dotm, e'_dotde, 
                                omegam_dot, omegade_dot,
                                z_dot,
                                dl_dot
    ready to be integrated with odeint.
    Uses same lambda for all fluids.
    """
    (a, a_dot, e_dashm, e_dashde, omegam, omegade, z, dl) = v
    (w_m, w_de) = w
    
    # fist derivatives of functions I want to find:
    f = [# a_dot (=scale factor)
         a_dot,
         # a_dotdot
         (-a/2) * (e_dashm * (1+3*w_m) + e_dashde * (1+3*w_de)), 
         # e'_dotm (=density(t) / crit density(t0))
         -3 * (a_dot/a) * e_dashm * (1 + w_m -lamb/3 * a/a_dot 
              * e_dashde/e_dashm ),
         # e'_dotde
         -3 * (a_dot/a) * e_dashde * (1 + w_de +lamb/3 * a/a_dot),
         # omegam_dot (=density(t) / crit density(t))
         (H0**2) * e_dashm * (2*a/a_dot) * (1 - a * (-a/2) 
             * (e_dashm * (1+3*w_m) + e_dashde * (1+3*w_de)) 
             / a_dot**2) + (a/a_dot)**2 * (-3 * (a_dot/a) * e_dashm 
             * (1 + w_m -lamb/3 * a/a_dot * e_dashde/e_dashm)),
         # omegade_dot        
         (H0**2) * e_dashde * (2*a/a_dot) * (1 - a * (-a/2) 
             * (e_dashm * (1+3*w_m) + e_dashde * (1+3*w_de)) 
             / a_dot**2) + (a/a_dot)**2 * (-3 * (a_dot/a) 
             * e_dashde * (1 + w_de +lamb/3 * a/a_dot)),
         # z_dot (=redshift)
         -a_dot/a**2,
         # dl_dot (=luminosty distance)
         -1/a
         ]
        
    return f



def odesolve(lamb,m,de):
    """
    Takes in:
        lamb = e_lamb(t)/ec(t0) at t=t0;
        m = e_m(t)/ec(t0) at t=t0;
        de = e_de(t)/ec(t0) at t=t0.
    Returns: 
        dlmpc = luminosity distance in Mpc;
        z = redshift under 2.
    
    """
    # Last value for a before results are to be considered 
    # invalid due to close proximity to big bang, a_d is set
    # somewhat arbitrarily - sometimes jumps over the result(?).
    a_d = 10e-6
    
    # Value to display redshift up to (and including). 
    ztrim = 6
    
    # Time (in 1/H0) to integrate until.  If this time isn't long enough for a to 
    # decrease to a_d then stoptime will be extended by time until a_d is reached.
    # 0.665 matter only, 0.96 standard m+de
    time = 0.6
    
    
    # Initial conditions at time = t0.
    a0 = 1.0        # scale factor
    a_dot0 = 1.0    # speed of expansion
    e_dash0m = m    # e_m(t)/ec(t0)
    e_dash0de = de  # e_de(t)/ec(t0)
    omega0m = m     # e_m(t)/ec(t)
    omega0de = de   # e_de(t)/ec(t)
    z0 = 0
    dl0 = 0
    
    # ODE solver parameters:
    abserr = 1.0e-8
    relerr = 1.0e-6
    numpoints = 250
    
    stoptime = -time # Integrating back in time as time now is t0.
    
    while True:
        # Create time samples for the ODE solver.
        t = [stoptime * tH * float(i) / (numpoints - 1) for i in range(numpoints)]
        
        # Pack up the initial conditions and eq of state parameters.
        v0 = [a0, a_dot0, e_dash0m, e_dash0de, omega0m, omega0de, z0, dl0]
        w = [w_m, w_de]
        
        # Call the ODE solver.
        vsol = odeint(firstderivs, v0, t, args=(w,lamb,), atol=abserr, rtol=relerr)
        
        # Remove unwanted results which are too close to big bang from the plot.
        # Separate results into their own arrays:
        a = vsol[:,0]
        a_dot = vsol[:,1]
        e_dashm = vsol[:,2]
        e_dashde = vsol[:,3]
        omegam = vsol[:,4]
        omegade = vsol[:,5]
        z = vsol[:,6]
        dl = vsol[:,7] * (1+z)
        dlmpc = dl * c_over_H0    # dl in Mega parsecs (= vsol[dl] * c/H0)
    
        
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
        z = z[:blowup]
        dl = dl[:blowup]
        dlmpc = dlmpc[:blowup]
        
        
        # Find when density of DE was equal to that of matter.  Rounding the 
        # omegas to avoid not finding an instance of equality.
        equiv = np.argmin(abs(omegam - omegade))    # Index of time of equivlence.
        equiv = np.asarray(equiv)                   # Converting to nparray.
        
        # Age of the universe.
        age = t_cut[np.argmin(t_cut)]
        age = -round(age, 2)
        
        # Trim redshift to remove meaningless values (z > 10) if such are present.
        if np.max(z) > ztrim:
            ztrim = np.where(z > ztrim)
            ztrim = np.asarray(ztrim)
            ztrim = ztrim[0,0]
            t_trim = t[:ztrim]
            z_trim = z[:ztrim]
            dl_trim = dl[:ztrim]
        
        # Plotting selected results:
        # a and a_dot vs time.
        while False:
            figure()
            xlabel('time in $H_0^{-1}$')
            grid(True)
            plot(t_cut, a, 'r', t_cut, a_dot, 'b', lw=1)
            legend((r'$a$', r'$\.a$'), prop=FontProperties(size=16))
            title('Cut results for $\omega$ = %s, $\lambda$ = %s, age = %s $H_0^{-1}$'
                  %(w,lamb,age))
            break
        
#        while False:
#            # Luminosity distance dl vs time.
#            figure()
#            xlabel('time in $H_0^{-1}$')
#            grid(True)
#            plot(t_cut, dl, 'tab:purple', lw=1)
#            title('$D_L$ vs time for $\omega$ = %s, $\lambda$ = %s,'
#              ' age = %s $H_0^{-1}$'%(w,lamb,age))
#            break
        
        # Luminosity distance dl vs redshift.
        while False:
            figure()
            xlabel('redshift $z$')
            grid(True)
            plot(z_trim, dl_trim, 'tab:green', lw=1)
            title('$D_L$ vs redshift for $\omega$ = %s, $\lambda$ = %s,'
                  ' age = %s $H_0^{-1}$'%(w,lamb,age))
            break
        
        while False:
            # Redshift vs time.
            figure()
            xlabel('time in $H_0^{-1}$')
            grid(True)
            plot(t_trim, z_trim, 'tab:pink', lw=1)
            title('Redshift evolution for $\omega$ = %s, $\lambda$ = %s,'
              ' age = %s $H_0^{-1}$'%(w,lamb,age))
            break
        
#        while False:
#            # Time vs redshift.
#            figure()
#            xlabel('z')
#            grid(True)
#            plot(z_trim, t_trim, 'tab:pink', lw=1)
#            title('Time evolution with $z$ for $\omega$ = %s, $\lambda$ = %s,'
#              ' age = %s $H_0^{-1}$'%(w,lamb,age))
#            break
        
#        while False:    # Looped to make it faster to switch plots on and off.
#            # e_dashm
#            figure()
#            xlabel('time in $H_0^{-1}$')
#            plot(t_cut, e_dashm, 'g', lw=1)
#            title('Cut $\epsilon_m \'$ for $\omega$ = %s'%(w))
#            
#            # e_dashde
#            figure()
#            xlabel('time in $H_0^{-1}$')
#            plot(t_cut, e_dashde, 'm', lw=1)
#            title('Cut $\epsilon_{DE} \'$ for $\omega$ = %s'%(w))
#            break
        
        # omegam and omegade calculated using the integrator:
        while False:
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
            break
        
#        while False:
#            # Verified omegam and omegade using e_dashm and e_dashde.
#            vomegam = e_dashm / (e_dashm + e_dashde)
#            vomegade = e_dashde / (e_dashm + e_dashde)
#            
#            figure()
#            xlabel('time in $H_0^{-1}$')
#            grid(True)
#            plot(t_cut,vomegam,t_cut,vomegade)
#            legend((r'$\Omega_m$', r'$\Omega_{DE}$'), prop=FontProperties(size=16))
#            title('Verified $\Omega_m$ and $\Omega_{DE}$ for'
#                  ' $\lambda$ = %s, age = %s $H_0^{-1}$'
#                  %(lamb,age))
#            break
            
        break
            
    
#    # Complete results with blow up resulting from a approaching big bang.
#    while False:  
#        figure()
#        xlabel('time in $H_0^{-1}$')
#        grid(True)
#        
#        # Plotting complete results.
#        plot(t, vsol[:,0], 'r', lw=1)
#        plot(t, vsol[:,1], 'b', lw=1)
#        
#        legend((r'$a$', r'$\.a$', r'$\'\epsilon$'), prop=FontProperties(size=16))
#        title('Complete results for $\omega$ = %s'%(w))
#    
#        break

    return dlmpc, z



def gnoise(mag, p):
    """
    adds p% noise to data given
    """
    sigma = p/100               # standard deviation
    mu = 0                      # mean of noise distribution
    noise = np.random.normal(mu,sigma,n)
    mag = mag + noise
    return mag



def msim(lamb, m, de, n, p):
    """
    Takes in:
        theta = 
            lamb (= e_lamb(t)/ec(t0) at t=t0),
            m (= e_m(t)/ec(t0) at t=t0),
            de (= e_de(t)/ec(t0) at t=t0),
            n (= dimensionless number of data points to be generated),
            p (= percentage of noise).
    Returns n apparent magnitudes mag and corresponding redshits z < 2, 
    offest by gaussian noise from gnoise.
    """
    dlmpc, z = odesolve(lamb,m,de)
    
    zmax = 2        # Largest meausured z for a supernovae is 2.
    maxindex_i = np.where(z > zmax) 
    maxindex_i = np.asarray(maxindex_i)   # Converting to np array.
       
    if maxindex_i.any():              # Check if instances of z > zmax exist.   
        max_index = maxindex_i[0,0]
    else:
        print('msim found no z above z = %s'%(zmax))

#    print('dlmpc before is: ',dlmpc)
    dlmpc = dlmpc[1:max_index]  # Avoiding the log(dlmpc=0) problem
    z = z[1:max_index]
#    print('dlmpc after [0:max_index] is: ',dlmpc)
    index_opts = range(len(z))
#    print('length of index_opts = ',len(index_opts))
#    print('index_opts are = ',index_opts)
    if len(index_opts) < n:
        print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
        print(' Decrease n, there are only %s z<2 points'%(len(z)))
        print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
    data_ind = random.sample(index_opts, n)

    listdlmpc = []
    listz = []
    
    for i in data_ind:
        listdlmpc.append(dlmpc[i])
        listz.append(z[i])

    dlmpc = listdlmpc#np.asarray(listdlmpc)
    z = listz#np.asarray(listz)
    
    # Calculating apparent magnitudes of supernovae at the simulated
    # luminosity distances using the distance modulus formula.
    mag = []
    for i in range(len(dlmpc)):
#        print('dl in Mpc = ',str(dlmpc[i]))
        mdistmod = 5 * log10(dlmpc[i]/10) + M
#        print('m from distance modulus = ',mdistmod)
        mag.append(mdistmod) 
#        print('mag = ',mag[i])
        
    mag = gnoise(mag, p)
    z = np.asarray(z)
    return z, mag

    
# Number of points to be simulated.
n = 100 # 100, 1000

# Percentage error on apparent magnitudes:
p = 10

# Standard cosmological parameters.
H0 = 1       # Hubble parameter at t=now
tH = 1.0/H0  # Hubble time
# Eq of state parameters for known fluids:
w_r = 1/3     # radiation
w_m = 0.0     # matter
w_de = -1.0   # cosmological constant (dark energy?)  

# Model specific parameters.  
# Interaction term, rate at which DE decays into matter.
lamb = 0
# t=t0 fraction of matter and dark energy compared to critical density.
m = 0.3
de = 0.7

# Empirical parameters.
M = -19                     # Absolute brightness of supernovae.
c_over_H0 = 4167 * 10**6    # c/H0 in parsecs

# Generating apparent magnitues mag at redshift z<2 (calculted from
# luminosity distances given by LambdaCMD with parameters stated above.
theta = lamb, m, de, n, p
z, mag = msim(lamb, m, de, n, p)
    

#z is x
#mag is y




###############################################################################
###############################################################################
###############################################################################
###############################    STATISTICS    ##############################
###############################################################################
###############################################################################
###############################################################################

try:
    timet0 = time.time()    # starting script timer
    
    
    # Input
    # "True" parameters are the cosmological ones.
    
    N = n        # number of datapoints
    sigma = p/100    # standard deviation
    mu = 0          # mean
    
    ndim, nwalkers = 3, 6
    nsteps = 1000
    burnin = 500
    
    
    # Functions
    def lnlike(theta, x, y, sigma):
        lamb, m, de, n, p = theta
        z, model = msim(lamb, m, de, n, p)
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
