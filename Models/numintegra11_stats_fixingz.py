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



3) What would happen if error were 1%? or 50%? How does it change my parameters
4) What if I generate a cosmology that has an interaction term? 
5) Choose another cosmology, generate fake data and try the same with omegam,
    omegade and lamb
(try the standard one with no lambda and see if you get lambda = 0 back)
6) Working towards distributions of omega lambda, omega matter and interaction term.
make numintegra8 give D_L in parsecs
put plots and calculation into separate functions




10 point sample does not generate enough points
plot results to see if they look right?



Make plot of everything vs everything through saving outcomes of variables and plot them after
change interaction term to gamma instead of labda, mindful of lambda teh cosm constant

Ask Sue Yang re debugging python to check type of error, underflow/overflow

Run msim with slightly different parameters, plot redshift vs dlmpc and see if 
anything looks super strange (some parameters might be unphysics)

NOT AN ASSIGNMENT
"""
# cosmo
import numpy as np
import random
from scipy.integrate import odeint
from math import log10
from pylab import figure, plot, scatter, xlabel, grid, legend, title
from matplotlib.font_manager import FontProperties
import corner
import emcee
import logging
import matplotlib.pyplot as pl
import scipy.optimize as op
import sys
import time

timet0 = time.time()    # starting script timer


def flist(start, stop, step):
    """
    Takes in:
        start, stop, step - integers or floats
    Returns:
        zlist - a list start to stop with step as increment
    """
#    print('-flist has been called')
    i = 0
    zlist = [start]
    
    while zlist[i] < stop:
        nextvalue = zlist[i] + step
        zlist.append(nextvalue)
        i += 1
        continue

    return zlist


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
#    print('@ firstderivs has been called')
    (a, a_dot, e_dashm, e_dashde, z, dl) = v #omegam, omegade, z, dl) = v
    (w_m, w_de) = w
    
    # fist derivatives of functions I want to find:
    f = [# a_dot (=scale factor)
         a_dot,
         # a_dotdot
         (-a/2) * (e_dashm * (1+3*w_m) + e_dashde * (1+3*w_de)), 
         
         # e'_dotm (=density(t) / crit density(t0))
#         -3 * (a_dot/a) * e_dashm * (1 + w_m -lamb/3 * a/a_dot * e_dashde/e_dashm ),
         -3 * (a_dot/a) * e_dashm + e_dashm * w_m - e_dashm* lamb/3 * e_dashm * a/a_dot * e_dashde,

         # e'_dotde
         -3 * (a_dot/a) * e_dashde * (1 + w_de +lamb/3 * a/a_dot),

         # z_dot (=redshift)
         -a_dot/a**2,
         # dl_dot (=luminosty distance)
         -1/a]
        
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
#    print('@@ odesolve has been called')
    # Last value for a before results are considered close enough to z = 2.
    a_d = 0.25
    
    # Time (in 1/H0) to integrate until.  If this time isn't long enough for a to 
    # decrease to a_d then stoptime will be extended by time until a_d is reached.
    # 0.665 matter only, 0.96 standard m+de
    time = 0.9
    
    
    # Initial conditions at time = t0.
    a0 = 1.0        # scale factor
    a_dot0 = 1.0    # speed of expansion
    e_dash0m = m    # e_m(t)/ec(t0)
    e_dash0de = de  # e_de(t)/ec(t0)
    z0 = 0
    dl0 = 0
    
    # ODE solver parameters:
    abserr = 1.0e-8
    relerr = 1.0e-6
    numpoints = 1000000
    
    stoptime = -time # Integrating back in time as time now is t0.
    
    while True:
        # Create time samples for the ODE solver.
        t = [stoptime * tH * float(i) / (numpoints - 1) for i in range(numpoints)]
#        print('time is : ',t[0])
        # Pack up the initial conditions and eq of state parameters.
        v0 = [a0, a_dot0, e_dash0m, e_dash0de, z0, dl0]
        w = [w_m, w_de]
        
        # Call the ODE solver. maxstep=5000000 added later to try and avoid 
        # ODEintWarning: Excess work done on this call (perhaps wrong Dfun type).
        vsol = odeint(firstderivs, v0, t, args=(w,lamb,), atol=abserr, rtol=relerr, mxstep=5000000)
        # vsol type is:  <class 'numpy.ndarray'>
                
        # Remove unwanted results which are too close to big bang from the plot.
        # Separate results into their own arrays:
        a = vsol[:,0]
        a_dot = vsol[:,1]
#        e_dashm = vsol[:,2]
#        e_dashde = vsol[:,3]
        z = vsol[:,4]
        dl = vsol[:,5] * (1+z)
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
        z = z[:blowup]
        dl = dl[:blowup]
        dlmpc = dlmpc[:blowup]
                
        
        # Age of the universe.
        age = t_cut[np.argmin(t_cut)]
        age = -round(age, 2)

        # Plotting selected results:
        # a and a_dot vs time.
        while True:
            figure()
            xlabel('time in $H_0^{-1}$')
            grid(True)
            plot(t_cut, a, 'r', t_cut, a_dot, 'b', lw=1)
            legend((r'$a$', r'$\.a$'), prop=FontProperties(size=16))
            title('Cut results for $\omega$ = %s, $\lambda$ = %s, age = %s $H_0^{-1}$'
                  %(w,lamb,age))
            break
        
        # Luminosity distance dl vs redshift.
        while True:
            figure()
            xlabel('redshift $z$')
            grid(True)
            plot(z, dl, 'tab:green', lw=1)
            title('$D_L$ vs redshift for $\omega$ = %s, $\lambda$ = %s,'
                  ' age = %s $H_0^{-1}$'%(w,lamb,age))
            break
        
        while True:
            # Redshift vs time.
            figure()
            xlabel('time in $H_0^{-1}$')
            pl.axis([0,-0.1,0,5])
            grid(True)
            plot(t_cut, z, 'tab:pink', lw=1)
            title('Redshift evolution for $\omega$ = %s, $\lambda$ = %s,'
              ' age = %s $H_0^{-1}$'%(w,lamb,age))
            break
            
        break
            
    
#    # Complete results with blow up resulting from a approaching big bang.
    while True:  
        figure()
        xlabel('time in $H_0^{-1}$')
        grid(True)
        
        # Plotting complete results.
        plot(t, vsol[:,0], 'r', lw=1)
        plot(t, vsol[:,1], 'b', lw=1)
        
        legend((r'$a$', r'$\.a$', r'$\'\epsilon$'), prop=FontProperties(size=16))
        title('Complete results for $\omega$ = %s'%(w))
        break

    return z, dlmpc


def msim(lamb, m, de, n, p, zpicks):
    """
    Takes in:
            lamb (= e_lamb(t)/ec(t0) at t=t0),
            m (= e_m(t)/ec(t0) at t=t0),
            de (= e_de(t)/ec(t0) at t=t0),
            n (= dimensionless number of data points to be generated),
            p (= percentage of noise)
            zpicks (=list of z to match the interpolated dlmpc to).
    Returns:
        mag (=list of n apparent magnitudes mag from corresponding redshits).
    """
#    print('@@@ msim has been called')
    z, dlmpc = odesolve(lamb,m,de)
    dlmpcinterp = np.interp(zpicks, z, dlmpc)

    # Calculating apparent magnitudes of supernovae at the simulated
    # luminosity distances using the distance modulus formula.
    mag = []
    for i in range(len(dlmpcinterp)):
        mdistmod = 5 * log10(dlmpcinterp[i]/10) + M
        mag.append(mdistmod) 
    return mag


def gnoise(mag, sigma, mu):
    """
    calculates and adds random p% Gaussian noise to each mag datapoint
    """
#    print('-gnoise has been called')
    noise = np.random.normal(mu,sigma,n)
#    print('noise from inside gnoise is = ', noise)
    mag = mag + noise
    return mag, noise


def lnlike(theta, n, p, zpicks, mag, sigma):
#    print('@@@@ lnlike has been called')
    lamb, m, de = theta
    model = msim(lamb, m, de, n, p, zpicks)
    inv_sigma2 = 1.0/(sigma**2)
    return -0.5*(np.sum((mag-model)**2*inv_sigma2 - np.log(inv_sigma2)))    

    
def lnprior(theta):
#    print(' lnprior has been called')
    lamb, m, de = theta
    if (-0.001 < lamb < 0.001 and 0 < m < 1.0 and 0.0 < de < 1.0):
        return 0.0
    return -np.inf
        

def lnprob(theta, n, p, zpicks, mag, sigma):
#    print('@@@@@ lnprob has been called')
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, n, p, zpicks, mag, sigma)


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

# Standard cosmological parameters.
H0 = 1       # Hubble parameter at t=now
tH = 1.0/H0  # Hubble time
# Eq of state parameters for known fluids:
w_r = 1/3     # radiation
w_m = 0.0     # matter
w_de = -1.0   # cosmological constant (dark energy?)  

# Model specific parameters.  
# Interaction term, rate at which DE decays into matter.
lamb_true = 0
# Fraction of matter and dark energy compared to critical density at t=t0.
m_true = 0.3
de_true = 0.7

# Empirical parameters.
M = -19                     # Absolute brightness of supernovae.
c_over_H0 = 4167 * 10**6    # c/H0 in parsecs





# Code:

# Picking redshifts to investigate.
zmin = 0.001
zmax = 2        # Largest meausured z for a supernovae is 2.
zinterval = (zmax - zmin) / (n*2)
z_opts = flist(zmin, zmax, zinterval)
zpicks = random.sample(z_opts, n)
zpicks = np.asarray(zpicks)


# Generating apparent magnitues mag at redshift z<2 (calculated from
# luminosity distances given by LambdaCMD with parameters stated above.
#theta = lamb, m, de
model = msim(lamb_true, m_true, de_true, n, p, zpicks)
model = np.asarray(model)

mag, noise = gnoise(model, sigma, mu)
#print('noise in code body is = ', noise)



try:
    # Finding a "good" place to start using alternative method to emcee.
    nll = lambda *args: -lnlike(*args)  # type of nll is: <class 'function'>
    result = op.minimize(nll, [lamb_true, m_true, de_true], 
                         args=(n, p, zpicks, mag, noise))
    lamb_ml, m_ml, de_ml = result["x"]    
    
        
    # Initializing walkers in a Gaussian ball around the max likelihood. 
    pos = [result["x"] + 1*np.random.randn(ndim) for i in range(nwalkers)]    
        
    
    # Sampler setup
    times0 = time.time()    # starting emcee timer
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(n, p, zpicks, mag, sigma))
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
    zl = zpicks
    figure()
    scatter(zl, model,color="r", lw=2, alpha=0.8)
    pl.errorbar(zpicks, model, yerr=sigma, fmt=".k")
    pl.show()
    
    # Best line of fit found by emcee.
    bi = np.argmax(sampler.lnprobability)   # index with highest post prob                                       
    lambbest = sampler.flatchain[bi,0]      # parameters with the highest 
    mbest = sampler.flatchain[bi,1]         # posterior probability
    debest = sampler.flatchain[bi,2]
    
    # plot of data with errorbars + model
    figure()
    pl.errorbar(zpicks, mag, yerr=sigma, fmt='o', alpha=0.3)
    modelt = msim(lamb_true, m_true, de_true, n, p, zpicks)
    model, = scatter(zpicks, modelt, lw='3', c='g')
    magbest = msim(lambbest, mbest, debest, n, p, zpicks)
    best_fit, = scatter(zpicks,magbest,lw='3', c='r')
    pl.legend([model, best_fit], ['Model', 'Best Fit'])
    pl.show
    
    
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