#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 19:05:26 2018

@author: BallBlueMeercat
"""
import numpy as np
import matplotlib.pyplot as plt
import time
import os
import zodesolve
import tools
import datasim
import results


res_recomb.py

#                mlfilename = f'maxlike_da_distrib_{test_key}_{noise}_{npoints}'
#                ml_da_list.append(ml_da[-1])
#                pickle.dump(ml_da_list, open(mlfilename, 'wb'))

#        plt.figure()
#        plt.title(f'$\sigma$ on data = {noise}, {npoints} SN Ia used')
#        plt.ylabel(r'$z = 1089$')
#        plt.xlabel('$(H_0 /c) * D_A$')
#        plt.grid(True)
#        plt.xlim(0.00285,0.00295)
#        for i in range(len(models)):
#            da_distrib = da_list[i]
#            z_array = np.ones(len(da_distrib))*1089
#            plt.scatter(da_distrib, z_array, s=msize[i], facecolors='none', edgecolors="C{}".format(i), label=models[i])
#        plt.locator_params(axis='x', nbins=4)
#        plt.yticks([])
#        plt.legend()
#        plt.show()
#
##        # Non-normalised histogram
##        plt.figure()
##        plt.title(f'$D_A$ at z = {zpicks[-1]},'
##                +f'\n $\sigma$ on data = {noise}, {npoints} SN Ia used')
##        plt.xlabel('$(H_0 /c) * D_A$')
##        for i in range(len(models)):
##            da_distrib = da_list[i]
##            plt.hist(da_distrib, normed=False, histtype='step',
##                     stacked=True, label=models[i])
##        plt.locator_params(axis='x', nbins=5)
##        plt.legend()
##        plt.show()
#
#        cutoff = 0.1
#        plt.figure()
#        smallest_da = 1
#        largest_da = 0
#        plt.title(f'scatter of x[1:] $D_A$ histogram $y > {cutoff}$'
#                  +'\n from all guessed parameter sets')
#        for i in range(len(models)):
#            face_color = 'none', "C{}".format(i)
#            da_distrib = da_list[i]
#            y, x = np.histogram(da_distrib, bins=n_bin)
##            print('x = ',x)
#            y_norm = y/max(y)
#            x = x[1:]
##            print('scatter x[1:] = ',x)
##            print('scatter y = ',y)
##            print('scatter y_norm = ',y_norm)
#            indices = np.where(y_norm > cutoff)
#            y_norm = y_norm[indices]
#            x = x[indices]
#            if min(x) < smallest_da:
#                smallest_da = min(x)
#            if max(x) > largest_da:
#                largest_da = max(x)
#            plt.scatter(x, y_norm, s=msize[i], facecolors=face_color[i], edgecolors="C{}".format(i), label=f'{models[i]}')
#        plt.xlim((smallest_da-0.000001),(largest_da+0.000001))
#        plt.locator_params(axis='x', nbins=5)
#        plt.legend()
#        plt.show()
#
#
#        plt.figure()
#        smallest_da = 1
#        largest_da = 0
#        plt.title(f'scatter of x[:-1] $D_A$ histogram $y > {cutoff}$'
#                  +'\n from all guessed parameter sets')
#        for i in range(len(models)):
#            face_color = 'none', "C{}".format(i)
#            da_distrib = da_list[i]
#            y, x = np.histogram(da_distrib, bins=n_bin)
##            print('x = ',x)
#            y_norm = y/max(y)
#            x = x[:-1]
##            print('scatter x[:-1] = ',x)
##            print('scatter y = ',y)
##            print('scatter y_norm = ',y_norm)
#            indices = np.where(y_norm > cutoff)
#            y_norm = y_norm[indices]
#            x = x[indices]
#            if min(x) < smallest_da:
#                smallest_da = min(x)
#            if max(x) > largest_da:
#                largest_da = max(x)
#            plt.scatter(x, y_norm, s=msize[i], facecolors=face_color[i], edgecolors="C{}".format(i), label=f'{models[i]}')
#        plt.xlim((smallest_da-0.000001),(largest_da+0.000001))
#        plt.locator_params(axis='x', nbins=5)
#        plt.legend()
#        plt.show()
#
#        plt.figure()
#        smallest_da = 1
#        largest_da = 0
#        plt.title('scatter of complete $D_A$ histogram'
#                  +'\n from all guessed parameter sets')
#        for i in range(len(models)):
#            face_color = 'none', "C{}".format(i)
#            da_distrib = da_list[i]
#            y, x = np.histogram(da_distrib, bins=n_bin)
##            print('scatter x = ',x)
#            y_norm = y/max(y)
#            x = x[:-1]
##            print('scatter x[:-1] = ',x)
##            print('scatter y = ',y)
##            print('scatter y_norm = ',y_norm)
#            if min(x) < smallest_da:
#                smallest_da = min(x)
#            if max(x) > largest_da:
#                largest_da = max(x)
#            plt.scatter(x, y_norm, s=msize[i], facecolors=face_color[i], edgecolors="C{}".format(i), label=f'{models[i]}')
#        plt.xlim((smallest_da-0.000001),(largest_da+0.000001))
#        plt.locator_params(axis='x', nbins=5)
#        plt.legend()
#        plt.show()
#
#
#
#
#
##        for i in range(len(models)):
##            plt.figure()
##            da_distrib = da_list[i]
##            y, x = np.histogram(da_distrib)
##            y_norm = y/max(y)
##            x = x[1:]
##            plt.title("$D_A$'s using all guessed parameters")
##            plt.plot(x, y_norm, label=f'{models[i]}')
##            plt.xlim(0.00285,0.003)
##            plt.locator_params(axis='x', nbins=5)
##            plt.legend()
##        plt.show()
#
##        plt.figure()
##        plt.title(f'$\mu$ and $\sigma$ of the $D_A$ distribution'
##                  +f'\n $\sigma$ on data = {noise}, {npoints} SN Ia used')
##        plt.xlabel('$(H_0 /c) * D_A$')
##        plt.ylabel(r'$z = 1089$')
##        for i in range(len(models)):
##            da_distrib = np.asarray(da_list[i])
##            da_mean = np.mean(da_distrib)
##            da_sd = np.std(da_distrib)
##            if da_sd < da_mean/100:
##                print(f'{test_key} da_sd (={da_sd}) < mean/100 (={da_mean/100})')
###            plt.errorbar(da_mean, yaxis_tick[i], xerr=da_sd, fmt='-o',
###        color=c[i], ecolor=ec[i], elinewidth=3, capsize=0, label=models[i])
##            plt.errorbar(da_mean, yaxis_tick[i], xerr=da_sd, fmt='o', label=models[i])
##        plt.locator_params(axis='x', nbins=5)
##        plt.ylim(-4,8)
##        plt.yticks([])
##        plt.legend()
##        plt.show()


datasim.py

#def magn(params, data, firstderivs_key, plot_key=False):
#    """
#    Finding matter density m, interaction gamma.
#
#    Takes in:
#            params = dictionary with true parameters;
#            zpicks = list of redshifts to integrate over, in accending order;
#            firstderivs_key = string, indicates which firstderivs to integrate;
#            plot_key = Boolean, to plot or not to plot model figures;
#    Returns:
#        mag = np.ndarray of apparent mag corresponding to input redshits.
#    """
##    print('@@@ magn has been called')
#    if firstderivs_key == 'LCDM':
#        params['gamma'] = 0
#        del params['gamma']
#
#    zpicks = data['zpicks']
#
#    # Absolute brightness of supernovae.
#    M = -19
#
#    dlpc, plot_var = zodesolve.zodesolve(params, zpicks, firstderivs_key)
#
#    # Calculating apparent magnitudes of supernovae at the simulated
#    # luminosity distances using the distance modulus formula.
#    mag = 5 * np.log10(dlpc/10) + M
#
#    if plot_key:
#        # Checking evolution of the model.
#        import plots
#        plots.modelcheck(mag, zpicks, plot_var, firstderivs_key)
#
#    return mag

#def magn(params, data, firstderivs_key, plot_key=False):
#    """
#    Finding matter density m, alpha, beta, interaction gamma.
#    Takes in:
#            params = dictionary with true parameters;
#            zpicks = list of redshifts to integrate over, in accending order;
#            firstderivs_key = string, indicates which firstderivs to integrate;
#            plot_key = Boolean, to plot or not to plot model figures;
#    Returns:
#        mag = np.ndarray of apparent mag corresponding to input redshits.
#    """
##    print('@@@ magn has been called')
#    if firstderivs_key == 'LCDM':
#        params['gamma'] = 0
#        del params['gamma']
#
#    zpicks = data['zpicks']
#    x1 = data['x1']
#    colour = data['colour']
#
#    # Absolute brightness of supernovae.
#    M = params['M']
#    alpha = params['alpha']
#    beta = params['beta']
#
#    dlpc, plot_var = zodesolve.zodesolve(params, zpicks, firstderivs_key)
#
#    # Calculating apparent magnitudes of supernovae at the simulated
#    # luminosity distances using the distance modulus formula.
#    mag = 5 * np.log10(dlpc/10) + M - alpha*x1 +beta*colour
#
#    if plot_key:
#        # Checking evolution of the model.
#        import plots
#        plots.modelcheck(mag, zpicks, plot_var, firstderivs_key)
#
#    return mag

# Slow mag calculation
#    # Calculating apparent magnitudes of supernovae at the simulated
#    # luminosity distances using the distance modulus formula.
#    mag = []
#    for i in range(len(dlpc)):
#        if dlpc[i] == 0:
#            magnitude = M
#        else:
#            # magnitude from the distance modulus formula
#            magnitude = 5 * math.log10(dlpc[i]/10) + M
#        mag.append(magnitude)

evaluator.py

#dataname = 'mag_z_LCDM_1000_sigma_'+str(sigma)

# Number of datapoints to be simulated
#npoints = 1000
#zmax = 2

#    Making data (mag and z)
#dataname = 'mag_z_'+data_key+'_'+str(npoints)+'_sigma_'+str(sigma)
#datasim.makensavemagnz(m_true, g_true, mu, sigma, zpicks, data_key, dataname)

#    Making redshifts to use in this script.
#zpicks = datasim.redshift_picks(0.005, zmax, npoints)


ln.py

#def lnlike(theta, data, sigma, firstderivs_key, ndim):
#    '''
#    Finding matter density m, interaction gamma.
#    '''
#    mag = data['mag']
#
#    params = {}
#    if ndim == 1:
#        params = {'m':theta}
#    elif ndim == 2:
#        params = {'m':theta[0],'gamma':theta[1]}
#
#    model = magn(params, data, firstderivs_key)
#    var = sigma**2
#    return -0.5*np.sum((mag-model)**2 /var +0.5*np.log(2*np.pi*var))


#def lnlike(theta, data, sigma, firstderivs_key, ndim):
#    '''
#    Finding matter density m, corrected absolute mag M, interaction gamma.
#    '''
#    mag = data['mag']
#
#    params = {}
#    if ndim == 2:
#        params = {'m':theta[0], 'M':theta[1]}
#    elif ndim == 3:
#        params = {'m':theta[0],'M':theta[1], 'gamma':theta[2]}
#
#    model = magn(params, data, firstderivs_key)
#    var = sigma**2
#    return -0.5*np.sum((mag-model)**2 /var +0.5*np.log(2*np.pi*var))

#def lnprior(theta, key):
#    '''
#    Finding matter density m, interaction gamma.
#    '''
#
#    if key == 'LCDM':
#        m = theta
#        if 0 < m < 1 or m == 1:
#            return 0.0
#    elif key == 'late_int' or 'heaviside_late_int' or 'late_intxde':
#        m, gamma = theta
#        if (0 < m < 1 or m == 1) and -1.45 < gamma < 0.2:
#            return 0.0
#    elif key == 'rdecay':
#        m, gamma = theta
#        if (0 < m < 1 or m == 1) and -10 < gamma < 0:
#            return 0.0
#    elif key == 'interacting':
#        m, gamma = theta
#        if (0 < m < 1 or m == 1) and abs(gamma) < 1.45:
#            return 0.0
#    elif key == 'expgamma':
#        m, gamma = theta
#        if (0 < m < 1 or m == 1) and abs(gamma) < 25:
#            return 0.0
#    elif key == 'zxxgamma' or 'gammaxxz':
#        m, gamma = theta
#        if (0 < m < 1 or m == 1) and 0 < gamma < 10:
#            return 0.0
#    else:
#        m, gamma = theta
#        if (0 < m < 1 or m == 1) and abs(gamma) < 10:
#            return 0.0
#
#    return -np.inf
#
#
#def lnprior(theta, key):
#    '''
#    Finding matter density m, corrected absolute mag M, interaction gamma.
#    '''
#
#    Mmin = -20
#
#    Mmax = -18
#
#    if key == 'LCDM':
#        m, M = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax:
#            return 0.0
#    elif key == 'late_int' or 'heaviside_late_int' or 'late_intxde':
#        m, M, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and -1.45 < gamma < 0.2:
#            return 0.0
#    elif key == 'rdecay':
#        m, M, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and -10 < gamma < 0 :
#            return 0.0
#    elif key == 'interacting':
#        m, M, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(gamma) < 1.45:
#            return 0.0
#    elif key == 'expgamma':
#        m, M, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(gamma) < 25 :
#            return 0.0
#    elif key == 'zxxgamma' or 'gammaxxz':
#        m, M, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and 0 < gamma < 10:
#            return 0.0
#    else:
#        m, M, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(gamma) < 10:
#            return 0.0
#
#    return -np.inf
#
#
#def lnprior(theta, key):
#    '''
#    Finding matter density m, absolute M, alpha, beta, interaction gamma.
#    '''
#
#    Mmin, Mmax = -20, -18
#    amax = 5
#    bmax = 5
#
#    print('key ln prior gets is = ',key)
#
#    if key == 'LCDM':
#        m, M, alpha, beta = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(alpha) < amax and abs(beta) < bmax:
#            return 0.0
#    elif key == 'late_int' or key == 'heaviside_late_int' or key == 'late_intxde':
#        m, M, alpha, beta, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(alpha) < amax and abs(beta) < bmax and -1.45 < gamma < 0.2:
#            return 0.0
#    elif key == 'rdecay':
#        m, M, alpha, beta, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(alpha) < amax and abs(beta) < bmax and -10 < gamma < 0 :
#            return 0.0
#    elif key == 'interacting':
#        m, M, alpha, beta, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(alpha) < amax and abs(beta) < bmax and abs(gamma) < 1.45:
#            return 0.0
#    elif key == 'expgamma':
#        m, M, alpha, beta, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(alpha) < amax and abs(beta) < bmax and abs(gamma) < 25 :
#            return 0.0
#    elif key == 'zxxgamma' or key == 'gammaxxz':
#        m, M, alpha, beta, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(alpha) < amax and abs(beta) < bmax and 0 < gamma < 10:
#            return 0.0
#    elif key == 'exotic':
#        m, M, alpha, beta, gamma, zeta = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(alpha) < amax and abs(beta) < bmax and 0 < gamma < 10 and 0 < zeta < 10:
#            return 0.0
#    else:
#        m, M, alpha, beta, gamma = theta
#        if (0 < m < 1 or m == 1) and Mmin < M < Mmax and abs(alpha) < amax and abs(beta) < bmax and abs(gamma) < 10:
#            return 0.0
#
#    return -np.inf

zodesolve.py

#    theta = dl/da
#    print('theta',theta[-1])
#    print(model,'redshift = ',z[-1],'da =',da[-1])
#    plt.figure()
#    plt.title('Angular diameter distance vs redshift')
#    plt.xlabel('z')
#    plt.ylabel('D_A')
#    plt.plot(z, da, label='da')
#    plt.legend()
#
#    plt.figure()
#    plt.title('Angular diameter vs redshift')
#    plt.xlabel('z')
#    plt.ylabel(r'$\theta$')
#    plt.plot(z, theta, label='theta')
#    plt.legend()
#    plt.show()

#        plt.figure()
#        plt.title('Angular diameter distance evolution')
#        plt.xlabel('z')
#        plt.ylabel(r'$ \left( \frac{H_0}{c} \right) d_A $', fontsize=15, labelpad=10)
#        plt.plot(z, da)

#        plt.figure()
#        plt.title(r'$d_V$ evolution')
#        plt.xlabel('z')
#        plt.ylabel(r'$ d_V (z)$ [Mpc]')
#        plt.grid(True)
#        plt.plot(z, dV)


#        Dv = ((1+z)**2 * daMpc**2 * c*z/Hz)**(1/3)
#        plt.figure()
#        plt.title(r'$D_v$ evolution')
#        plt.xlabel('z')
#        plt.ylabel(r'$ D_v (z)$ [Mpc]')
#        plt.grid(True)
#        plt.plot(z, Dv)

#        # Calculating the sound horizon
#        ombar_m = vsol[1:,2]
#        ombar_baryon =  ombar_m*0.04 #0.0125
##        ombar_baryon = values[(fluid_in-1)]
#        s = 44.5 * np.log(9.83 / ombar_m) / (1 +10 * ombar_baryon**(3/4))**(1/2)
#        plt.figure()
#        plt.title('Sound horizon')
#        plt.xlabel(r'$z$')
#        plt.ylabel('Physical length in Mpc')
#        plt.grid(True)
#        plt.plot(z, s, label=r'$s_H$')
#        plt.legend()

#        plt.show()


#def odesolve(params, zpicks, firstderivs_key):
#    """
#    Takes in:
#        gamma = interaction constant;
#        m = e_m(t)/ec(t0) at t=t0;
##        de = e_de(t)/ec(t0) at t=t0.
#    Returns:
#        z = numpoints number of redshifts zmin<z<zmax;
#        dlpc = luminosity distance in pc.
#
#    """
##    print('@@ zodesolve has been called')
#
#    # Inserting 0 at the front of redshifts to allow initial conditions.
#    zpicks = [0.0] + zpicks
#
#    # Standard cosmological parameters.
#    H0 = 1
#    c_over_H0 = 4167 * 10**6    # c/H0 in parsecs
#
#    # Initial conditions at z = 0 (now).
#    dl0 = 0             # luminosity distance
#    rho_c0 = H0**2      # critical density
#    ombar_m0 = params.get('m', 0)                        # e_m(z)/ec(z=0)
#    gamma = params.get('gamma',0)
#    ombar_de0 = params.get('de', rho_c0/rho_c0 -ombar_m0)# e_de(z)/ec(z=0)
#
#    # ODE solver parameters:
#    abserr = 1.0e-8
#    relerr = 1.0e-6
#
#    # Pack up the initial conditions and eq of state parameters.
#    v0 = [ombar_m0, ombar_de0, dl0]
#
#    # Extracting the parsed mode of interaction.
#    firstderivs_function = firstderivs_functions.get(firstderivs_key,0)
#    if firstderivs_function == 0:
#        print("firstderivs_functions dict didn't have the key zodeosolve asked for")
#
#    # Call the ODE solver.
#    vsol = odeint(firstderivs_function, v0, zpicks, args=(gamma,H0),
#                  atol=abserr, rtol=relerr)
#
#    # Separate results into their own arrays:
#    z = np.asarray(zpicks)
#    z = z[1:]
#    ombar_m = vsol[1:,0]
#    ombar_de = vsol[1:,1]
#    dl = vsol[1:,2] * (1+z)  # in units of dl*(H0/c)
#    dlpc = dl * c_over_H0    # dl in parsecs (= vsol[dl] * c/H0)
#
#    plot_var = dlpc, dl, ombar_m, gamma, ombar_de, ombar_m0, ombar_de0
#
#    return dlpc, plot_var

Bfactor.py

#    def randh(self):
#        """
#        Generate from the heavy-tailed distribution.
#        """
#        a = np.random.randn()
#        b = np.random.rand()
#        t = a/np.sqrt(-np.log(b))
#        n = np.random.randn()
#        return 10.0**(1.5 - 3*np.abs(t))*n
#
#    def wrap(self, x, a, b):
#        assert b > a
#        return (x - a)%(b - a) + a

#import six
#import sys
## Run the postprocessing to get marginal likelihood and generate posterior
#samples logZdnest4, infogaindnest4, plot = dnest4.postprocess()
#
#postsamples = np.loadtxt('posterior_sample.txt')
#
#print(six.u('Marginalised evidence is {}'.format(logZdnest4)))
#
#print('Number of posterior samples is {}'.format(postsamples.shape[0]))
#
## plot posterior samples (if corner.py is installed)
#try:
#    import matplotlib as mpl
#    mpl.use("Agg") # force Matplotlib backend to Agg
#    import corner # import corner.py
#except ImportError:
#    sys.exit(1)
#
#m = 0.3
#g=0
#fig = corner.corner(postsamples, labels=[r"$m$", r"$c$"], truths=[m, g])
#fig.savefig('DNest4.png')

# LCDM
#log(Z) = -1622866.8534441872
#Information = 14.078678027261049 nats.
#Effective sample size = 129.22232212112772
#time 297min 50s

#log(Z) = -1622866.790641218
#Information = 13.905435690656304 nats.
#Effective sample size = 167.73507536834273
#time 34 min

#rdecay
#log(Z) = -1622866.8177826053
#Information = 13.970533838961273 nats.
#Effective sample size = 85.54638980461822
#Sampling time:   37min 5s

############ 0.01 sigma data
#Hdecay
#Sampling time:   38min 57s
#log(Z) = -1158842.6212481956
#Information = 26.434626991627738 nats.
#Effective sample size = 116.96489141639181

#edecay
#Sampling time:   45min 57s
#log(Z) = -49925.259544267705
#Information = 19.683044903278642 nats.
#Effective sample size = 162.7283801030449

#LCDM
#Sampling time:   31min 52s
#log(Z) = -1622866.7230921672
#Information = 13.870062695583329 nats.
#Effective sample size = 178.67158154325102

############ 0.1 sigma data

#Hdecay
#Sampling time:   26min 26s
#data = mag_z_LCDM_1000_sigma_0.1
#sigma = 0.1
#log(Z) = -11392.938034458695
#Information = 16.85219457607309 nats.
#Effective sample size = 216.9365844057018

#rdecay
#Sampling time:   25min 4s
#data = mag_z_LCDM_1000_sigma_0.1
#sigma = 0.1
#log(Z) = -16069.573635539238
#Information = 8.730470507740392 nats.
#Effective sample size = 172.4071834775586

#LCDM
#Sampling time:   23min 45s
#data = mag_z_LCDM_1000_sigma_0.1
#sigma = 0.1
#log(Z) = -16070.356294581907
#Information = 9.449718869756907 nats.
#Effective sample size = 142.47418654118337

plot.py multi_modelcheck

#    # Evolution of the angular diameter distance.
#    plt.figure()
##    plt.title('Angular diameter distance evolution'+f'\n Models: {model_names}')
#    plt.xlabel('z')
#    plt.ylabel(r'$ (H_0/c) \ D_A $')
#    for i in range(len(da)):
#        plt.plot(zpicks, da[i], label='$\gamma$ = %d, %s'%(int_terms[i], keys[i]))
#    plt.legend()

#     # Evolution of dV.
#    plt.figure()
##    plt.title(r'$d_V$ evolution, models: %s.'%(model_names))
#    plt.xlabel('z')
#    plt.ylabel(r'$ d_V (z)$ [Mpc]')
#    for i in range(len(dV)):
#        plt.plot(zpicks, dV[i], label='$\gamma$ = %d, %s'%(int_terms[i], keys[i]))
#    plt.legend()

#    # Luminosity distance vs redshift.
#    plt.figure()
#    plt.xlabel('$z$')
#    plt.ylabel('$D_L$*($H_0$/c)')
#    for i in range(len(dl)):
#        plt.plot(zpicks, dl[0], label='$\gamma$ = %d, %s'%(int_terms[i], keys[i]))
##    plt.title('$d_L$ evolution, models: %s.'%(model_names))
#    plt.legend()
#
#    # H vs redshift.
#    plt.figure()
#    plt.xlabel('$z$')
#    plt.ylabel('H')
#    for i in range(len(Hz)):
#        plt.plot(zpicks, Hz[i], label='$\gamma$ = %d, %s'%(int_terms[i], keys[i]))
##    plt.title('H evolution, models: %s.'%(model_names))
#    plt.legend()
#
#    # Scale factor vs redshift.
#    plt.figure()
#    plt.xlabel('$z$')
#    plt.ylabel('a')
#    for i in range(len(a)):
#        plt.plot(zpicks, a[0], label='$\gamma$ = %d, %s'%(int_terms[i], keys[i]))
##    plt.title('Scale factor evolution, models: %s.'%(model_names))
#    plt.legend()
#
#    # Scale factor vs age.
#    plt.figure()
#    plt.xlabel('Age')
#    plt.ylabel('a')
#    for i in range(len(age)):
#        plt.plot(age[i], a[i], label='$\gamma$ = %d, %s'%(int_terms[i], keys[i]))
##    plt.title('Scale factor evolution, models: %s.'%(model_names))
#    plt.legend()
#
#    # Redshift vs age.
#    plt.figure()
#    plt.xlabel('Age')
#    plt.ylabel('$z$')
#    for i in range(len(age)):
#        plt.plot(age[i], zpicks, label='$\gamma$ = %d'%(int_terms[i]))
##    plt.title('Redshift evolution, models: %s.'%(model_names))
#    plt.legend()