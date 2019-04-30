#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 20:20:03 2018

@author: BallBlueMeercat
"""
import time
import dnest4
import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt
import pickle
import datasim
import results
import tools

import matplotlib as mpl
mpl.style.use('default') # has to be switched on to set figure size
mpl.style.use('fivethirtyeight')
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['grid.color'] = 'white'

timed = False
plot = False

print('Bfactor')

# Loading pantheon SN Ia data:
dataname = 'pantheon'
#dataname = 'synth'

if dataname == 'pantheon':
    import pandas as pd
    print('-----Using pantheon')
    # Pantheon data:
    pantheon = pd.read_csv('./data/lcparam_full_long.txt', sep=" ")
    pantheon.set_index('name', inplace=True)
    pantheon.sort_values('zhel', inplace=True)
    mag = pantheon.mb.values
    sigma = pantheon.dmb.values
    zpicks = pantheon.zhel.values
    data_dic = {'mag':mag, 'zpicks':zpicks}
#    plt.figure()
#    plt.title('Pantheon')
#    plt.scatter(zpicks, mag)
#    plt.show(block=False)
elif dataname == 'synth':
    mu, sigma = 0.0, 0.001    # Mean and standard deviation of the noise on the data.
    npoints = 1048
    min_z = 0.01012
    max_z = 2.26
    # Loading artificial LCDM SN Ia data:
    from pathlib import Path
    dataname = f'data/{npoints}_{max_z}_sigma_{sigma}.p'
    my_file = Path(dataname)
    if my_file.is_file():
        with open(dataname,'rb') as rfp: zpicks, mag = pickle.load(rfp)
        plt.figure()
        plt.title(f'Artificial data being used N={len(zpicks)}, $\sigma$={sigma}')
        plt.scatter(zpicks, mag)
        plt.show(block=False)
        data_dic = {'mag':mag, 'zpicks':zpicks}
    else:
        print(f'failed to get zpicks, mag from {dataname}')
        print('generating zpicks and mag')
        # Generating redshifts.
        zpicks = np.random.uniform(low=min_z, high=max_z, size=(npoints,))
        zpicks = np.sort(zpicks, axis=None)
        if zpicks[-1] != max_z:
            zpicks[-1] = max_z
        data_dic = {'zpicks':zpicks}
        # Generating LCDM mag and da.
        mag, da = datasim.magn(['Mcorr', 'matter'], np.array([-19.3, 0.3]), data_dic, 'LCDM')
        # Adding noise to LCDM mag.
        mag = datasim.gnoise(mag, mu, sigma)

        plt.figure()
        plt.title(f'Artificial data N={len(zpicks)}, $\sigma$={sigma}')
        plt.scatter(zpicks, mag)
        plt.show(block=False)

        data = zpicks, mag
        pickle.dump(data, open(dataname, 'wb'))
        data_dic = {'mag':mag, 'zpicks':zpicks}


class Model(object):
    """
    Specify the model in Python.
    """
    def __init__(self, names, int_lim, fluid_number):
        """
        Parameter values *are not* stored inside the class
        """
        # Prior on corrected magnitude.
        self.M_min = -20
        self.M_max = -18

        self.fluid_number = fluid_number
        self.names = names
        self.int_lim = int_lim

    def from_prior(self):
        """
        Unlike in C++, this must *return* a numpy array of parameters.
        """
        M = 1E3*rng.rand()
        M = dnest4.wrap(M, self.M_min, self.M_max)
        M = np.array([M])
        # Sampling fluids from prior (uniform distribution between 0 and 1).
        fluids = [rng.rand() for i in range(0,fluid_number)]

        if self.int_lim:
            int_terms = np.zeros(len(self.int_lim))
            for i in range(len(self.int_lim)):
                term = 1E3*rng.rand()
                term = dnest4.wrap(term, self.int_lim[i][0], self.int_lim[i][1])
                int_terms[i] = term
            return np.concatenate((M, fluids, int_terms))
        return np.concatenate((M, fluids))

    def perturb(self, theta):
        """
        Unlike in C++, this takes a numpy array of parameters as input,
        and modifies it in-place. The return value is still logH.
        """
        logH = 0.0
        pic = rng.randint(len(theta))
        # Note the difference between dnest4.wrap in Python and
        # DNest4::wrap in C++. The former *returns* the wrapped value.
        if pic == 0:
            theta[pic] += dnest4.randh()
            theta[pic] = dnest4.wrap(theta[pic], self.M_min, self.M_max)
        elif 0 < pic < (fluid_number+1):
            theta[pic] += dnest4.randh()
            theta[pic] = dnest4.wrap(theta[pic], 0.0, 1.0)
        elif fluid_number < pic:
            i = pic - fluid_number - 1   # index of interaction term
            theta[pic] += dnest4.randh()
            theta[pic] = dnest4.wrap(theta[pic],
                  self.int_lim[i][0], self.int_lim[i][1])
        return logH

    def log_likelihood(self, theta):
        """
        Gaussian sampling distribution.
        """
        model_mag, model_da = datasim.magn(self.names, theta, data_dic, key)
        var = sigma**2.0
        like = -0.5*np.sum((mag-model_mag)**2.0 /var +np.log(2.0*np.pi*var))
        return like

firstderivs_functions = [None
#            ,'rainbow' # speed 5
#            ,'niagara' # speed 4
#            ,'kanangra' # speed 4
#            ,'waterfall' # speed 3
#            ,'stepfall' # speed 3
#            ,'exotic'
#            ,'late_intxde'
#            ,'heaviside_late_int'
#            ,'heaviside_sudden' # didn't converge
#            ,'late_int'
#            ,'expgamma'
#            ,'txgamma'
#            ,'zxgamma'
#            ,'gamma_over_z'
#            ,'zxxgamma'        # gamma forced positive in firstderivs
#            ,'gammaxxz'        # gamma forced positive in firstderivs
#            ,'rdecay_m'
#            ,'rdecay_de'
#            ,'rdecay_mxde'
#            ,'rdecay'           # didn't converge
#            ,'interacting'
            ,'LCDM'
#            ,'rLCDM'
            ]

for key in firstderivs_functions:
    if key:
        print(key)
        names, int_lim, speed = tools.names_intlim_speed(key)

        speed = 0
        if int_lim:
            fluid_number = len(names) - 1 - len(int_lim)
        else:
            fluid_number = len(names) - 1

        # Create a model object and a sampler.
        model = Model(names, int_lim, fluid_number)
        sampler = dnest4.DNest4Sampler(model,
                            backend=dnest4.backends.CSVBackend(".",sep=" "))

        if speed == 6: #  extra long
            max_lvl,nstep,new_lvl,n_per_step,th_step = 60,1000,10000,10000,100

        elif speed == 5: #  extra long
            max_lvl,nstep,new_lvl,n_per_step,th_step = 50,1000,10000,10000,100

        elif speed == 4: # extra long
            max_lvl,nstep,new_lvl,n_per_step,th_step = 40,1000,10000,10000,100

        elif speed == 3: # LONG
            max_lvl,nstep,new_lvl,n_per_step,th_step = 30,1000,10000,10000,100

        elif speed == 2: # MEDIUM
            max_lvl,nstep,new_lvl,n_per_step,th_step = 30,1000,1000,1000,100

        elif speed == 1: # SHORT
            max_lvl,nstep,new_lvl,n_per_step,th_step = 30,100,100,100,10

        elif speed == 0: # sampling from prior
            max_lvl,nstep,new_lvl,n_per_step,th_step = 1,1000,100,100,10

        # Set up the sampler. num_per_step can be down to a few thousand.
        gen = sampler.sample(max_num_levels=max_lvl, num_steps=nstep,
                        new_level_interval=new_lvl, num_per_step=n_per_step,
                        thread_steps=th_step, num_particles=5,
                        lam=10, beta=100, seed=1234)

        if timed:
            import cProfile, pstats, io
            pr = cProfile.Profile()
            pr.enable()

        ti = time.time()

        # Do the sampling (one iteration here = one particle save).
        for i, sample in enumerate(gen):
#            print("# Saved {k} particles.".format(k=(i+1)))
            pass
        tf = time.time()

        if timed:
            pr.disable()
            s = io.StringIO()
            sortby = 'tottime'
            ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
            ps.print_stats()
            print (s.getvalue())

        dnest_time = tools.timer('Bfactor', ti, tf)

        print('testing =',key)
        print('data =', dataname)
        print('sigma =', sigma)

        # Histogram of parameters found by DNest4.
        array = np.loadtxt('sample.txt')

        DNest_distr = {}

        if plot: # ['light red', 'berry', 'purple', 'black']
            hue = ['coral', 'amber', 'apple', 'aquamarine',
                'raspberry', 'green blue', 'deep blue', 'emerald',
                'blue violet', 'dark violet', 'yellow orange',
                'light red', 'berry', 'coral', 'amber', 'apple',
                'aquamarine', 'raspberry', 'green blue','deep blue',
                'emerald', 'blue violet', 'dark violet', 'black']
#            ndim = len(array[0,:])
            ndim = len(names)
            for i in range(ndim):
                name = names[i]
                plt.figure()
#                plt.title(name)
                plt.hist(array[:,i], color='xkcd:'+hue[i])
                plt.locator_params(axis='x', nbins=5)
                distribution = array[:,i]
                # Standard deviation and mean of the DNest distribution.
                DNest_distr[name+'_sd'] = np.std(distribution)
                DNest_distr[name+'_mean'] = np.mean(distribution)
                DNest_distr[name] = array[:,i]
        plt.show(block=False)

        # Run the postprocessing
        info = dnest4.postprocess()

        if speed > 1:

            f = open('brief.txt','w')
            f.write(dnest_time +'\n'
                    +'model = '+key +'\n'
                    +'data = '+ dataname +'\n'
                    +'sigma = '+str(sigma) +'\n'
                    +'log(Z) = '+str(info[0]) +'\n'
                    +'Information = '+str(info[1]) +'\n'
                    +'speed = '+str(speed))
            f.close()

            pickle.dump(info[0], open('evidence.p', 'wb'))
            # Moving output .txt files into a run specific folder.
            results.relocate('evidence.p', speed, key)
            results.relocate('levels.txt', speed, key)
            results.relocate('posterior_sample.txt', speed, key)
            results.relocate('sample_info.txt', speed, key)
            results.relocate('sample.txt', speed, key)
            results.relocate('sampler_state.txt', speed, key)
            results.relocate('weights.txt', speed, key)
            results.relocate('brief.txt', speed, key)
            results.relocate('plot_1.pdf', speed, key)
            results.relocate('plot_2.pdf', speed, key)
            results.relocate('plot_3.pdf', speed, key)

if speed > 1:
    import os
    os.system("say 'the Befactor script has finished my love' &")