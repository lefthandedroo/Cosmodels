#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 20:20:03 2018

@author: BallBlueMeercat
"""
import time
import dnest4
import numpy as np
import pandas as pd
import numpy.random as rng
import matplotlib.pyplot as plt
import pickle
import datasim
import results
import tools

speed = 2       # From prior = 0, short = 1, medium = 2, long = 3.
timed = False
plot = True
mu, sigma = 0.0, 0.07    # Mean and standard deviation of the noise on the data.

# Loading pantheon SNIa data:
dataname = './data/lcparam_full_long.txt'
pantheon = pd.read_csv(dataname, sep=" ")
# Reading each txt file column of interest as numpy.ndarray.
mag = pantheon.mb.values
zpicks = pantheon.zhel.values
# Stacking np.arrays together and sorting by accending redshift.
data = np.stack((mag,zpicks), axis=0)
data.sort(axis=-1)
mag = data[0]
zpicks = data[-1]
#zpicks = zpicks.tolist()
data_dic = {'mag':mag, 'zpicks':zpicks}

## Generating LCDM data.
#names = ['Mcorr', 'matter']
#values = np.array([-19.3, 0.3])
#mag = datasim.noisy_mag(mu, sigma, names, values, data_dic, 'LCDM')
#data_dic = {'mag':mag, 'zpicks':zpicks}

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
        return -0.5*np.sum((mag-model_mag)**2.0 /var +0.5*np.log(2.0*np.pi*var))


firstderivs_functions = [None
#            ,'stepfall'
#            ,'waterfall'
#            ,'exotic'
            ,'late_intxde'
            ,'heaviside_late_int'
            ,'late_int'
            ,'expgamma'
            ,'txgamma'         # doesn't converge
            ,'zxgamma'
            ,'gamma_over_z'    # doesn't converge
            ,'zxxgamma'        # gamma forced positive in firstderivs
            ,'gammaxxz'        # gamma forced positive in firstderivs
            ,'rdecay_m'
            ,'rdecay_de'
            ,'rdecay_mxde'
            ,'rdecay'
            ,'interacting'
#            ,'LCDM'
             ]


for key in firstderivs_functions:
    if key:
        if key == 'stepfall':
            int_lim = [[-0.1, 0.1], [-0.1, 0.1], [-0.1, 0.1]]
            int_lim = [[-1, 1], [-1, 1], [-0.8, 1]]
            names = ['Mcorr','matter','radiation','a_ombar',
                     'v_in','w_in','x_in']
        elif key =='waterfall':
            int_lim = [[-1, 1], [-1, 1], [-1, 1],[-1, 1], [-1, 1]]
            names = ['Mcorr','matter','radiation','a_ombar','b_ombar','c_ombar',
                     'v_in','w_in','x_in','y_in','z_in']
        elif key == 'exotic':
            names = ['Mcorr','matter','radiation','gamma','zeta']
            int_lim = [[-2, 0.1],[-1.5, 2.5]]
        elif key == 'LCDM':
            int_lim = None
            names = ['Mcorr','matter']
        else:
            names = ['Mcorr','matter','gamma']
            if  key == 'late_intxde':
                int_lim = [[-2, 0.1]]
            elif key == 'heaviside_late_int':
                int_lim = [[-1.45, 0.1]]
            elif key == 'late_int':
                int_lim = [[-15, 0.1]]
            elif key == 'expgamma':
                int_lim = [[-0.1, 1.5]]
            elif key == 'txgamma':
                int_lim = [[-0.5, 0.1]]
            elif key == 'zxgamma':
                int_lim = [[-10, 0.1]]
            elif key == 'zxxgamma':
                int_lim = [[-0.1, 12]]
            elif key == 'gammaxxz':
                int_lim = [[-1, 1]]
            elif key == 'rdecay_m':
                int_lim = [[-3, 0]]
            elif key == 'rdecay':
                int_lim = [[-2, 0]]
            elif key == 'interacting':
                int_lim = [[-1.5, 0.1]]
            else:
                int_lim = [[-10,10]]

        if int_lim:
            fluid_number = len(names) - 1 - len(int_lim)
        else:
            fluid_number = len(names) - 1

        # Create a model object and a sampler.
        model = Model(names, int_lim, fluid_number)
        sampler = dnest4.DNest4Sampler(model,
                            backend=dnest4.backends.CSVBackend(".",sep=" "))

        if speed == 3: # LONG
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
            sortby = 'cumulative'
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

        if plot:
            hue = ['light red', 'berry', 'coral', 'amber', 'apple',
                        'aquamarine', 'raspberry', 'green blue', 'deep blue',
                        'emerald', 'blue violet', 'dark violet', 'yellow orange']
#            ndim = len(array[0,:])
            ndim = len(names)
            for i in range(ndim):
                name = names[i]
                plt.figure()
                plt.title(name)
                plt.hist(array[:,i], color='xkcd:'+hue[i])
                distribution = array[:,i]
                # Standard deviation and mean of the DNest distribution.
                DNest_distr[name+'_sd'] = np.std(distribution)
                DNest_distr[name+'_mean'] = np.mean(distribution)
                DNest_distr[name] = array[:,i]
        plt.show()

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
