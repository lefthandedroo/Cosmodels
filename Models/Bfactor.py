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
#from numba import jitclass, int32
import datasim
import results
import tools
#from scipy.special import erf

# from prior = 0, short = 1, medium = 2, long = 3
speed = 1
timed = False
plot = True
# Sigma of the noise on data.
sigma = 0.07

# Loading data:
#dataname = 'mag_z_LCDM_1000_sigma_'+str(sigma)
#mag, zpicks = results.load('./data', dataname)

dataname = './data/lcparam_full_long.txt'
pantheon = pd.read_csv(dataname, sep=" ")

# Reading each txt file column of interest as numpy.ndarray
mag = pantheon.mb.values
x1 = pantheon.x1.values
colour = pantheon.color.values
zpicks = pantheon.zhel.values

# Stacking them together and sorting by accending redshift.
data = np.stack((mag,x1,colour,zpicks), axis=0)
data.sort(axis=-1)

mag = data[0]
zpicks = data[3]
zpicks = zpicks.tolist()
data_dict = {'mag':mag, 'zpicks':zpicks}
ombar_names = ['matter', 'radiation', 'de', 'a', 'b', 'c', 'd', 'e', 'f']
int_names = ['p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']

class Model(object):
    """
    Specify the model in Python.
    """
    def __init__(self, int_lim=False, ombar_lim=False):
        """
        Parameter values *are not* stored inside the class
        """
        self.M_min = -20
        self.M_max = -18
        self.ombar_lim = ombar_lim
        if int_lim:
            for i in range(len(int_lim)):
                if i == 1:
                    self.v_min = int_lim[i][0]
                    self.v_max = int_lim[i][1]
                    if i == 2:
                        self.w_min = int_lim[i][0]
                        self.w_max = int_lim[i][1]
                        if i == 3:
                            self.x_min = int_lim[i][0]
                            self.x_max = int_lim[i][1]
                            if i == 3:
                                self.y_min = int_lim[i][0]
                                self.y_max = int_lim[i][1]
                                if i == 4:
                                    self.z_min = int_lim[i][0]
                                    self.z_max = int_lim[i][1]
#        pass

    def from_prior(self):
        """
        Unlike in C++, this must *return* a numpy array of parameters.
        """
        m = rng.rand()
        M = 1E3*rng.rand()
        M = dnest4.wrap(M, self.M_min, self.M_max)

        if self.ombar_lim:
            radiation = rng.rand()
            a_ombar = rng.rand()
            b_ombar = rng.rand()
            c_ombar = rng.rand()
                         
        if int_lim:
            for i in range (len(int_lim)):
                if i == 1:
                    g = 1E3*rng.rand()
                    g = dnest4.wrap(g, self.g_min, self.g_max)
                    return np.array([m, M, radiation, a_ombar, b_ombar, c_ombar, g])
                
                    if len(int_lim) == 2:
                        z = 1E3*rng.rand()
                        z = dnest4.wrap(z, self.z_min, self.z_max)
                        return np.array([m, M, g, z])
                    
                    elif len(int_lim) == 5:
                        z = 1E3*rng.rand()
                        z = dnest4.wrap(z, self.z_min, self.z_max)
                        return np.array([m, M, g, z])
                              
        return np.array([m, M])

    def perturb(self, params):
        """
        Unlike in C++, this takes a numpy array of parameters as input,
        and modifies it in-place. The return value is still logH.
        """
        logH = 0.0
        which = rng.randint(len(params))
        # Note the difference between dnest4.wrap in Python and
        # DNest4::wrap in C++. The former *returns* the wrapped value.
        
        if which == 0:
            params[which] += dnest4.randh()
            params[which] = dnest4.wrap(params[which], 0.0, 1.0)
        elif which == 1:
            params[which] += dnest4.randh()
            params[which] = dnest4.wrap(params[which], self.M_min, self.M_max)            
        elif which == 2:
            params[which] += dnest4.randh()
            params[which] = dnest4.wrap(params[which], self.g_min, self.g_max)
        elif which == 3:
            params[which] += dnest4.randh()
            params[which] = dnest4.wrap(params[which], self.z_min, self.z_max)
        return logH

    def log_likelihood(self, params):
        """
        Gaussian sampling distribution.
        """
        if len(params) == 11:
            m, M, r, a, b, c, v, w, x, y, z = params
            theta = [{'matter':m},{'Mcorr':M},{'radiation':r},
                          {'a_ombar':a},{'b_ombar':b},{'c_ombar':c},
                          {'v_in':v},{'w_in':w},{'x_in':x},{'y_in':y},
                          {'z_in':z}]
        elif len(params) == 4:
            m, M, g, z = params
            theta = [{'matter':m}, {'Mcorr':M}, {'gamma':g}, {'zeta':z}]
        elif len(params) == 3:
            m, M, g = params
            theta = [{'matter':m}, {'Mcorr':M}, {'gamma':g}]
        else:
            m, M = params
            theta = [{'matter':m}, {'Mcorr':M}]
        
        model = datasim.magn(theta, data_dict, key)
        
        var = sigma**2.0
        return -0.5*np.sum((mag-model)**2.0 /var +0.5*np.log(2.0*np.pi*var))
    
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


firstderivs_functions = [None
#            ,'waterfall'
            ,'exotic'
#            ,'late_intxde'
#            ,'heaviside_late_int'
#            ,'late_int'
#            ,'expgamma'
#            ,'txgamma'         # doesn't converge
#            ,'zxgamma'
#            ,'gamma_over_z'    # doesn't converge
#            ,'zxxgamma'        # gamma forced positive in firstderivs
#            ,'gammaxxz'        # gamma forced positive in firstderivs
#            ,'rdecay_m'
#            ,'rdecay_de'
#            ,'rdecay_mxde'
#            ,'rdecay'               
#            ,'interacting'
#            ,'LCDM'
             ]

for key in firstderivs_functions:
    if key:
        if key =='waterfall':
            params_dic = [{'matter':0.3},{'Mcorr':-19.3},{'radiation':0.025},
                          {'a_ombar':0.1},{'b_ombar':0.1},{'c_ombar':0.1},
                          {'v_in':0.0},{'w_in':0.0},{'x_in':0.0},
                          {'y_in':0.0},{'z_in':0.0}]
        elif key == 'exotic':
            int_lim = [[-2, 0.1],[-1.5, 2.5]]
            params_dic = [{'matter':0.3},{'Mcorr':-19.3},
                          {'gamma':0},{'zeta':0}]
        elif key == 'late_intxde':
            int_lim = [[-2, 0.1]]
            params_dic = [{'matter':0.3},{'Mcorr':-19.3},{'gamma':0}]
        elif key == 'heaviside_late_int':
            int_lim = [[-1.45, 0.1]]
            params_dic = [{'matter':0.3},{'Mcorr':-19.3},{'gamma':0}]
        elif key == 'late_int':
            int_lim = [[-15, 0.1]]
            params_dic = [{'matter':0.3},{'Mcorr':-19.3},{'gamma':0}]
        elif key == 'expgamma':
            int_lim = [[-0.1, 1.5]]
            params_dic = [{'matter':0.3},{'Mcorr':-19.3},{'gamma':0}]
        elif key == 'txgamma':
            int_lim = [[-0.5, 0.1]]
            params_dic = [{'matter':0.3},{'Mcorr':-19.3},{'gamma':0}]
        elif key == 'zxgamma':
            int_lim = [[-10, 0.1]]
            params_dic = [{'matter':0.3},{'Mcorr':-19.3},{'gamma':0}]
        elif key == 'zxxgamma':
            int_lim = [[-0.1, 12]]
            params_dic = [{'matter':0.3},{'Mcorr':-19.3},{'gamma':0}]
        elif key == 'gammaxxz':
            int_lim = [[-1, 1]]
            params_dic = [{'matter':0.3},{'Mcorr':-19.3},{'gamma':0}]
        elif key == 'rdecay_m':
            int_lim = [[-3, 0]]
            params_dic = [{'matter':0.3},{'Mcorr':-19.3},{'gamma':0}]
        elif key == 'rdecay':
            int_lim = [[-2, 0]]
            params_dic = [{'matter':0.3},{'Mcorr':-19.3},{'gamma':0}]
        elif key == 'interacting':
            int_lim = [[-1.5, 0.1]]
            params_dic = [{'matter':0.3},{'Mcorr':-19.3},{'gamma':0}]
        elif key == 'LCDM':
            int_lim = None
            params_dic = [{'matter':0.3},{'Mcorr':-19.3}]
        else:
            int_lim = [[-10,10]]
        
        # Create a model object and a sampler
        model = Model(int_lim, ombar_lim)
        sampler = dnest4.DNest4Sampler(model,
                                       backend=dnest4.backends.CSVBackend(".",
                                                                      sep=" "))
        
        if speed == 3:   
            # LONG Set up the sampler. The first argument is max_num_levels
            gen = sampler.sample(max_num_levels=30, num_steps=1000, 
                                 new_level_interval=10000, num_per_step=10000, 
                                 thread_steps=100, num_particles=5, 
                                 lam=10, beta=100, seed=1234)
        elif speed == 2:
            # MEDIUM num_per_step can be down to a few thousand 
            gen = sampler.sample(max_num_levels=30, num_steps=1000, 
                                 new_level_interval=1000, num_per_step=1000, 
                                  thread_steps=100, num_particles=5, 
                                  lam=10, beta=100, seed=1234)
        elif speed == 1:
            # SHORT
            gen = sampler.sample(max_num_levels=30, num_steps=100, 
                                 new_level_interval=100, num_per_step=100, 
                                 thread_steps=10, num_particles=5, 
                                 lam=10, beta=100, seed=1234)
        elif speed == 0:
            # SHORT, sampling from prior
            gen = sampler.sample(max_num_levels=1, num_steps=1000, 
                                 new_level_interval=100, num_per_step=100, 
                                 thread_steps=10, num_particles=5, 
                                 lam=10, beta=100, seed=1234)    

        if timed:       
            import cProfile, pstats, io
            pr = cProfile.Profile()
            pr.enable()
        
        ti = time.time()
        
        # Do the sampling (one iteration here = one particle save)
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
            ndim = len(array[0,:])            
            for i in range(ndim):
                for key in params_dic[i]:
                    plt.figure()
                    plt.title(key)
                    plt.hist(array[:,i], color='xkcd:'+hue[i])
                    distribution = array[:,i]
                    # Standard deviation and mean of the DNest distribution.
                    DNest_distr[key+'_sd'] = np.std(distribution)
                    DNest_distr[key+'_mean'] = np.mean(distribution)
                    DNest_distr[key] = array[:,i]
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
            