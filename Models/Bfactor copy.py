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
import pickle
#from numba import jitclass, int32
import datasim
import results
import tools
#from scipy.special import erf

# from prior = 0, short = 1, medium = 2, long = 3
speed = 1

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
x1 = data[1]
colour = data[2]
zpicks = data[3]
zpicks = zpicks.tolist()
data_dict = {'mag':mag, 'x1':x1, 'colour':colour, 'zpicks':zpicks}

class Model(object):
    """
    Specify the model in Python.
    """
    def __init__(self, g_lim=None):
        """
        Parameter values *are not* stored inside the class
        """
        self.M_min = -20
        self.M_max = -18
        self.a_min = -20
        self.a_max = 20
        self.b_min = -20
        self.b_max = 20
        
        if g_lim != None:
            self.g_min = g_lim[0]
            self.g_max = g_lim[1]
#        pass

    def from_prior(self):
        """
        Unlike in C++, this must *return* a numpy array of parameters.
        """
        m = rng.rand()
        M = 1E3*rng.rand()
        M = dnest4.wrap(M, self.M_min, self.M_max)
        a = 1E3*rng.rand()
        a = dnest4.wrap(a, self.a_min, self.a_max)
        b = 1E3*rng.rand()
        b = dnest4.wrap(b, self.b_min, self.b_max)
        if g_lim != None:
            g = 1E3*rng.rand()
            g = dnest4.wrap(g, self.g_min, self.g_max)
            return np.array([m, M, a, b, g])
        return np.array([m, M, a, b])

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
            params[which] = dnest4.wrap(params[which], self.a_min, self.a_max)
        elif which == 3:
            params[which] += dnest4.randh()
            params[which] = dnest4.wrap(params[which], self.b_min, self.b_max)            
        elif which == 4:
            params[which] += dnest4.randh()
            params[which] = dnest4.wrap(params[which], self.g_min, self.g_max)

        return logH

    def log_likelihood(self, params):
        """
        Gaussian sampling distribution.
        """
        if len(params) > 4:
            m, M, a, b, g = params
            theta = {'m':m, 'M':M, 'a':a, 'b':b, 'gamma':g}
        else:
            m, M, a, b = params
            theta = {'m':m, 'M':M, 'a':a, 'b':b}
        
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
            ,'exotic'
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
            ,'LCDM'
             ]

for key in firstderivs_functions:
     if key:   
        if key == 'exotic':
            g_lim = [-1.5, 0.1]
            
        elif key == 'late_intxde':
            g_lim = [-2, 0.1]
            
        elif key == 'heaviside_late_int':
            g_lim = [-1.45, 0.1]
            
        elif key == 'late_int':
            g_lim = [-15, 0.1]
            
        elif key == 'expgamma':
            g_lim = [-0.1, 1.5]
            
        elif key == 'txgamma':
            g_lim = [-0.5, 0.1]
            
        elif key == 'zxgamma':
            g_lim = [-10, 0.1]
            
        elif key == 'zxxgamma':
            g_lim = [-0.1, 12]
            
        elif key == 'gammaxxz':
            g_lim = [-1, 1]
            
        elif key == 'rdecay_m':
            g_lim = [-3, 0]
            
        elif key == 'rdecay':
            g_lim = [-2, 0]
            
        elif key == 'interacting':
            g_lim = [-1.5, 0.1]
            
        elif key == 'LCDM':
            g_lim = None
        
        else:
            g_lim = [-10,10]
        
        # Create a model object and a sampler
        model = Model(g_lim)
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
        
    #    import cProfile, pstats, io
    #    pr = cProfile.Profile()
    #    pr.enable()
        
        ti = time.time()
        # Do the sampling (one iteration here = one particle save)
        for i, sample in enumerate(gen):
    #        print("# Saved {k} particles.".format(k=(i+1)))
            pass
        tf = time.time()
        
    #    pr.disable()
    #    s = io.StringIO()
    #    sortby = 'cumulative'
    #    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    #    ps.print_stats()
    #    print (s.getvalue())
        
        dnest_time = tools.timer('Bfactor', ti, tf)
        
        print('testing =',key)
        print('data =', dataname)
        print('sigma =', sigma)
           
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
            
        else:
            # Histogram of parameters found by DNest4.
            array = np.loadtxt('sample.txt')
            import matplotlib.pyplot as plt
            if key == 'LCDM':
                plt.figure()
                plt.title('matter')
                plt.hist(array[:,0])
                plt.show()
            
                plt.figure()
                plt.title('M_b')
                plt.hist(array[:,1])
                plt.show()   
                
                plt.figure()
                plt.title('alpha')
                plt.hist(array[:,2])
                plt.show()
                
                plt.figure()
                plt.title('beta')
                plt.hist(array[:,3])
                plt.show()            
            else:
                plt.figure()
                plt.title('matter')
                plt.hist(array[:,0])
                plt.show()
            
                plt.figure()
                plt.title('M_b')
                plt.hist(array[:,1])
                plt.show()   
                
                plt.figure()
                plt.title('alpha')
                plt.hist(array[:,2])
                plt.show()
                
                plt.figure()
                plt.title('beta')
                plt.hist(array[:,3])
                plt.show()
                
                plt.figure()
                plt.title('gamma')
                plt.hist(array[:,4])
                plt.show()




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