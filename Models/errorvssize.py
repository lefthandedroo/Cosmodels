#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 18:31:51 2018

@author: BallBlueMeercat

Runs through given sigma (error on data) and npoints (dataset size), use pre-made artificial mag or create new mag using pre-made zpicks. Saves sampler.p for each run.
"""
import matplotlib.pyplot as plt
import numpy as np
import time
import datasim
from pathlib import Path
import tools
import os.path
import stats
import pickle

import sys
import emcee
from emcee.utils import MPIPool

try:
    pool = MPIPool()
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
except:
    print('pool = None')
    pool = None

# Script timer.
timet0 = time.time()

nsteps = 100#100000   # Number of emcee steps per?
mu = 0.0        # Mean of noise added to LCDM to simulate data.

#Model, errors on data and dataset sizes to iterate through:
test_keys = [None
#            ,'rainbow'
#            ,'kanangra'
#            ,'waterfall'
#            ,'stepfall'
#            ,'exotic'
#            ,'late_intxde'
#            ,'heaviside_late_int'
#            ,'heaviside_sudden'
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
            ,'LCDM'
#            ,'rLCDM'
            ]
min_z = 0.01012
max_z = 2.26 # highest expected redshift for a type Ia supernova from pantheon

sigma_options = None, 0.001 #0.14, 0.2 #0.0001, 0.005, 0.007
npoints_options = None,  10480, #1048, 10480, 104800, 1048000, 10480000
run = 0

for key in test_keys:
    if key:
        names, values = tools.names_values(key)
        # Folder for saving output.
        directory = f'{int(time.time())}_{key}'
        # Relative path of output folder.
        save_path = os.path.join('results_error_vs_data',directory)
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        for sigma in sigma_options:
            if sigma:
                pass
            else:
                continue
            for npoints in npoints_options:
                if npoints:
                    pass
                else:
                    continue
                run += 1
                data_path = f'data/{npoints}_{max_z}_sigma_{sigma}.p'
                my_file = Path(data_path)
                if my_file.is_file():
                    with open(data_path,'rb') as rfp: zpicks, nmag = pickle.load(rfp)
                else:
                    print(f'failed to get zpicks, nmag from {data_path}')
                    # Generating redshifts.
                    zpicks = np.random.uniform(low=min_z, high=max_z, size=(npoints,))
                    zpicks = np.sort(zpicks, axis=None)
                    if zpicks[-1] != max_z:
                        zpicks[-1] = max_z
                    data_dic = {'zpicks':zpicks}
                    # Generating LCDM mag and da.
                    mag, da = datasim.magn(['Mcorr', 'matter'], np.array([-19.3, 0.3]), data_dic, 'LCDM')
                    # Adding noise to LCDM mag.
                    nmag = datasim.gnoise(mag, mu, sigma)

                    plt.figure()
                    plt.title(f'Artificial data N={len(zpicks)}, $\sigma$={sigma}')
                    plt.scatter(zpicks, nmag)
                    plt.show()

                    data = zpicks, nmag
                    pickle.dump(data, open(data_path, 'wb'))

                data_dic = {'mag':nmag, 'zpicks':zpicks}

                print(f'--- {key} --------- run number {run}')
                propert, sampler = stats.stats(names, values, data_dic,
                                               sigma, nsteps, save_path,
                                               key, pool=None, plot=False)
                try:
                    pool.close()
                except:
                    pass
#                output = propert, sampler
                output_path = os.path.join(save_path, f'{key}_sigma{sigma}_npoints{npoints}.txt')
#                pickle.dump(output, open(output_path, 'wb'))
                output = sampler.flatchain[:, :]
                f = open(output_path, 'w+')
                f.write(str(output))
                f.close()


print('directory:',directory)

# Time taken by script.
timet1=time.time()
tools.timer('errorvsdatasize', timet0, timet1)

#plots.precise_runs(key, names, values, 9, 1.2)