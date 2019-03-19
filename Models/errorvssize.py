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
from emcee.utils import MPIPool

import argparse

parser = argparse.ArgumentParser(description="This program does stuff", usage="name_of_script --nsteps 2000 --xwalkers 10 --datapoints 10480", epilog="this is (optional) text added to the end of the help")

# argument: files to process
# short flag (one character) followed by verbose flag
parser.add_argument("-n", "--nsteps", type=int,
                    help="number of emcee steps per walker")
parser.add_argument("-x", "--xwalkers", type=int,
                    help="number to multiply ndim by to get number of walkers")
parser.add_argument("-d", "--datapoints", type=int,
                    help="SN Ia data set size")
parser.add_argument("-f", "--filename",
                    help="name of a file to save printing of progress to")
args = parser.parse_args()

timed = False
plot = False

node = 'master'
try:
    pool = MPIPool()
    if not pool.is_master():
        node = 'worker'
        print('worker')
        pool.wait()
        sys.exit(0)
except Exception as e:
    print('exception is = ',e)
    print('pool = None')
    pool = None
    plot = True

if timed:
    import cProfile, pstats, io
    pr = cProfile.Profile()
    pr.enable()

filename = args.filename
if filename:
        pass
else:
    filename = 'print.txt'

f = open(filename,"w+")
f.write('pool = '+str(pool)+'\n')
f.close()

# Script timer.
timet0 = time.time()

nsteps = args.nsteps
xwalkers = args.xwalkers
datapoints = args.datapoints
if nsteps:
    pass
else:
    nsteps = 2000  # Number of emcee steps per walker
print('nsteps',nsteps)

if xwalkers:
    pass
else:
   xwalkers = 10
print('xwalkers',xwalkers)

if datapoints:
    pass
else:
   datapoints = 1048
print('datapoints',datapoints)

mu = 0.0        # Mean of noise added to LCDM to simulate data.

#Model, errors on data and dataset sizes to iterate through:
test_keys = [None
#            ,'rainbow'
#             ,'niagara'
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
npoints_options = None,  datapoints, #1048, 10480, 104800, 1048000, 10480000
run = 0

for key in test_keys:
    if key:
        names, values = tools.names_values(key)
        # Folder for saving output.
        directory = f'{int(time.time())}_{key}'
        # Relative path of output folder.
        save_path = os.path.join('results_error_vs_data',directory)
        if not os.path.exists(save_path):
            os.makedirs(save_path, exist_ok=True)
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
                    plt.show(block=False)

                    data = zpicks, nmag
                    pickle.dump(data, open(data_path, 'wb'))

                data_dic = {'mag':nmag, 'zpicks':zpicks}

                f=open(filename, "a+")
                f.write('got to before stats'+'\n')
                f.close()

                print(f'--- {key} --------- run number {run}')
                propert, sampler = stats.stats(names, values, data_dic,
                    sigma, nsteps, save_path, key, xwalkers=xwalkers, pool=None,
                    plot=plot, filename=filename)

                f=open(filename, "a+")
                f.write('after stats'+'\n')
                f.close()

                try:
                    pool.close()
                except:
                    pass
#                output = propert, sampler
                output_path = os.path.join(save_path, f'{key}_sigma{sigma}_npoints{npoints}.txt')
#                pickle.dump(output, open(output_path, 'wb'))
                print('-------------------------------- f.writing')
                print(node,' node')
                output = sampler.flatchain[:, :]
                np.savetxt(output_path, output)

f=open(filename, "a+")
f.write('after the loop')
f.close()

print('directory:',directory)

# Time taken by script.
timet1=time.time()
tools.timer('errorvsdatasize', timet0, timet1)
print(' ')

if timed:
    pr.disable()
    s = io.StringIO()
    sortby = 'tottime'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print (s.getvalue())
#plots.precise_runs(key, names, values, 9, 1.2)