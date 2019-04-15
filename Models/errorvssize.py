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

import sys
from emcee.utils import MPIPool

import argparse

parser = argparse.ArgumentParser(description="This program runs emcee on syntheitc SN Ia dataset of specified size", usage="name_of_script --model LCDM --nsteps 2000 --xwalkers 10 --datapoints 10480 --filename output.txt", epilog="this is (optional) text added to the end of the help")

# argument: files to process
# short flag (one character) followed by verbose flag
parser.add_argument("-m", "--model",
                    help="model to be fitted to data")
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

node= 'master'
try:
    pool = MPIPool()
    if not pool.is_master():
        node = 'worker'
        print(node)
        pool.wait()
        sys.exit(0)
#except ValueError as e:
except Exception as e:
    print('exception is = ',e)
    pool = None
    plot = False
print('pool =', pool)
print(node)

if timed:
    import cProfile, pstats, io
    pr = cProfile.Profile()
    pr.enable()

filename = args.filename
if not filename:
    filename = 'print.txt'

f = open(filename,"w+")
f.write('pool = '+str(pool)+'\n')

# Script timer.
timet0 = time.time()
args_key = args.model
nsteps = args.nsteps
xwalkers = args.xwalkers
datapoints = args.datapoints
if not args_key:
    args_key = 'LCDM'
if not nsteps:
    nsteps = 1000  # Number of emcee steps per walker
if not xwalkers:
    xwalkers = 10
if not datapoints:
    datapoints = 10480
print('nsteps',nsteps)
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
             ,args_key
#            ,'LCDM'
#            ,'rLCDM'
            ]

min_z = 0.01012
max_z = 2.26 # highest expected redshift for a type Ia supernova from pantheon

sigma_options = None, 0.001 #0.14, 0.2 #0.0001, 0.005, 0.007
npoints_options = None,  datapoints, #1048, 10480, 104800, 1048000, 10480000
#sigma_options = None, 0.001, #0.3, 0.15, 0.07, 0.03, 0.015, 0.007 #0.001
#npoints_options = None, 100, 1000, 10000, 30000, 50000, 70000, 100000, 300000, 500000, 700000, 1000000 #datapoints

run = 0

for key in test_keys:
    if not key:
        continue
    names, values = tools.names_values(key)
    # Folder for saving output.
    directory = f'{int(time.time())}_{key}'
    # Relative path of output folder.
    save_path = os.path.join('results_error_vs_data',directory)
    if not os.path.exists(save_path):
        os.makedirs(save_path, exist_ok=True)
    for sigma in sigma_options:
        if not sigma:
            continue
        for npoints in npoints_options:
            if not npoints:
                continue
            run += 1
            data_path = f'data/{npoints}_{max_z}_sigma_{sigma}.txt'
            my_file = Path(data_path)
            if my_file.is_file():
                data = np.loadtxt(data_path)
                zpicks, nmag = np.hsplit(data, 2)
                plt.figure()
                plt.scatter(zpicks, nmag)
                plt.title('loaded data')
                plt.show(block=False)
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

                data = np.column_stack((zpicks, nmag))
                np.savetxt(data_path, data)
            nmag = nmag.flatten()
            zpicks = zpicks.flatten()
            data_dic = {'mag':nmag, 'zpicks':zpicks}

            f = open(filename,"a+")
            f.write('got to before stats'+'\n')

            print(f'--- {key} --------- run number {run}')
            propert, sampler = stats.stats(names, values, data_dic,
                sigma, nsteps, save_path, key, xwalkers=xwalkers, pool=None,
                plot=plot, filename=filename)

            f.write('after stats'+'\n')

            try:
                pool.close()
            except:
                pass
#            import pickle
#            print('saving sampler with pickle as .p')
#            output = propert, sampler
#            output_path = os.path.join(save_path, f'{key}_sigma{sigma}_npoints{npoints}.p')
#            pickle.dump(output, open(output_path, 'wb'))
            output_path = os.path.join(save_path, f'{key}_sigma{sigma}_npoints{npoints}.txt')
            print('----------------------------- np.savetxt saving output')
            output = sampler.flatchain[:, :]
            np.savetxt(output_path, output)

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