#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 18:06:28 2018

@author: usyd

This script collects propert, sampler for the relevant model in
results_erro_vs_data and integrates the equations describing the model with
every set of parameters from emcee steps stored inside the sampler. Integration
is over a single z=1089 (z=0 is incerted as initial condition in zodesolve.py).
This integration results in mag and angular diameter distance at recombination.
The distribution of angular diameter distances is saved in current directory.

"""
from pathlib import Path
import pickle
import numpy as np
import datasim


zpicks = np.array([1089])
data_dic = {'zpicks':zpicks}

timed = False

if timed:
    import cProfile, pstats, io
    pr = cProfile.Profile()
    pr.enable()

models = None, 'stepfall'
noise_options = 0.001, None      #0.001, 0.01, 0.07, 0.14, 0.2
npoints_options = None, 1048000  #1048, 10480, 104800, 1048000

for npoints in npoints_options:
    if not npoints: # skipping None to allow iteratng over a list of 1
        continue
    for noise in noise_options:
        if not noise:
            continue
        for test_key in models:
            if test_key:
                file_path = f'results_error_vs_data/{test_key}/{test_key}_sigma{noise}_npoints{npoints}.p'
                my_file = Path(file_path)
                if my_file.is_file():
                    # results of emcee runs
                    with open(file_path,'rb') as rfp: propert, sampler = pickle.load(rfp)
                    file_open_indicator = 1
                    print(f'{file_path} opened successfully')
                else:
                    print(f"Couldn't open {file_path}")
                    file_open_indicator = 0
                # assert ensures only sampler from intended file_path is used
                assert file_open_indicator > 0, f"{file_path} didn't have a sampler"

                da_distrib = []
                for i in range(0, len(sampler.flatchain)):
                    values = sampler.flatchain[i, :]
                    # not parsing real names as not needed
                    names = np.zeros(len(values)).tolist()
                    mag, da = datasim.magn(names, values, data_dic,
                                           test_key, plot_key=False)
                    da_distrib.append(da[-1]) # adding da at z[-1]=1089
                print(f'last da = {da[-1]}, z = {zpicks}')
                filename = f'da_distribs/da_distrib_{test_key}_{noise}_{npoints}.p'
                pickle.dump(da_distrib, open(filename, 'wb'))

if timed:
    pr.disable()
    s = io.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print (s.getvalue())