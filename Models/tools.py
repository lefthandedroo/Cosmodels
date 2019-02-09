#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 17:06:45 2018

@author: BallBlueMeercat
"""
import numpy as np
import time
import os.path

def flist(start, stop, step):
    """
    Takes in:
        start, stop, step = integers or floats
    Returns:
        zlist = list start to stop with step as increment
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


def timer(string,i,f):
    """
    Takes in:
        i = starting time;
        f = finishing time.
    Returns: Time taken in full minutes and seconds.
    """
    sec = f - i     # Total time to run.
    mins, sec= divmod(sec, 60.0)
    time = string+' time:  '+str(int(mins))+'min '+str(int(sec))+'s'
    print(time)

    return time


def runcount(test_keys, sigma, sigma_max, sigma_step,
              npoints_min, npoints_max, npoints_step):
    """
    Returns number of interations that evaluator will run.
    """
    run = 1
    for key in test_keys:
        if key:
            while sigma < sigma_max:
                npoints = npoints_min
                while npoints < npoints_max:
                    npoints += npoints_step
                    run += 1
                sigma += sigma_step
    return run

def path():
    # Folder for saving output.
    directory = 'run'+str(int(time.time()))
    # Relative path of output folder.
    save_path = './'+directory
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    return save_path, directory

def names_values(test_key):
    if test_key =='rainbow':
        names = ['Mcorr',
                 'm_ombar', 'r_ombar', 'a_ombar', 'b_ombar', 'c_ombar',
                 'd_ombar', 'e_ombar', 'f_ombar', 'g_ombar', 'h_ombar',
                 'i_ombar',
                 'p_in', 'q_in', 'r_in', 's_in', 't_in', 'u_in',
                 'v_in', 'w_in', 'x_in', 'y_in', 'z_in']
        values = np.array([-19.3,
                           0.3, 0.025, 0.01, 0.01, 0.01, 0.01,
                           0.01, 0.01, 0.01, 0.01, 0.01,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0])
    elif test_key == 'kanangra':
        names = ['Mcorr',
                 'm_ombar', 'r_ombar', 'a_ombar', 'b_ombar', 'c_ombar',
                 'd_ombar', 'e_ombar',
                 't_in', 'u_in', 'v_in', 'w_in', 'x_in', 'y_in', 'z_in']
        values = np.array([-19.3,
                           0.3, 0.025, 0.01, 0.01, 0.01,
                           0.01,0.01,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    elif test_key == 'waterfall':
        names = ['Mcorr',
                 'm_ombar', 'r_ombar', 'a_ombar', 'b_ombar', 'c_ombar',
                 'v_in', 'w_in', 'x_in', 'y_in', 'z_in']
        values = np.array([-19.3,
                           0.3, 0.025, 0.1, 0.1, 0.1,
                           0.0, 0.0, 0.0, 0.0, 0.0])
    elif test_key == 'stepfall':
        names = ['Mcorr', 'm_ombar', 'r_ombar', 'a_ombar',
                 'v_in', 'w_in', 'x_in']
        values = np.array([-19.3, 0.3, 0.025, 0.1, 0.0, 0.0, 0.0])
    elif test_key == 'exotic':
        names = ['Mcorr', 'm_ombar', 'r_ombar', 'gamma', 'zeta']
        values = np.array([-19.3, 0.3, 0.025, 0.0, 0.0])
    elif test_key == 'rLCDM':
        names = ['Mcorr', 'm_ombar', 'r_ombar']
        values = np.array([-19.3, 0.3, 0.025])
    elif test_key == 'LCDM':
        names = ['Mcorr', 'm_ombar']
        values = np.array([-19.3, 0.3])
    else:
        names = ['Mcorr', 'm_ombar','gamma']
        values = np.array([-19.3, 0.3, -0.0])
    return names, values