#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('default') # has to be switched on to set figure size
mpl.style.use('fivethirtyeight')
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['grid.color'] = 'white'

steps=[25, 250, 2500, 25000, 250000]
times = [2, 4, 31, 287, 2760]

fig, ax = plt.subplots()
ax.scatter(times, steps, color="C{}".format(6), s=70, zorder=2)
for i, txt in enumerate(steps):
    txt = str(txt)
    ax.annotate(txt, (times[i], steps[i]))
plt.plot(times, steps, zorder=1, color="C{}".format(7))
plt.xscale('log')
plt.ylabel('MCMC steps')
plt.yticks([])
plt.xlabel('Time (s)')
plt.show()
#
#cores = [
#1,1,1,1,1,2,2,2,2,2,2,2,2,4,4,4,4,7,7,7,7,7,7,10,10,10,10,13,13,13,13,15,15,15,15,18,18,18,18,19,19,19,19,19,20,20,20,20
#]
#
#times =[
#'00:03:08','00:03:27','00:03:17','00:02:32','00:03:19','00:09:37',None,None,None,'00:07:49',None,None,'00:06:42','00:31:48','00:34:07','00:32:06','00:28:57','00:55:02','00:56:28','00:53:30','00:52:05',None,None,'00:58:47','00:19:20','00:59:00','00:58:55','01:37:53','01:02:26','01:38:16','01:00:29','01:11:35','01:04:18','01:03:42','01:03:49','01:24:09','01:51:12','01:51:48','00:58:29','00:28:37','00:38:28',None,'01:07:49','00:41:56','00:47:51','00:49:21','00:59:21','00:48:16'
#]
#
#cores_dic = {}
#key = 1
#time_list = []
#for i, core in enumerate(cores):
#    if core != key:
#        cores_dic[key] = time_list
#        key = core
#        time_list = []
#    str_time = times[i]
#    if str_time:
#        min_time = int(str_time[:2])*60 +int(str_time[3:5]) +float(str_time[6:])/60
#        time_list.append(min_time)
#
#error_list = []
#
#for core in cores_dic:
#    time_list = cores_dic[core]
#    mean_time = sum(time_list)/len(time_list)
#    time_min = min(time_list)
#    time_max = max(time_list)
#    cores_dic[core] = mean_time
#    error_list.append((time_max-time_min)/2)
#
#times = cores_dic.values()
#cores = cores_dic.keys()
#
#ax = plt.figure().gca()
#plt.scatter(cores, times)
#plt.errorbar(cores, times, yerr=error_list, color="C{}".format(5), fmt='.', alpha=0.3)
#ax.xaxis.set_major_locator(MaxNLocator(integer=True))
#plt.ylabel('Minutes')
#plt.xlabel('Cores')
#
## Pantheon data:
#import pandas as pd
#pantheon = pd.read_csv('./data/lcparam_full_long.txt', sep=" ")
#pantheon.set_index('name', inplace=True)
#pantheon.sort_values('zhel', inplace=True)
#mag = pantheon.mb.values
#zpicks = pantheon.zhel.values
#dmag = pantheon.dmb.values
#
#plt.figure()
#data = plt.errorbar(zpicks, mag, yerr=dmag, fmt='.', lw=1, alpha=0.5)
#plt.ylabel('Mag')
#plt.xlabel(r'$z$')


#run = 0
#from pathlib import Path
#import matplotlib.pyplot as plt
#
#sigma_options = None, 0.001#, 0.01, 0.07, 0.14, 0.2
#npoints_options = None, 104800 #1048, 10480, 104800, 1048000
#for sigma in sigma_options:
#    if sigma:
#        pass
#    else:
#        continue
#    for npoints in npoints_options:
#        if npoints:
#            pass
#        else:
#            continue
#        data_path = f'data/{npoints}_3.0_sigma_{sigma}.p'
#        my_file = Path(data_path)
#        if my_file.is_file():
#            run += 1
#            with open(data_path,'rb') as rfp: zpicks, nmag = pickle.load(rfp)
#            plt.figure()
#            plt.title(data_path)
#            plt.scatter(zpicks, nmag)
#            plt.show()
#        else:
#            print(data_path, "didn't open")
#            pass
#
#new_zpicks = np.linspace(0.001, 2.26, num=104800, endpoint=True)


#all_zpicks = zpicks
#print('zpicks',zpicks)
#print('zpicks[-1]',zpicks[-1])
#if len(zpicks) > 2000:
#    zpicks_cut = np.linspace(0, zpicks[-1], num=1000, endpoint=True)
#    print('zpicks_cut[-1]',zpicks_cut[-1])

#a = np.arange(5)
#hist, bin_edges = np.histogram(a, density=True)

#def wrap(x, a, b):
#    assert b > a
#    return (x - a)%(b - a) + a
#
#
#N = 3
#
#ombar_wm = 0.3    # starting from 0.3
#ombar_wr = 0.025    # starting from 0.025
#
#density_left = 1.0 - ombar_wm - ombar_wr       # e_total/e_crit = 1
#
#lower = 0.0
#upper = density_left/ N
#
#ombar_w01 = rng.rand()      # random number between 0 and 1
#ombar_w01 = wrap(ombar_w01, lower, upper)
#ombar_w05 = rng.rand()
#ombar_w05 = wrap(ombar_w05, lower, upper)
#ombar_w08 = rng.rand()
#ombar_w08 = wrap(ombar_w08, lower, upper)
#ombar_wde = density_left - ombar_w01 - ombar_w05 - ombar_w08
#
##print('ombar_m = ',ombar_wm)
##print('ombar_r = ',ombar_wr)
##print('ombar_w(-0.1) = ',ombar_w01)
##print('ombar_w(-0.5) = ',ombar_w05)
##print('ombar_w(-0.8) = ',ombar_w08)
##print('ombar_de = ',ombar_wde)
#
#total = ombar_wm +ombar_wr +ombar_w01 +ombar_w05 +ombar_w08 +ombar_wde
#print('total =',total)
#
#import datasim
#import numpy as np
#try:
#    with open('zpicks_1089.p','rb') as rfp: zpicks = pickle.load(rfp)
#except:
#    print("zpicks_1089.p didnt't open")
#print(len(zpicks))
#data_dic = {'zpicks':zpicks}
#mag, da = datasim.magn(['Mcorr', 'matter'], np.array([-19.3, 0.3]), data_dic, 'LCDM', plot_key=False)
#plt.figure()
#plt.title(f'Artificial data, zpicks_1089')
#plt.scatter(zpicks, mag)
#


#try:
#    with open('data/zpicks_1048_3.p','rb') as rfp: all_zpicks = pickle.load(rfp)
#except:
#    print("zpicks_1048_3.p didnt't open")
#
#data_dic = {'zpicks':zpicks}
#mag, da = datasim.magn(['Mcorr', 'matter'], np.array([-19.3, 0.3]), data_dic, 'LCDM', plot_key=False)
#plt.figure()
#plt.title(f'Artificial data, zpicks_10890')
#plt.scatter(zpicks, mag)
#
#plt.show()
#
## Generating and saving redshifts.
#zpicks = np.random.uniform(low=0.0001, high=3, size=(1048,))
#zpicks = np.sort(zpicks, axis=None)
#zpicks[-1] = 3
#pickle.dump(zpicks, open(f'zpicks_{len(zpicks)}_{zpicks[-1]}.p', 'wb'))