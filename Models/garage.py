#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""
import pandas as pd
xl = pd.ExcelFile("hpc_times.xlsx")
df = xl.parse("Sheet1")
df_cores = df.cores.values
df_time = df.time.values
print(df_cores)
print(df_time)
#df_time_s = (df_time.dt.hour*60+df_time.dt.minute)*60 + df_time.dt.second
#print(df_time_s)




#import matplotlib.pyplot as plt
#import matplotlib as mpl
#mpl.style.use('default') # has to be switched on to set figure size
#mpl.style.use('fivethirtyeight')
#plt.rcParams['axes.facecolor'] = 'white'
#plt.rcParams['figure.facecolor'] = 'white'
#plt.rcParams['grid.color'] = 'white'
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