#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from results import load
firstderivs_key = 'exotic'
params_dic = [{'matter':0.3},{'Mcorr':-19.3}, {'gamma':0.0}, {'zeta':0.0}]
p = 10
x = 2

directory = os.path.join('./results_error_vs_data/'+firstderivs_key)

sd_list = []
mean_list = []
vc_list = []
sigma_list = []
npoints_list = []

folders = []
for d in os.walk(directory):
    folders.append(d[0])
folders.pop(0)
    
for folder in folders:
    sd_list += load(folder, 'sd_list.p')
    mean_list += load(folder, 'mean_list.p')
    vc_list += load(folder, 'vc_list.p')
    sigma_list += load(folder, 'sigma_list.p')
    npoints_list += load(folder, 'npoints_list.p')

n_param = len(params_dic)
n_inter = int(len(vc_list) / len(params_dic))

for j in range(n_param):
    sd = []
    mean = []
    vc = []
    for i in range(n_inter):
        index = i*len(params_dic)+j
        vc_name = vc_list[index][0]
        initial = vc_name[0]
        sd.append(sd_list[index][1])
        mean.append(mean_list[index][1])
        vc.append(vc_list[index][1])            
        
        i+=1

#    print('vc_name',vc_name)     
#    print('initial',initial)
#    print('sd',sd)
#    print('mean', mean) 
#    print('vc',vc)

    sigma = np.asarray(sigma_list)
    npoints = np.asarray(npoints_list)
    
    if initial == 'm' or initial == 'M':
        # Narrowing down to points with variance coefficient below p%
        vc = np.asarray(vc)
        
        all_p_index = np.where(abs(vc) < p)   # Indicies of rows with vc < p%.
        all_p_vc = vc[all_p_index]
        all_p_npoints = npoints[all_p_index]
        all_p_sigma = sigma[all_p_index]
    
        # Removing all but the noisiest run for each dataset size.
        stack = np.stack((all_p_npoints,all_p_sigma), axis=1)
        print('stack',stack)
        
        for i in range(len(stack)-1):
            j = i+1 # next line
            np = stack[i][0]
            if np == stack[j][0]:
                print('np',np)
                print('np+1',stack[j][0])
                print('sigma=', str(stack[i][1]), ',sigma+1=', str(stack[j][1]))
                print('i =',i,'j =',j)
                sigma_i = stack[i][1]
                sigma_j = stack[j][1]
                if sigma_i > sigma_j:
                    print('stack[j] = ', str(stack[j]))                    
                    stack = np.delete(stack, j, 0)
                elif sigma_i < sigma_j:
                    print('stack[i] = ', str(stack[i]))
                    stack = np.delete(stack, i, 0)
            i+=1
        
        print('stack',stack)
        
        fig, ax = plt.subplots()
        ax.scatter(all_p_npoints, all_p_sigma, c='m', 
                   label=('within %s percent'%(p)))
        
        # Plotting SD vs dataset size for vc < p.
        for i, txt in enumerate(all_p_vc):
            txt = str(round(txt,2))
            ax.annotate(txt, (all_p_npoints[i], all_p_sigma[i]))
            
        plt.xlabel('Dataset size')
        plt.ylabel('Sigma of noise added to data')
        plt.title('Runs where '+initial+' was found within '+
                  str(p)+' percent'+'\n vc values annotated')  
        plt.legend()      
    
    else:
        # Collecting results foudn with sd < x.
        sd = np.asarray(sd)

        all_x_index = np.where(abs(sd) < x)   # Indicies of rows with sd < x.
        all_x_sd = sd[all_x_index]
        all_x_npoints = npoints[all_x_index]
        all_x_sigma = sigma[all_x_index]
    
         # Removing all but the noisiest run for each dataset size.
         
        fig, ax = plt.subplots()
        ax.scatter(all_x_npoints, all_x_sigma, c='c', 
                   label=('standard deviation < %s'%(x)))
        
        # Plotting SD vs dataset size for vc < p.
        for i, txt in enumerate(all_x_sd):
            txt = str(round(txt,2))
            ax.annotate(txt, (all_x_npoints[i], all_x_sigma[i]))
            
        plt.xlabel('Dataset size')
        plt.ylabel('Sigma of noise added to data')
        plt.title('Runs where sd on '+initial+' found was <'+
                  str(x)+'\n sd values annotated')  
        plt.legend()         
            
        
    j+=1

plt.show()