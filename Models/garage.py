#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""
N = 2
sd_list = [['m_sd', 0.09237825912653962], 
           ['M_sd', 0.0308032706586774], 
           ['a_sd', 2.946726410979336], 
           ['b_sd', 2.7984382158521295], 
           ['g_sd', 0.4967545852017745], 
           ['m_sd', 0.13301192817647778], 
           ['M_sd', 0.03821243057710596], 
           ['a_sd', 3.1558132604018363], 
           ['b_sd', 2.758724508013216], 
           ['g_sd', 0.8494551025469079]]

params_dic = [{'matter':0.3},{'Mcorr':-19.3},{'alpha':0},
                   {'beta':0},{'gamma':0.0}]

print('len(sd_list)',len(sd_list))

#for i in range(N):
#    index = i*len(params_dic)
#    print('index',index)
#    print('sd value',sd_list[index][1])
#    i+=1
    
for j in range(len(params_dic)):  
    for i in range(N):
        index = i*len(params_dic)+j
        print('index',index)
        print('j',j)
        print('sd name',sd_list[index][0])
        print('sd value',sd_list[index][1])
        i+=1
    j+=1