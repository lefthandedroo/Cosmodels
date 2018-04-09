#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 14:49:28 2018

@author: BallBlueMeercat
"""

import pickle
import os
from pylab import figure, plot, xlabel, ylabel, title, show


def rsltplot(filename):
    
    if os.path.exists(filename):
        with open(filename,'rb') as rfp: 
            results = pickle.load(rfp)
    
    ttdd = [item[0] for item in results]
    print(ttdd)
    
    pick = input('Enter ID of entry you want to plot: ')
    
    if not pick:
        gamma = results[-1][1]
        m = results[-1][2]
        de = results[-1][3]
        slnprob = results[-1][4]
    else:
        pick = int(pick)
        gamma = results[pick][1]
        m = results[pick][2]
        de = results[pick][3]
        slnprob = results[pick][4]
       
    
    figure()
    xlabel('parameter value')
    ylabel('step number')
    plot(gamma, slnprob, '.')
    title('slnprob for gamma')
    show()
    
    figure()
    xlabel('parameter value')
    ylabel('step number')
    plot(m, slnprob, '.', color='red')
    title('slnprob for m')
    show()
    
    figure()
    xlabel('parameter value')
    ylabel('step number')
    plot(de, slnprob, '.', color='green')
    title('slnprob for de')
    show()
    
    return gamma, m, de, slnprob

rsltplot('resultlist.p')