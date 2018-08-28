#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

amanullah = pd.read_csv('/Users/BallBlueMeercat/Documents/Study/MPhil/Geraint/Cosmodels/Models/data/Amanullah.txt', sep=" ")

mag = amanullah.mb.values
zpicks = amanullah.zhel.values

data = np.stack((mag,x1,colour,zpicks), axis=0)

data.sort(axis=-1)



mag = data[0]
x1 = data[1]
colour = data[2]
zpicks = data[3]

plt.plot(zpicks, mag)
