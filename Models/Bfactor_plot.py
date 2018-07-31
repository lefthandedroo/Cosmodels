#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 12:10:06 2018

@author: usyd
"""
import numpy as np
import matplotlib.pyplot as plt

filename = 'sample.txt'
array = np.loadtxt(filename)

plt.figure()
plt.title('matter')
plt.hist(array[:, 0])

plt.figure()
plt.title('gamma')
plt.hist(array[:, 1])

plt.show()