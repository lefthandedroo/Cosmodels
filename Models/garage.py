#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""
import time
import math
import numpy as np
import scipy as sp

ti = time.time()   
i = 0
while i < 1000000:
    np.log10(10)
    i+=1
tf = time.time()
print(tf-ti)

ti = time.time()   
i = 0
while i < 1000000:
    sp.log10(10)
    i+=1
tf = time.time()
print(tf-ti)

ti = time.time()
i = 0
while i < 1000000:
    math.log10(10)
    i+=1
tf = time.time()
print(tf-ti)



