#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 14:49:28 2019

@author: BallBlueMeercat
"""
import sys
import emcee
from emcee.utils import MPIPool

import time
from example import test_func

try:
    pool = MPIPool()
    if not pool.is_master():
        print('worker')
        pool.wait()
        sys.exit(0)
except:
    pool=None
    pass

# Choose the "true" parameters.
m_true = -0.9594
b_true = 4.294
f_true = 0.534

samples = test_func(m_true, b_true, f_true, pool)

print(samples)

f = open("new_output.txt","w+")
f.write(str(samples))
#f.write("helloo world. testing from login folder")
f.close()

try:
    pool.close()
except:
    pass

print('time is ',time.time())