#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 16:10:05 2018

@author: BallBlueMeercat
"""
import numpy as np
import random

import flist

def zpicks(zmin, zmax, n):
    """
    Takes in:
        zmin = integer lowest redshift;
        zmax = integer highest redshift;
        n = integer number of redshift picks to be generated.
    Rwturns:
        zpicks = array of randomly selected redshifts between zmin and zmax.
    """
#    print('-zpicks has been called')
    zinterval = (zmax - zmin) / (n*2)
    z_opts = flist.flist(zmin, zmax, zinterval)
    zpicks = random.sample(z_opts, n)
    zpicks = np.asarray(zpicks)
    return zpicks