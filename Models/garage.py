#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""
import numpy.random as rng

def wrap(x, a, b):
    assert b > a
    return (x - a)%(b - a) + a


N = 3

ombar_wm = 0.3    # starting from 0.3
ombar_wr = 0.025    # starting from 0.025

density_left = 1.0 - ombar_wm - ombar_wr       # e_total/e_crit = 1

lower = 0.0
upper = density_left/ N

ombar_w01 = rng.rand()      # random number between 0 and 1
ombar_w01 = wrap(ombar_w01, lower, upper)
ombar_w05 = rng.rand()
ombar_w05 = wrap(ombar_w05, lower, upper)
ombar_w08 = rng.rand()
ombar_w08 = wrap(ombar_w08, lower, upper)
ombar_wde = density_left - ombar_w01 - ombar_w05 - ombar_w08

print('ombar_m = ',ombar_wm)
print('ombar_r = ',ombar_wr)
print('ombar_w(-0.1) = ',ombar_w01)
print('ombar_w(-0.5) = ',ombar_w05)
print('ombar_w(-0.8) = ',ombar_w08)
print('ombar_de = ',ombar_wde)

total = ombar_wm +ombar_wr +ombar_w01 +ombar_w05 +ombar_w08 +ombar_wde
print('total =',total)