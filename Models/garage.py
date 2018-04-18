#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""

    # Create z samples for the ODE solver.
#    import random
#    numpoints = 100
#    zpicks = random.sample(range(0, numpoints*8), numpoints-1)
#    zpicks.sort()
#    denom = numpoints*10
#    zpicks = [time / denom for time in zpicks]
#    zpicks = [0] + zpicks
#    print('len zpicks after creation', len(zpicks))
#    print('zpicks:')
#    print(zpicks)

def priortest(m, de): 
    flatness = m + de
    if flatness == 1:
        print('flat')
    else:
        print('not flat')
    return

priortest(0.3,-0.7)