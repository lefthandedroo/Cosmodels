#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 17:06:45 2018

@author: BallBlueMeercat
"""

def flist(start, stop, step):
    """
    Takes in:
        start, stop, step = integers or floats
    Returns:
        zlist = list start to stop with step as increment
    """
#    print('-flist has been called')
    i = 0
    zlist = [start]
    
    while zlist[i] < stop:
        nextvalue = zlist[i] + step
        zlist.append(nextvalue)
        i += 1
        continue

    return zlist


def timer(string,i,f):
    """
    Takes in:
        i = starting time;
        f = finishing time.
    Returns: Time taken in full minutes and seconds.
    """
    sec = f - i     # Total time to run.
    mins, sec= divmod(sec, 60)
#    mins = round((sec / 60),1)  # minutes
#    sec = round((sec % 60),2)   # seconds
    print(str(string),'time:  ',str(int(mins))+'min',str(int(sec))+'s')
    return


def runcount(sigma, sigma_max, sigma_step,
              npoints_min, npoints_max, npoints_step):
    """
    Returns number of interations that evaluator will run.
    """
    run = 0
    
    while sigma < sigma_max:
        npoints = npoints_min
        while npoints < npoints_max:
            npoints += npoints_step
            run += 1
        sigma += sigma_step

    print('Iterations to come:',run)
    return