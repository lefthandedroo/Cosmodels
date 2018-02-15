#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:20:36 2018

@author: BallBlueMeercat
"""

def flist(start, stop, step):
    """
    Takes in:
        start, stop, step - integers or floats
    Returns:
        zlist - a list start to stop with step as increment
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