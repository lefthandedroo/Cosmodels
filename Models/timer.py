#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 17:06:45 2018

@author: BallBlueMeercat
"""

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