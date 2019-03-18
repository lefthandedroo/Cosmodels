#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 11:15:20 2019

@author: BallBlueMeercat
"""

import os

rootdir = './results_Bfactor/recent/0.07'

for subdir, dirs, files in os.walk(rootdir):
    for directory in dirs:
        print(directory)
        brief_directory = os.path.join(rootdir, directory)
        brief_path = os.path.join(brief_directory, 'brief.txt')
        f=open(brief_path, "r")
        contents = f.read()
        print(contents)
        print('')