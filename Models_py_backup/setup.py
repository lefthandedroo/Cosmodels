#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 13:23:20 2018

@author: BallBlueMeercat
"""

from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize('firstderivs_cython.pyx'))