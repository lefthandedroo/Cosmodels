#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:50:29 2017

@author: BallBlueMeercat
"""
import pandas as pd
data = pd.read_csv('./data/lcparam_full_long.txt', sep=" ")#,
#                   names = ['#name','zcmb','zhel','dz','mb','dmb','x1','dx1',
#                            'color', 'dcolor','3rdvar','d3rdvar','cov_m_s',
#                            'cov_m_c','cov_s_c','set','ra','dec','biascor'])