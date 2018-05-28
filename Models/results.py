#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 10:31:30 2018

@author: BallBlueMeercat

"""

import pickle
import os.path

def process(propert):
    
    
    return

def save(save_path, output_name, output): 
    # Saving results to run directory.
    filename = output_name+'.p'
    filename = os.path.join(save_path, filename)
    
    pickle.dump(output, open(filename, 'wb'))

    return


def load(save_path, filename):
    # Opening files with name=filename from given directory=save_path    
    load_path = os.path.join(save_path, filename)    
    content = []
    
    try:
        with open(load_path,'rb') as rfp:
            content = pickle.load(rfp)        
    except:
        print("didnt't open",load_path)
        
    return content

