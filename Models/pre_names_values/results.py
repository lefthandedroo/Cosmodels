#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 10:31:30 2018

@author: BallBlueMeercat

"""

import pickle
import os.path
import time


def save(save_path, output_name, output): 
    # Saving results to run directory.
    filename = output_name+'.p'
    filename = os.path.join(save_path, filename)
    
    pickle.dump(output, open(filename, 'wb'))

    return


def load(save_path, filename):
    
    if filename[-2:] != '.p':
        print('.p added to',filename,'in',save_path,'- assuming binary file...')
        filename += '.p'
    
    # Opening files with name=filename from given directory=save_path    
    load_path = os.path.join(save_path, filename)    
    content = []
    
    try:
        with open(load_path,'rb') as rfp:
            content = pickle.load(rfp)        
    except:
        print("didnt't open",load_path)
        
    return content


def relocate(filename, speed, firstderivs_key):
    
    filename = filename
    
    # Create directory to move files to.
    directory = './results_Bfactor/'+str(int(time.time()))+'_'+str(speed)+'_model_'+firstderivs_key
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    # Add filename to the new directory    
    new_filename = os.path.join(directory, filename)
    
    os.rename(filename, new_filename)