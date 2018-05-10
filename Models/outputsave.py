#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 10:31:30 2018

@author: BallBlueMeercat

Saving cheat sheet:
    
# Saving corner plot to file.
numpoints = len(zpicks)
fig.savefig('zz_Day'+str(time.strftime("%j"))+'_nsteps_'+
            str(nsteps)+'_'+'nwalkers_'+str(nwalkers)+'numpoints_'+
            str(numpoints)+str(time.strftime("%c"))+'.png')

"""
import time
import pickle
import os


def filesave(directory): 

    return

def rslt(m, de, gamma, slnprob):
    
    filename = 'datawalkers.p' #  testresults
    
    if os.path.exists(filename):
        with open(filename,'rb') as rfp: 
            results = pickle.load(rfp)
     
    runid = len(results)
    runkey = 'ID:' + str(runid) + ' ' + str(time.strftime("%x %X"))

    rundata = [runkey, gamma, m, de, slnprob]
    
    results.append(rundata)
#    results = []
    pickle.dump(results, open(filename, 'wb'))
    
    walkerdata = pickle.load(open(filename, 'rb'))

    return walkerdata

#walkerdata = rslt(0,0.3,0.7,9999)
    
def samplersave(sampler):
    
    filename = 'datasampler.p' #  testresults
    
    if os.path.exists(filename):
        with open(filename,'rb') as rfp: 
            results = pickle.load(rfp)
     
    runid = len(results)
    runkey = 'ID:' + str(runid) + ' ' + str(time.strftime("%x %X"))

    rundata = [runkey, sampler]
    
    results.append(rundata)
#    results = []
    pickle.dump(results, open(filename, 'wb'))
    
    samplerdata = pickle.load(open(filename, 'rb'))
#    print(samplerdata)

    return samplerdata

#samplerdata = samplersave(0)
