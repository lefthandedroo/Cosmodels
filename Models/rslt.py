#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 10:31:30 2018

@author: BallBlueMeercat
"""
import time
import pickle
import os

#def rslt(mag, zpicks, z, dlpc, dl, gamma, e_dash0m, e_dash0de, t, a, 
# a_dot, t_cut, a_cut, a_dotcut, e_dashm, e_dashde):

def rslt(gamma, m, de, slnprob):
    
    filename = 'resultlist.p' #  testresults
    
    if os.path.exists(filename):
        with open(filename,'rb') as rfp: 
            results = pickle.load(rfp)
     
    runid = len(results)
    runkey = 'ID:' + str(runid) + ' ' + str(time.strftime("%x %X"))

    rundata = [runkey, gamma, m, de, slnprob]
    
    results.append(rundata)
#    results = []
    pickle.dump(results, open(filename, 'wb'))
    
    fromfileresults = pickle.load(open(filename, 'rb'))
#    print('fromfileresults is: ',fromfileresults)


    #    # Saving results for plotting later.
    #    f = open('results.txt','a')
    #    f.write(str(mag),)
    #    f.write(str(zpicks),)
    #    f.write(str(z),)
    #    f.write(str(dlpc),)
    #    f.write(str(dl),)
    #    f.write(str(gamma),)
    #    f.write(str(e_dash0m),)
    #    f.write(str(e_dash0de),)
    #    f.write(str(t),)
    #    f.write(str(a),)
    #    f.write(str(a_dot),)
    #    f.write(str(t_cut),)
    #    f.write(str(a_cut),)
    #    f.write(str(a_dotcut),)
    #    f.write(str(e_dashm),)
    #    f.write(str(a_cut),)
    #    f.write(str(e_dashde) + '\n')
    #    f.close()
    
#    obj = mag, zpicks, z, dlpc, dl, gamma, e_dash0m, e_dash0de, t, a, 
# a_dot, t_cut, a_cut, a_dotcut, e_dashm, e_dashde
    
    return fromfileresults

#fromfileresults = rslt(0,0.3,0.7,9999)