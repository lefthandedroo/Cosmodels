#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 16:17:52 2018

@author: BallBlueMeercat
"""
import numpy as np
import lnprior

# prior:  if -0.1 < gamma < 0.1 and 0.299 < m < 0.301 and 0.699 < de < 0.701:

a = [0.05, 0.01, 0, -10]
b = [0.3, 0.3009, -10, 0.2999]
c = [0.7, -10, 0.6999, 0.7009]
prob = [0.0, -np.inf, -np.inf, -np.inf]

def thetainf(gamma, m, de, slnprob):
    theta = np.column_stack((gamma, m, de))
    
    badtheta_finite_slnprob = 0
    goodtheta_inf_slnprob = 0
    
    badtheta_inf_slnprob = 0
    goodtheta_finite_slnprob = 0
    
    for i in range(len(slnprob)):
        lp = lnprior.lnprior(theta[i])
        
        if not np.isfinite(lp):             # if theta is outside prior
            if np.isfinite(slnprob[i]):     # and has finite slnprob
                print('theta = %s, slnprob = %s'%(theta,slnprob))
                badtheta_finite_slnprob +=1
                
            if not np.isfinite(slnprob[i]): # and has inf slnprob as it should
                badtheta_inf_slnprob +=1
                
        if np.isfinite(lp):                 # if theta is within prior
            if not np.isfinite(slnprob[i]): # but has inf slnprob
                print('theta = %s, slnprob = %s'%(theta,slnprob))
                goodtheta_inf_slnprob +=1
            
            if np.isfinite(slnprob[i]):     # & has finite slnprob as it should
                goodtheta_finite_slnprob +=1

#    wrong = badtheta_finite_slnprob + goodtheta_inf_slnprob
#    right = badtheta_inf_slnprob + goodtheta_finite_slnprob
#    print('Bad theta with finite slnprob found: ',badtheta_finite_slnprob)
#    print('Good theta with inf slnprob found: ',goodtheta_inf_slnprob)
#    print('Bad theta with inf slnprob found: ',badtheta_inf_slnprob)
    print('Good theta with finite slnprob found: ',goodtheta_finite_slnprob)
#    print('wrong = %s, right = %s'%(wrong, right))
    return


import pylab as pl

def sanity(t_rslt, z_rzlt, zpicks, gamma, e_dash0m, e_dash0de):
    
    t, mag, dlpc, dl, a, e_dashm, e_dashde = t_rslt
    zt, zmag, zdlpc, zdl, za, ze_dashm, ze_dashde = z_rzlt

    
    print('len t %s, len mag %s, dlpc %s, dl %s, a %s, e_dashm %s, e_dashde %s'%
          (len(t), len(mag), len(dlpc), len(dl), len(a), len(e_dashm), len(e_dashde)))
    print('lenzt %s, len zmag %s, zdlpc %s, zdl %s, za %s, ze_dashm %s, ze_dashde %s'%
          (len(zt), len(zmag), len(zdlpc), len(zdl), len(za), len(ze_dashm), len(ze_dashde)))

    pl.figure()
    pl.title('a wrt time (blue) and a wrt z (red) vs time')
    pl.xlabel('time')
    pl.ylabel('magnitude')
    pl.plot(t, a, 'b')
    pl.plot(zt, za, 'r')
    pl.show()
    

    pl.figure()
    pl.title('e_dashm wrt time (blue) and e_dashm wrt z (red) vs time')
    pl.xlabel('time')
    pl.ylabel('scale factor')
    pl.plot(t, e_dashm, 'b')
    pl.plot(zt, ze_dashm, 'r')
    pl.show()

    pl.figure()
    pl.title('e_dashde wrt time (blue) and e_dashde wrt z (red) vs time')
    pl.xlabel('time')
    pl.ylabel('scale factor')
    pl.plot(t, e_dashde, 'b')
    pl.plot(zt, ze_dashde, 'r')
    pl.show() 

    pl.figure()
    pl.title('dlpc wrt time (blue) and dlpc wrt z (red) vs time')
    pl.xlabel('time')
    pl.ylabel('scale factor')
    pl.plot(t, dlpc, 'b')
    pl.plot(zt, zdlpc, 'r')
    pl.show()    
    
    
    print(type(mag))
    if not sorted(mag) == mag:
        mag = mag.sort()
        print('mag sorted to accending')
        
    print('mag')
    print(mag)
    
    pl.figure()
    pl.title('mag wrt time (blue) and mag wrt z (red) vs time')
    pl.xlabel('time')
    pl.ylabel('magnitude')
    pl.plot(t, mag, 'b')
    pl.plot(zt, zmag, 'r')
    pl.show()
        

    # Difference in scale factors a:
#    plist = []
#    for i in range(len(zpicks)):
#        
#        difference = a[i] - za[i]
#        plist.append((zpicks[i],difference))
#        
#    plist = sorted(plist,key=lambda x: x[0])
#    
#    parray = np.asarray(plist)
#    zpicks = parray[:,0]
#    diff = parray[:,1]
#    pl.plot(zpicks, diff)
#    pl.xlabel('z')
#    pl.ylabel('a - za')
#    pl.title('(a w.r.t. time - a w.r.t. z) vs redshift')
#    pl.show()
#    
#        # Difference in magnitudes m:
#    plist = []
#    for i in range(len(zpicks)):
#        
#        difference = mag[i] - zmag[i]
#        plist.append((zpicks[i],difference))
#        
#    plist = sorted(plist,key=lambda x: x[0])
#    
#    parray = np.asarray(plist)
#    zpicks = parray[:,0]
#    diff = parray[:,1]
#    pl.plot(zpicks, diff)
#    pl.xlabel('z')
#    pl.ylabel('mag - zmag')
#    pl.title('(mag w.r.t. time - mag w.r.t. z) vs redshift')
#    pl.show()
    return