#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 12:40:48 2018

@author: BallBlueMeercat


"""

import matplotlib.pyplot as plt
import numpy as np
import os
import time

from results import load

def setsizevspread():
    
    return

def stat(hue, var, var_true, var_name, slnprob, zpicks, 
          mag, sigma, nsteps, nwalkers, save_path):
    
    name_l = var_name.lower()
    initial = name_l[:1]
    name_true = initial + '_true'
    hue_name = hue
    hue = 'xkcd:'+hue
    
    # Marginalised distribution histogram.
    plt.figure()
#    plt.xlabel(r'$\{}$'.format(var_name))
    plt.xlabel(var_name)
    plt.title('Marginalised distribution for '+var_name+' \n nsteps: '+str(nsteps)+
          ', noise: '+str(sigma)+', npoints: '+str(len(zpicks)))
    plt.hist(var, 50)
    stamp = str(int(time.time()))
    filename = str(stamp)+'_'+initial+'_mhist__nsteps_'+str(nsteps) \
    +'_nwalkers_'+str(nwalkers)+'_noise_'+str(sigma) \
    +'_numpoints_'+str(len(zpicks))+'.png'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    plt.show()
    
    # Walker steps.
    plt.figure()
    plt.title('flatlnprobability for '+var_name+' \n nsteps: '
          +str(nsteps)+', noise: '+str(sigma)+', npoints: '+str(len(zpicks)))
    plt.plot(var, slnprob, '.', color=hue)
    stamp = str(int(time.time()))
    filename = str(stamp)+'_'+initial+'_steps__nsteps_'+str(nsteps) \
    +'_nwalkers_'+str(nwalkers)+'_noise_'+str(sigma) \
    +'_numpoints_'+str(len(zpicks))+'.png'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    plt.show()
    
    # Chains.
    plt.figure()
    plt.xlabel('step number')
#    plt.ylabel(r'$\{}$'.format(var_name))
    plt.ylabel(var_name)
    plt.title('flatChains with '+name_true+' in '+hue_name+' \n nsteps: '
          +str(nsteps)+', noise: '+str(sigma)+', npoints: '+str(len(zpicks)))
    plt.plot(var.T, '-', color='k', alpha=0.3)
    plt.axhline(var_true, color=hue)
    stamp = str(int(time.time()))
    filename = str(stamp)+'_'+initial+'_chain__nsteps_'+str(nsteps) \
    +'_nwalkers_'+str(nwalkers)+'_noise_'+str(sigma) \
    +'_numpoints_'+str(len(zpicks))+'.png'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    plt.show()
    
    return

def onepercent():
    
    direclist = []
    for d in os.walk('./results/'):
        direclist.append(d[0])
    direclist.pop(0)
    
    m_vc = []
    g_vc = []
    sigma =[]
    npoints = []
    for directory in direclist:
            m_vc += load(directory, 'm_vc.p')
            g_vc += load(directory, 'g_vc.p')
            sigma += load(directory, 'sigma.p')
            npoints += load(directory, 'npoints.p')

    m_vc = np.asarray(m_vc)
    g_vc = np.asarray(g_vc)
    sigma = np.asarray(sigma)
    npoints = np.asarray(npoints)
      
    m_pi = np.where(m_vc < 1)   # Indicies of rows with m_vc < 1%.
    m_pinpoints = npoints[m_pi]
    m_pisigma = sigma[m_pi]

    g_pi = np.where(g_vc < 1.5)   # Indicies of rows with g_vc < 1%.
    g_pinpoints = npoints[g_pi]
    g_pisigma = sigma[g_pi]
    
    m_sinpoints = []  # Single results, removing doubles of nsteps.
    m_sisigma = []
    
    g_sinpoints = []  # Single results, removing doubles of nsteps.
    g_sisigma = []
    
    for i in range(len(m_pinpoints)):
        if m_pinpoints[i] in m_sinpoints:
           index = np.where(m_sinpoints == m_pinpoints[i])
           index = int(index[0])
           if m_pisigma[i] > m_sisigma[index]:
              m_sisigma[index] = m_pisigma[i]
        else:
            m_sinpoints.append(m_pinpoints[i])
            m_sisigma.append(m_pisigma[i])

    for i in range(len(g_pinpoints)):
        if g_pinpoints[i] in g_sinpoints:
           index = np.where(g_sinpoints == g_pinpoints[i])
           index = int(index[0])
           if g_pisigma[i] > g_sisigma[index]:
              g_sisigma[index] = g_pisigma[i]

        else:
            g_sinpoints.append(g_pinpoints[i])
            g_sisigma.append(g_pisigma[i])
    
#    ind = np.ones((10,), bool)
#    ind[n] = False
#    A1 = A[ind,:]
    plt.figure()
    plt.xlabel('dataset size')
    plt.ylabel('sigma of noise added to data')
    plt.title('noisiest runs where m was found within 1%')       
    plt.scatter(m_sinpoints, m_sisigma, c='m', label='1% sd on m')
    plt.scatter(npoints, sigma, c='c', marker='x', label='all runs')
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.xlabel('dataset size')
    plt.ylabel('sigma of noise added to data')
    plt.title('noisiest runs where gamma was found within 1.5%')       
    plt.scatter(g_sinpoints, g_sisigma, c='g', label='sd on gamma')
    plt.scatter(npoints, sigma, c='c', marker='x', label='all runs')
    plt.legend()
    plt.show()
    
    return m_vc, g_vc, sigma, npoints

#m_vc, g_vc, sigma, npoints = onepercent()

def modelcheck(mag, zpicks, plot_var, firstderivs_key):
    
    t, dlpc, dl, a, ombar_m, gamma, ombar_de, ombar_m0, ombar_de0 = plot_var
    
    print('plots.modelcheck:',firstderivs_key, ', gamma =',str(gamma))
    
#    a = list(reversed(a))
#    zpicks = list(reversed(zpicks))
#    ombar_de = list(reversed(ombar_de))
#    ombar_m = list(reversed(ombar_m))
    
    print()
    print('a goes from ',a[-1], 'to ',a[0])
    print('ombar_de goes from ',ombar_de[-1], 'to ',ombar_de[0])
    print('ombar_m goes from ',ombar_m[-1], 'to ',ombar_m[0])
    print('z goes from ',zpicks[-1], 'to ',zpicks[0])
    
    if min(ombar_m) < 0:
        print('unphysical ombar_m', str(min(ombar_m)))
        print(ombar_m)
        plt.figure()
        plt.title('ombar_m')
        plt.plot(ombar_m)
        plt.show()
        index = np.argmin(ombar_m)
        print('ombar_m close to -ve',ombar_m[index-1], ombar_m[index], ombar_m[index+1])
        
    elif min(ombar_de) < 0:
        print('unphysical ombar_de', str(min(ombar_de)))
        print(type(ombar_de))
        
        pos_signal = ombar_de.copy()
        neg_signal = ombar_de.copy()
        
        pos_signal[pos_signal <= 0] = np.nan
        neg_signal[neg_signal > 0] = np.nan
        
        #plotting
        plt.figure()
        plt.title('ombar_de')
        plt.plot(pos_signal, color='r')
        plt.plot(neg_signal, color='b')
        plt.show()
        
#        plt.figure()
#        plt.title('ombar_de')
#        plt.plot(ombar_de)
#        plt.show()
#        negatives = ombar_de[ombar_de < 0]
#        print(negatives)
#        print(len(negatives))
#        index = np.argmin(ombar_de)
        
#    # Scale factor vs redshift.
#    plt.figure()
#    plt.xlabel('redshift $z$')
#    plt.ylabel('a')
#    plt.grid(True)
#    plt.plot(zpicks, a, 'xkcd:crimson', lw=1)
#    plt.title(r'Scale factor evolution, IC: $\bar \Omega_{m0}$ = %s, $\bar \Omega_{DE0}$ =%s, $\gamma$ = %s'
#          %(ombar_m0, ombar_de0, gamma))
#    plt.show()
    
    Hz = []
    for i in range(len(ombar_m)):
        H = (ombar_m[i] + ombar_de[i])**(1/2)
#        if np.isnan(H):
#            print('plots.modelcheck got NaN value for H')
        Hz.append(H)
    print('Hz goes from ',Hz[-1], 'to ',Hz[0])

    # Scale factor vs redshift and expansion rate H vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel('a')
    plt.grid(True)
    plt.plot(zpicks, a, color='xkcd:crimson', lw=1, label='a = scale factor')
    plt.plot(zpicks, Hz, color='xkcd:blue', lw=1, label='H(z)')
    plt.legend()
    plt.title(r'Scale factor evolution, IC: $\bar \Omega_{m0}$ = %s, $\bar \Omega_{DE0}$ =%s, $\gamma$ = %s'
          %(ombar_m0, ombar_de0, gamma))
    plt.show()

    # ombar_m, ombar_de vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel(r'$\bar \Omega $')
    plt.grid(True)
    plt.plot(zpicks, ombar_m, label=r'$\bar \Omega_{m}$', color='xkcd:coral', linestyle=':')
    plt.plot(zpicks, ombar_de, label=r'$\bar \Omega_{DE}$', color='xkcd:aquamarine')
    plt.legend()
    plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, IC: $\bar \Omega_{m0}$ = %s, $\bar \Omega_{DE0}$ =%s, $\gamma$ = %s'
          %(ombar_m0, ombar_de0, gamma))
    plt.show()
    
    # ombar_m, ombar_de vs redshift log x axis.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel(r'$\bar \Omega $')
    plt.grid(True)
    plt.semilogx(zpicks, ombar_m, label=r'$\bar \Omega_{m}$', color='xkcd:coral', linestyle=':')
    plt.semilogx(zpicks, ombar_de, label=r'$\bar \Omega_{DE}$', color='xkcd:aquamarine')
    plt.legend()
    plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, IC: $\bar \Omega_{m0}$ = %s, $\bar \Omega_{DE0}$ =%s, $\gamma$ = %s'
          %(ombar_m0, ombar_de0, gamma))
    plt.show()

#    # ombar_m vs redshift.
#    plt.figure()
#    plt.xlabel('redshift $z$')
#    plt.ylabel('$\Omega_{m0}$')
#    plt.grid(True)
#    plt.plot(zpicks, ombar_m, 'xkcd:coral', lw=1)
#    plt.title('$\Omega_{m}$ evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
#          %(ombar_m0, ombar_de0, gamma))
#    plt.show()
#
#    # ombar_de vs redshift.
#    plt.figure()
#    plt.xlabel('redshift $z$')
#    plt.ylabel('$\Omega_{DE}$')
#    plt.grid(True)
#    plt.plot(zpicks, ombar_de, 'xkcd:aquamarine', lw=1)
#    plt.title('$\Omega_{DE}$ evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
#          %(ombar_m0, ombar_de0, gamma))
#    plt.show()

    # Luminosity distance vs redshift.
    while False:
        plt.figure()
        plt.xlabel('redshift $z$')
        plt.ylabel('$d_L$*($H_0$/c)')
        plt.grid(True)
        plt.plot(zpicks, dl, 'xkcd:lightgreen', lw=1)
        plt.title('$d_L$ evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
              %(ombar_m0, ombar_de0, gamma))
        plt.show()
        break

#    plt.figure()
#    plt.xlabel('redshift $z$')
#    plt.ylabel('pc')
#    plt.grid(True)
#    plt.plot(zpicks, dlpc, 'xkcd:green', lw=1)
#    plt.title('$d_L$ evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
#          %(ombar_m0, ombar_de0, gamma))
#    plt.show()
#
#    dlgpc = dlpc /10**9
#    plt.figure()
#    plt.xlabel('redshift $z$')
#    plt.ylabel('Gpc')
#    plt.grid(True)
#    plt.plot(zpicks, dlgpc, 'xkcd:darkgreen', lw=1)
#    plt.title('$d_L$ evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
#          %(ombar_m0, ombar_de0, gamma))
#    plt.show()

    while True:
        # Redshift vs time.
        plt.figure()
        plt.xlabel('time')
        plt.ylabel('redshift $z$')
    #        plt.axis([0,-0.8,0,5])
        plt.grid(True)
        plt.plot(t, zpicks, 'm', lw=1)
        plt.title('Redshift evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
              %(ombar_m0, ombar_de0, gamma))
        plt.show()
        
        # Scale factor, H vs time.
        plt.figure()
        plt.xlabel('time')
        plt.ylabel('a')
        plt.grid(True)
        plt.plot(t, a, color='xkcd:crimson', lw=1, label='a = scale factor')
        plt.plot(t, Hz, color='xkcd:blue', lw=1, label='H(t)')
        plt.legend()
        plt.title(r'Scale factor evolution, IC: $\bar \Omega_{m0}$ = %s, $\bar \Omega_{DE0}$ =%s, $\gamma$ = %s'
              %(ombar_m0, ombar_de0, gamma))
        plt.show()
    
#        # Magnitude vs redshift.
#        plt.figure()
#        plt.xlabel('redshift $z$')
#        plt.ylabel('magnitude')
#        plt.title('Magnitude evolution')
#        plt.scatter(zpicks, mag, marker='.', lw='1', c='xkcd:tomato')
#        plt.show()
        break

    return