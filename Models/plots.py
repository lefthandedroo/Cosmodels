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


def stat(hue, var, var_true, var_name, slnprob, zpicks, 
          mag, sigma, nsteps, nwalkers, save_path, firstderivs_key):
    
    name_l = var_name.lower()
    initial = name_l[:1]
    name_true = initial + '_true'
    hue_name = hue
    hue = 'xkcd:'+hue
    
    # Marginalised distribution histogram.
    plt.figure()
#    plt.xlabel(r'$\{}$'.format(name_l))
    plt.xlabel(name_l)
    plt.title('model: '+firstderivs_key+'\n Marginalised distribution of '
              +name_l+' \n nsteps: '+str(nsteps)+', noise: '+str(sigma)
              +', npoints: '+str(len(zpicks))+' '+firstderivs_key)
    plt.hist(var, 50)
    stamp = str(int(time.time()))
    filename = str(stamp)+'_'+initial+'_mhist__nsteps_'+str(nsteps) \
    +'_nwalkers_'+str(nwalkers)+'_noise_'+str(sigma) \
    +'_numpoints_'+str(len(zpicks))+'.pdf'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    
    # Walker steps.
    plt.figure()
    plt.xlabel(name_l)
    plt.title('model: '+firstderivs_key+'\n lnprobability of '+name_l
              +' \n nsteps: '+str(nsteps)+', noise: '+str(sigma)
              +', npoints: '+str(len(zpicks)))
    plt.plot(var, slnprob, '.', color=hue)
    stamp = str(int(time.time()))
    filename = str(stamp)+'_'+initial+'_steps__nsteps_'+str(nsteps) \
    +'_nwalkers_'+str(nwalkers)+'_noise_'+str(sigma) \
    +'_numpoints_'+str(len(zpicks))+'.pdf'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    
    # Chains.
    plt.figure()
    plt.xlabel('step number')
#    plt.ylabel(r'$\{}$'.format(name_l))
    plt.ylabel(name_l)
    plt.title('model: '+firstderivs_key+'\n flatchains, '+name_true+
              ' in '+hue_name+' \n nsteps: '+str(nsteps)+', noise: '
              +str(sigma)+', npoints: '+str(len(zpicks)))
    plt.plot(var.T, '-', color='k', alpha=0.3)
    plt.axhline(var_true, color=hue)
    stamp = str(int(time.time()))
    filename = str(stamp)+'_'+initial+'_chain__nsteps_'+str(nsteps) \
    +'_nwalkers_'+str(nwalkers)+'_noise_'+str(sigma) \
    +'_numpoints_'+str(len(zpicks))+'.pdf'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    
    plt.show(block=False)
    
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
    
def ivcdmcheck(mag, zpicks, firstderivs_key, plot_var_dict):
    firstderivs_key_3 = firstderivs_key[0]   
    firstderivs_key_1 = firstderivs_key[1]
    firstderivs_key_2 = firstderivs_key[2]


    
    t_1, dlpc_1, dl_1, a_1, ombar_m_1, gamma_1, ombar_de_1, ombar_m0_1, ombar_de0_1 = plot_var_dict['plot_var_1']
    t_2, dlpc_2, dl_2, a_2, ombar_m_2, gamma_2, ombar_de_2, ombar_m0_2, ombar_de0_2 = plot_var_dict['plot_var_2']
    t_3, dlpc_3, dl_3, a_3, ombar_m_3, gamma_3, ombar_de_3, ombar_m0_3, ombar_de0_3 = plot_var_dict['plot_var_3']

    # Changing time into age
    t_1, t_2, t_3 = -t_1, -t_2, -t_3

    mag_1 = plot_var_dict['mag_1']
    mag_2 = plot_var_dict['mag_2']
    mag_3 = plot_var_dict['mag_3']

    
    # Scale factor vs redshift.
    # Create plots with pre-defined labels.
    fig, ax = plt.subplots()
    plt.xlabel('redshift $z$')
    plt.ylabel('a')
    plt.grid(True)
    ax.plot(zpicks, a_1, 'r:', label=firstderivs_key_1)
    ax.plot(zpicks, a_2, 'g-.', label=firstderivs_key_2)
    ax.plot(zpicks, a_3, 'b-', label=firstderivs_key_3)
    plt.title('Scale factor evolution, model = %s, $\gamma$ = %s'
          %(firstderivs_key, gamma_1))
    ax.legend()

    Hz_1 = (ombar_m_1 + ombar_de_1)**(1/2)
    Hz_2 = (ombar_m_2 + ombar_de_2)**(1/2)
    Hz_3 = (ombar_m_3 + ombar_de_3)**(1/2)


    # H vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel('H')
    plt.grid(True)
    plt.plot(zpicks, Hz_1, 'r:',
             zpicks, Hz_2, 'g-.', 
             zpicks, Hz_3, 'b-')
    plt.title(r'H evolution, model = %s'
          %(firstderivs_key))

    # ombar_m, ombar_de vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel(r'$\bar \Omega $')
    plt.grid(True)
    plt.plot(zpicks, ombar_m_1, 'r:', zpicks, ombar_de_1, 'r:', 
             zpicks, ombar_m_2, 'g-.', zpicks, ombar_de_2, 'g-.', 
             zpicks, ombar_m_3, 'b-', zpicks, ombar_de_3, 'b-')
    plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, model = %s, $\gamma$ = %s'
          %(firstderivs_key, gamma_1))
    plt.show()
    
    # ombar_m, ombar_de vs redshift log x axis.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel(r'$\bar \Omega $')
    plt.grid(True)
    plt.semilogx(zpicks, ombar_m_1, 'r:', zpicks, ombar_de_1, 'r:', 
                 zpicks, ombar_m_2, 'g-.', zpicks, ombar_de_2, 'g-.', 
                 zpicks, ombar_m_3, 'b-', zpicks, ombar_de_3, 'b-')
    plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, model = %s, $\gamma$ = %s'
          %(firstderivs_key, gamma_1))

    # ombar_m vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel('$\Omega_{m0}$')
    plt.grid(True)
    plt.plot(zpicks, ombar_m_1, 'r:',
             zpicks, ombar_m_2, 'g-.',
             zpicks, ombar_m_3, 'b-')
    plt.title('$\Omega_{m}$ evolution, $\gamma$ = %s'
          %(gamma_1))
    # ombar_de vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel('$\Omega_{DE}$')
    plt.grid(True)
    plt.plot(zpicks, ombar_de_1, 'r:',
             zpicks, ombar_de_2, 'g-.',
             zpicks, ombar_de_3, 'b-')
    plt.title('$\Omega_{de}$ evolution, $\gamma$ = %s'
          %(gamma_1))
    # Luminosity distance vs redshift.
    while True:
        plt.figure()
        plt.xlabel('redshift $z$')
        plt.ylabel('$d_L$*($H_0$/c)')
        plt.grid(True)
        plt.plot(zpicks, dl_1, 'r:',
                 zpicks, dl_2, 'g-.', 
                 zpicks, dl_3, 'b-')
        plt.title('$d_L$ evolution, model = %s, $\gamma$ = %s'
          %(firstderivs_key, gamma_1))
        break


    while True:
        # Redshift vs -time.
        plt.figure()
        plt.xlabel('age')
        plt.ylabel('redshift $z$')
    #        plt.axis([0,-0.8,0,5])
        plt.grid(True)
        plt.plot(t_1, dl_1, 'r:',
                 t_2, dl_2, 'g-.', 
                 t_3, dl_3, 'b-')
        plt.title('Redshift evolution, model = %s, $\gamma$ = %s'
          %(firstderivs_key, gamma_1))
        
        # Scale factor vs -time.
        plt.figure()
        plt.xlabel('age')
        plt.ylabel('a')
        plt.grid(True)
        plt.plot(t_1, a_1, 'r:',
                 t_2, a_2, 'g-.', 
                 t_3, a_3, 'b-')
        plt.title('Scale factor evolution, model = %s, $\gamma$ = %s'
          %(firstderivs_key, gamma_1))
    
        # Magnitude vs redshift.
        plt.figure()
        plt.xlabel('redshift $z$')
        plt.ylabel('magnitude')
        plt.title('Magnitude evolution, $\gamma$ = %s'%(gamma_1))
        plt.scatter(zpicks, mag_1, c='r', marker=',', linewidths=0.1)
        plt.scatter(zpicks, mag_2, c='g', marker='x', linewidths=0.1) 
        plt.scatter(zpicks, mag_3, c='b', marker='.', linewidths=0.2)
        break
    
    plt.show()
    return

def modelcheck(mag, zpicks, plot_var_1, firstderivs_key):
    
    t, dlpc, dl, a, ombar_m, gamma, ombar_de, ombar_m0, ombar_de0 = plot_var_1
    t = -t
    
    print()
    print('a:',a[-1], '---->',a[0])
    print('ombar_de:',ombar_de[-1], '---->',ombar_de[0])
    print('ombar_m:',ombar_m[-1], '---->',ombar_m[0])
    print('z:',zpicks[-1], '---->',zpicks[0])
    
    if min(ombar_m) < 0:
        print('unphysical ombar_m', str(min(ombar_m)))
        print(ombar_m)
        plt.figure()
        plt.title('ombar_m')
        plt.plot(ombar_m)
        index = np.argmin(ombar_m)
        print('ombar_m close to -ve',ombar_m[index-1], ombar_m[index], ombar_m[index+1])
        
    elif min(ombar_de) < 0:
        print('unphysical ombar_de', str(min(ombar_de)))
        
        pos_signal = ombar_de.copy()
        neg_signal = ombar_de.copy()
        
        pos_signal[pos_signal <= 0] = np.nan
        neg_signal[neg_signal > 0] = np.nan
        
        #plotting
        plt.figure()
        plt.title('ombar_de')
        plt.plot(pos_signal, color='r')
        plt.plot(neg_signal, color='b')
        
        plt.figure()
        plt.title('ombar_de')
        plt.plot(ombar_de)
        negatives = ombar_de[ombar_de < 0]
        print('Negative ombar_de: ')
        print(negatives)
        print(len(negatives))
        index = np.argmin(ombar_de)
        
#    # Scale factor vs redshift.
#    plt.figure()
#    plt.xlabel('redshift $z$')
#    plt.ylabel('a')
#    plt.grid(True)
#    plt.plot(zpicks, a, 'xkcd:crimson', lw=1)
#    plt.title('Scale factor evolution, model = %s, $\gamma$ = %s'
#              %(firstderivs_key, gamma))
##    plt.title(r'Scale factor evolution, IC: $\bar \Omega_{m0}$ = %s, $\bar \Omega_{DE0}$ =%s, $\gamma$ = %s'
#          %(ombar_m0, ombar_de0, gamma))
    
    Hz = []
    for i in range(len(ombar_m)):
        H = (ombar_m[i] + ombar_de[i])**(1/2)
#        if np.isnan(H):
#            print('plots.modelcheck got NaN value for H')
        Hz.append(H)
    print('Hz:',Hz[-1], '---->',Hz[0])

    # Scale factor vs redshift and expansion rate H vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel('a')
    plt.grid(True)
    plt.plot(zpicks, a, color='xkcd:crimson', lw=1, label='a = scale factor')
    plt.plot(zpicks, Hz, color='xkcd:blue', lw=1, label='H(z)')
    plt.legend()
    plt.title(r'Scale factor evolution, model = %s, $\gamma$ = %s'
          %(firstderivs_key, gamma))
#    plt.title(r'Scale factor evolution, %s \n IC: $\bar \Omega_{m0}$ = %s, $\bar \Omega_{DE0}$ =%s, $\gamma$ = %s'
#          %(firstderivs_key, ombar_m0, ombar_de0, gamma))

    # ombar_m, ombar_de vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel(r'$\bar \Omega $')
    plt.grid(True)
    plt.plot(zpicks, ombar_m, label=r'$\bar \Omega_{m}$', color='xkcd:coral', linestyle=':')
    plt.plot(zpicks, ombar_de, label=r'$\bar \Omega_{DE}$', color='xkcd:aquamarine')
    plt.legend()
    plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, model = %s, $\gamma$ = %s'
          %(firstderivs_key, gamma))
#    plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, IC: $\bar \Omega_{m0}$ = %s, $\bar \Omega_{DE0}$ =%s, $\gamma$ = %s'
#          %(ombar_m0, ombar_de0, gamma))
    
    # ombar_m, ombar_de vs redshift log x axis.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel(r'$\bar \Omega $')
    plt.grid(True)
    plt.semilogx(zpicks, ombar_m, label=r'$\bar \Omega_{m}$', color='xkcd:coral', linestyle=':')
    plt.semilogx(zpicks, ombar_de, label=r'$\bar \Omega_{DE}$', color='xkcd:aquamarine')
    plt.legend()
    plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, model = %s, $\gamma$ = %s'
          %(firstderivs_key, gamma))    
#    plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, IC: $\bar \Omega_{m0}$ = %s, $\bar \Omega_{DE0}$ =%s, $\gamma$ = %s'
#          %(ombar_m0, ombar_de0, gamma))

#    # ombar_m vs redshift.
#    plt.figure()
#    plt.xlabel('redshift $z$')
#    plt.ylabel('$\Omega_{m0}$')
#    plt.grid(True)
#    plt.plot(zpicks, ombar_m, 'xkcd:coral', lw=1)
#    plt.title('$\Omega_{m}$ evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
#          %(ombar_m0, ombar_de0, gamma))
#
#    # ombar_de vs redshift.
#    plt.figure()
#    plt.xlabel('redshift $z$')
#    plt.ylabel('$\Omega_{DE}$')
#    plt.grid(True)
#    plt.plot(zpicks, ombar_de, 'xkcd:aquamarine', lw=1)
#    plt.title('$\Omega_{DE}$ evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
#          %(ombar_m0, ombar_de0, gamma))

    # Luminosity distance vs redshift.
    while False:
        plt.figure()
        plt.xlabel('redshift $z$')
        plt.ylabel('$d_L$*($H_0$/c)')
        plt.grid(True)     
        plt.plot(zpicks, dl, 'xkcd:lightgreen', lw=1)
        plt.title('$d_L$ evolution, model = %s, $\gamma$ = %s'
              %(firstderivs_key, gamma))
#        plt.title('$d_L$ evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
#              %(ombar_m0, ombar_de0, gamma))
        break

#    plt.figure()
#    plt.xlabel('redshift $z$')
#    plt.ylabel('pc')
#    plt.grid(True)
#    plt.plot(zpicks, dlpc, 'xkcd:green', lw=1)
#    plt.title('$d_L$ evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
#          %(ombar_m0, ombar_de0, gamma))
#
#    dlgpc = dlpc /10**9
#    plt.figure()
#    plt.xlabel('redshift $z$')
#    plt.ylabel('Gpc')
#    plt.grid(True)
#    plt.plot(zpicks, dlgpc, 'xkcd:darkgreen', lw=1)
#    plt.title('$d_L$ evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
#          %(ombar_m0, ombar_de0, gamma))

    while True:
        # Redshift vs -time.
        plt.figure()
        plt.xlabel('age')
        plt.ylabel('redshift $z$')
    #        plt.axis([0,-0.8,0,5])
        plt.grid(True)
        plt.plot(t, zpicks, 'm', lw=1)
        plt.title('Redshift evolution, model = %s, $\gamma$ = %s'
              %(firstderivs_key, gamma))
#        plt.title('Redshift evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
#              %(ombar_m0, ombar_de0, gamma))
        
        # Scale factor vs -time.
        plt.figure()
        plt.xlabel('age')
        plt.ylabel('a')
        plt.grid(True)
        plt.plot(t, a, color='xkcd:crimson', lw=1, label='a = scale factor')
        plt.legend()
        plt.title('Scale factor evolution, model = %s, $\gamma$ = %s'
              %(firstderivs_key, gamma))
#        plt.title(r'Scale factor evolution, IC: $\bar \Omega_{m0}$ = %s, $\bar \Omega_{DE0}$ =%s, $\gamma$ = %s'
#              %(ombar_m0, ombar_de0, gamma))
    
        # Magnitude vs redshift.
        plt.figure()
        plt.xlabel('redshift $z$')
        plt.ylabel('magnitude')
        plt.title('Magnitude evolution')
        plt.scatter(zpicks, mag, marker='.', lw='1', c='xkcd:tomato')
        break
    plt.show()
    return

def gammacheck(mag, zpicks, firstderivs_key, plot_var_dict):
    
    t_1, dlpc_1, dl_1, a_1, ombar_m_1, gamma_1, ombar_de_1, ombar_m0_1, ombar_de0_1 = plot_var_dict['plot_var_1']
    t_2, dlpc_2, dl_2, a_2, ombar_m_2, gamma_2, ombar_de_2, ombar_m0_2, ombar_de0_2 = plot_var_dict['plot_var_2']
    t_3, dlpc_3, dl_3, a_3, ombar_m_3, gamma_3, ombar_de_3, ombar_m0_3, ombar_de0_3 = plot_var_dict['plot_var_3']

    # Changing time into age
    t_1, t_2, t_3 = -t_1, -t_2, -t_3

    mag_1 = plot_var_dict['mag_1']
    mag_2 = plot_var_dict['mag_2']
    mag_3 = plot_var_dict['mag_3']

    
    # Scale factor vs redshift.
    # Create plots with pre-defined labels.
    fig, ax = plt.subplots()
    plt.xlabel('redshift $z$')
    plt.ylabel('a')
    plt.grid(True)
    ax.plot(zpicks, a_1, 'r:', label= str(gamma_1))
    ax.plot(zpicks, a_2, 'g-.', label=str(gamma_2))
    ax.plot(zpicks, a_3, 'b-', label=str(gamma_3))
    plt.title('Scale factor evolution, model = %s' %(firstderivs_key))
    ax.legend()

    Hz_1 = (ombar_m_1 + ombar_de_1)**(1/2)
    Hz_2 = (ombar_m_2 + ombar_de_2)**(1/2)
    Hz_3 = (ombar_m_3 + ombar_de_3)**(1/2)


    # H vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel('H')
    plt.grid(True)
    plt.plot(zpicks, Hz_1, 'r:',
             zpicks, Hz_2, 'g-.', 
             zpicks, Hz_3, 'b-')
    plt.title(r'H evolution, model = %s'
          %(firstderivs_key))

    # ombar_m, ombar_de vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel(r'$\bar \Omega $')
    plt.grid(True)
    plt.plot(zpicks, ombar_m_1, 'r:', zpicks, ombar_de_1, 'r:', 
             zpicks, ombar_m_2, 'g-.', zpicks, ombar_de_2, 'g-.', 
             zpicks, ombar_m_3, 'b-', zpicks, ombar_de_3, 'b-')
    plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, model = %s'
          %(firstderivs_key))
    plt.show()
    
    # ombar_m, ombar_de vs redshift log x axis.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel(r'$\bar \Omega $')
    plt.grid(True)
    plt.semilogx(zpicks, ombar_m_1, 'r:', zpicks, ombar_de_1, 'r:', 
                 zpicks, ombar_m_2, 'g-.', zpicks, ombar_de_2, 'g-.', 
                 zpicks, ombar_m_3, 'b-', zpicks, ombar_de_3, 'b-')
    plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, model = %s'
          %(firstderivs_key))    

    # ombar_m vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel('$\Omega_{m0}$')
    plt.grid(True)
    plt.plot(zpicks, ombar_m_1, 'r:',
             zpicks, ombar_m_2, 'g-.',
             zpicks, ombar_m_3, 'b-')
    plt.title('$\Omega_{m}$ evolution')

    # ombar_de vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel('$\Omega_{DE}$')
    plt.grid(True)
    plt.plot(zpicks, ombar_de_1, 'r:',
             zpicks, ombar_de_2, 'g-.',
             zpicks, ombar_de_3, 'b-')
    plt.title('$\Omega_{de}$ evolution')

    # Luminosity distance vs redshift.
    while True:
        plt.figure()
        plt.xlabel('redshift $z$')
        plt.ylabel('$d_L$*($H_0$/c)')
        plt.grid(True)
        plt.plot(zpicks, dl_1, 'r:',
                 zpicks, dl_2, 'g-.', 
                 zpicks, dl_3, 'b-')
        plt.title('$d_L$ evolution, model = %s'
              %(firstderivs_key))
        break


    while True:
        # Redshift vs -time.
        plt.figure()
        plt.xlabel('age')
        plt.ylabel('redshift $z$')
    #        plt.axis([0,-0.8,0,5])
        plt.grid(True)
        plt.plot(t_1, dl_1, 'r:',
                 t_2, dl_2, 'g-.', 
                 t_3, dl_3, 'b-')
        plt.title('Redshift evolution, model = %s'
              %(firstderivs_key))
        
        # Scale factor vs -time.
        plt.figure()
        plt.xlabel('age')
        plt.ylabel('a')
        plt.grid(True)
        plt.plot(t_1, a_1, 'r:',
                 t_2, a_2, 'g-.', 
                 t_3, a_3, 'b-')
        plt.title('Scale factor evolution, model = %s'
              %(firstderivs_key))
    
        # Magnitude vs redshift.
        plt.figure()
        plt.xlabel('redshift $z$')
        plt.ylabel('magnitude')
        plt.title('Magnitude evolution')
        plt.scatter(zpicks, mag_1, c='r', marker=',', linewidths=0.1)
        plt.scatter(zpicks, mag_2, c='g', marker='x', linewidths=0.1) 
        plt.scatter(zpicks, mag_3, c='b', marker='.', linewidths=0.2)
        break
    
    plt.show()
    return