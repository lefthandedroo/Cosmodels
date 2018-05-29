#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 12:40:48 2018

@author: BallBlueMeercat


"""

from pylab import figure, plot, xlabel, ylabel, title, show, grid, scatter
from pylab import hist, savefig, axhline
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
    figure()
#    xlabel(r'$\{}$'.format(var_name))
    xlabel(var_name)
    title('Marginalised distribution for '+var_name+' \n nsteps: '+str(nsteps)+
          ', noise: '+str(sigma)+', npoints: '+str(len(zpicks)))
    hist(var, 50)
    stamp = str(int(time.time()))
    filename = str(stamp)+'_'+initial+'_mhist__nsteps_'+str(nsteps) \
    +'_nwalkers_'+str(nwalkers)+'_noise_'+str(sigma) \
    +'_numpoints_'+str(len(zpicks))+'.png'
    filename = os.path.join(save_path, filename)
    savefig(filename)
    show()
    
    # Walker steps.
    figure()
    title('slnprob for '+var_name+' \n nsteps: '+str(nsteps)+', noise: '
          +str(sigma)+', npoints: '+str(len(zpicks)))
    plot(var, slnprob, '.', color=hue)
    stamp = str(int(time.time()))
    filename = str(stamp)+'_'+initial+'_steps__nsteps_'+str(nsteps) \
    +'_nwalkers_'+str(nwalkers)+'_noise_'+str(sigma) \
    +'_numpoints_'+str(len(zpicks))+'.png'
    filename = os.path.join(save_path, filename)
    savefig(filename)
    show()
    
    # Chains.
    figure()
    xlabel('step number')
#    ylabel(r'$\{}$'.format(var_name))
    ylabel(var_name)
    title('flatChains with '+name_true+' in '+hue_name+' \n nsteps: '
          +str(nsteps)+', noise: '+str(sigma)+', npoints: '+str(len(zpicks)))
    plot(var.T, '-', color='k', alpha=0.3)
    axhline(var_true, color=hue)
    stamp = str(int(time.time()))
    filename = str(stamp)+'_'+initial+'_chain__nsteps_'+str(nsteps) \
    +'_nwalkers_'+str(nwalkers)+'_noise_'+str(sigma) \
    +'_numpoints_'+str(len(zpicks))+'.png'
    filename = os.path.join(save_path, filename)
    savefig(filename)
    show()
    
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
    figure()
    xlabel('dataset size')
    ylabel('sigma of noise added to data')
    title('noisiest runs where m was found within 1%')       
    plt.scatter(m_sinpoints, m_sisigma, c='m', label='1% sd on m')
    plt.scatter(npoints, sigma, c='c', marker='x', label='all runs')
    plt.legend()
    plt.show()
    
    figure()
    xlabel('dataset size')
    ylabel('sigma of noise added to data')
    title('noisiest runs where gamma was found within 1.5%')       
    plt.scatter(g_sinpoints, g_sisigma, c='g', label='sd on gamma')
    plt.scatter(npoints, sigma, c='c', marker='x', label='all runs')
    plt.legend()
    plt.show()
    
    return m_vc, g_vc, sigma, npoints

m_vc, g_vc, sigma, npoints = onepercent()

def modelcheck(mag, zpicks, plot_var):
    
    t, dlpc, dl, a, ombar_m, gamma, ombar_de, ombar_m0, ombar_de0 = plot_var
    
    # Scale factor vs redshift.
    figure()
    xlabel('redshift $z$')
    ylabel('a')
    grid(True)
    plot(zpicks, a, 'xkcd:crimson', lw=1)
    title(r'Scale factor evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
          %(ombar_m0, ombar_de0, gamma))
    show()

    # ombar_m, ombar_de vs redshift.
    figure()
    xlabel('redshift $z$')
    ylabel('$\Omega_{m0}$')
    grid(True)
    plot(zpicks, ombar_m, 'xkcd:coral', lw=1)
    plot(zpicks, ombar_de, 'xkcd:aquamarine', lw=1)
    title('$\Omega_{m}$, $\Omega_{DE0}$ evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
          %(ombar_m0, ombar_de0, gamma))
    show()

#    # ombar_m vs redshift.
#    figure()
#    xlabel('redshift $z$')
#    ylabel('$\Omega_{m0}$')
#    grid(True)
#    plot(zpicks, ombar_m, 'xkcd:coral', lw=1)
#    title('$\Omega_{m}$ evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
#          %(ombar_m0, ombar_de0, gamma))
#    show()
#
#    # ombar_de vs redshift.
#    figure()
#    xlabel('redshift $z$')
#    ylabel('$\Omega_{DE}$')
#    grid(True)
#    plot(zpicks, ombar_de, 'xkcd:aquamarine', lw=1)
#    title('$\Omega_{DE}$ evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
#          %(ombar_m0, ombar_de0, gamma))
#    show()

    # Luminosity distance vs redshift.
    while True:
        figure()
        xlabel('redshift $z$')
        ylabel('$d_L$*($H_0$/c)')
        grid(True)
        plot(zpicks, dl, 'xkcd:lightgreen', lw=1)
        title('$d_L$ evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
              %(ombar_m0, ombar_de0, gamma))
        show()
        break

#    figure()
#    xlabel('redshift $z$')
#    ylabel('pc')
#    grid(True)
#    plot(zpicks, dlpc, 'xkcd:green', lw=1)
#    title('$d_L$ evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
#          %(ombar_m0, ombar_de0, gamma))
#    show()
#
#    dlgpc = dlpc /10**9
#    figure()
#    xlabel('redshift $z$')
#    ylabel('Gpc')
#    grid(True)
#    plot(zpicks, dlgpc, 'xkcd:darkgreen', lw=1)
#    title('$d_L$ evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
#          %(ombar_m0, ombar_de0, gamma))
#    show()

    # Redshift vs time.
    figure()
    xlabel('time')
    ylabel('redshift $z$')
#        axis([0,-0.8,0,5])
    grid(True)
    plot(t, zpicks, 'm', lw=1)
    title('Redshift evolution, IC: $\Omega_{m0}$ = %s, $\Omega_{DE0}$ =%s, $\gamma$ = %s'
          %(ombar_m0, ombar_de0, gamma))
    show()

    # Magnitude vs redshift.
    figure()
    xlabel('redshift $z$')
    ylabel('magnitude')
    title('Magnitude evolution')
    scatter(zpicks, mag, marker='.', lw='1', c='xkcd:tomato')
    show()

    return