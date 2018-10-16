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
    
    initial = var_name.lower()[:1]
    name_true = var_name[:1] + '_true'
    hue_name = hue
    hue = 'xkcd:'+hue
    
    # Marginalised distribution histogram.
    plt.figure()
#    plt.xlabel(r'$\{}$'.format(name_l))
    plt.xlabel(var_name)
    plt.title('model: '+firstderivs_key+'\n Marginalised distribution of '
              +var_name+' \n nsteps: '+str(nsteps)+', noise: '+str(sigma)
              +', npoints: '+str(len(zpicks))+' '+firstderivs_key)
    plt.hist(var, 50, facecolor=hue)
    stamp = str(int(time.time()))
    filename = str(stamp)+'_'+initial+'_mhist__nsteps_'+str(nsteps) \
    +'_nwalkers_'+str(nwalkers)+'_noise_'+str(sigma) \
    +'_numpoints_'+str(len(zpicks))+'.png'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    
    # Walker steps.
    plt.figure()
    plt.xlabel(var_name)
    plt.title('model: '+firstderivs_key+'\n lnprobability of '+var_name
              +' \n nsteps: '+str(nsteps)+', noise: '+str(sigma)
              +', npoints: '+str(len(zpicks)))
    plt.plot(var, slnprob, '.', color=hue)
    stamp = str(int(time.time()))
    filename = str(stamp)+'_'+initial+'_steps__nsteps_'+str(nsteps) \
    +'_nwalkers_'+str(nwalkers)+'_noise_'+str(sigma) \
    +'_numpoints_'+str(len(zpicks))+'.png'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    
    # Chains.
    plt.figure()
    plt.xlabel('step number')
#    plt.ylabel(r'$\{}$'.format(name_l))
    plt.ylabel(var_name)
    plt.title('model: '+firstderivs_key+'\n flatchains, '+name_true+
              ' in '+hue_name+' \n nsteps: '+str(nsteps)+', noise: '
              +str(sigma)+', npoints: '+str(len(zpicks)))
    plt.plot(var.T, '-', color='k', alpha=0.3)
    plt.axhline(var_true, color=hue)
    stamp = str(int(time.time()))
    filename = str(stamp)+'_'+initial+'_chain__nsteps_'+str(nsteps) \
    +'_nwalkers_'+str(nwalkers)+'_noise_'+str(sigma) \
    +'_numpoints_'+str(len(zpicks))+'.png'
    filename = os.path.join(save_path, filename)
    plt.savefig(filename)
    
    plt.show(block=False)
    
    return

firstderivs_key = 'exotic'
params_dic = [{'matter':0.3},{'Mcorr':-19.3}, {'gamma':0.0}, {'zeta':0.0}]
def p_percent(p):
    '''
    Takes in:
        p = int/float, percent, precision to which you want parameters;
        firstderivs_key = string, model used;
        params_dic = parameters searched for.
    '''
    
    directory = os.path.join('./results_error_vs_data/'+firstderivs_key)
    
    vc_list = []
    sd_list = []
    mean_list = []
    sigma_list = []
    npoints_list = []
    
    folders = []
    for d in os.walk(directory):
        folders.append(d[0])
    folders.pop(0)
        
    for folder in folders:
        vc_list += load(folder, 'vc_list.p')
        sd_list += load(folder, 'sd_list.p')
        mean_list += load(folder, 'mean_list.p')
        sigma_list += load(folder, 'sigma_list.p')
        npoints_list += load(folder, 'npoints_list.p')
    
    n_param = len(params_dic)
    n_inter = int(len(sd_list) / len(params_dic))
    
    for j in range(n_param):
        sd = []
        mean = []
        vc = []
        for i in range(n_inter):
            index = i*len(params_dic)+j
            
            sd_name = sd_list[index][0]
#            sd_initial = sd_name[0]
            sd.append(sd_list[index][1])
            
#            mean_name = mean_list[index][0]
#            mean_initial = mean_name[0]
            mean.append(mean_list[index][1])
            
#            vc_name = vc_list[index][0]
#            vc_initial = vc_name[0]
            vc.append(vc_list[index][1])            
            
            i+=1
        
        # Narrowing down to points with variance coefficient below p%
        vc = np.asarray(vc)
        sigma = np.asarray(sigma_list)
        npoints = np.asarray(npoints_list)
        
        p1_index = np.where(vc < p)   # Indicies of rows with vc < p%.
        p1_npoints = npoints[p1_index]
        p1_sigma = sigma[p1_index]
        # Removing doubles of nsteps.
        p1_snpoints = []  
        p1_ssigma = []
        
        for i in range(len(p1_npoints)):
            if p1_npoints[i] in p1_snpoints:
               index = np.where(p1_snpoints == p1_npoints[i])
               index = int(index[0])
               if p1_sigma[i] > p1_ssigma[index]:
                  p1_ssigma[index] = p1_sigma[i]
            else:
                p1_snpoints.append(p1_npoints[i])
                p1_ssigma.append(p1_sigma[i])

        plt.figure()
        plt.xlabel('Dataset size')
        plt.ylabel('Sigma of noise added to data')
        plt.title('Noisiest runs where %s was found within %s'%(sd_name, p))       
        plt.scatter(p1_snpoints, p1_ssigma, c='m', label='1% sd')
        plt.scatter(p1_snpoints, sigma, c='c', marker='x', label='all runs')
        plt.legend()
        
#        fig, ax = plt.subplots()
#        ax.scatter(npoints_list, sd, c='r')
#        
#        # Plotting SD vs dataset size.
#        for i, txt in enumerate(sigma_list):
#            txt = 'sd = '+ str(txt)
#            ax.annotate(txt, (npoints_list[i], sd[i]))
#            
#        plt.xlabel('Dataset size')
#        plt.ylabel('s.d. of a marginalised distribution')
#        plt.title(sd_name+' vs dataset size'+
#                  '\n s.d. of noise labeled, model '+firstderivs_key)
##        stamp = str(int(time.time()))
##        filename = str(stamp)+'_sd_of_'+sd_initial+'_.png'
##        filename = os.path.join(save_path, filename)
##        plt.savefig(filename)
#        
#        # Plotting mean vs dataset size.
#        fig, ax = plt.subplots()
#        ax.scatter(npoints_list, mean, c='c')
#        for i, txt in enumerate(sigma_list):
#            txt = 'sd = '+ str(txt)
#            ax.annotate(txt, (npoints_list[i], mean[i]))
#            
#        plt.xlabel('Dataset size')
#        plt.ylabel('Mean of a marginalised distribution')
#        plt.title(mean_name+' vs dataset size'+
#                  '\n s.d. of noise labeled, model '+firstderivs_key)
##        stamp = str(int(time.time()))
##        filename = str(stamp)+'_mean_of_'+mean_initial+'_.png'
##        filename = os.path.join(save_path, filename)
##        plt.savefig(filename)
#        
#        # Plotting variance coefficient vs dataset size.
#        if len(vc) == n_inter:
#            fig, ax = plt.subplots()
#            ax.scatter(npoints_list, vc, c='g')
#            for i, txt in enumerate(sigma_list):
#                txt = 'sd = '+ str(txt)
#                ax.annotate(txt, (npoints_list[i], vc[i]))
#            
#            plt.xlabel('Dataset size')
#            plt.ylabel('s.d. / mean of a marginalised distribution')
#            plt.title(vc_name+' vs dataset size'+
#                      '\n s.d. of noise labeled, model '+firstderivs_key)
##            stamp = str(int(time.time()))
##            filename = str(stamp)+'_vc_of_'+vc_initial+'_.png'
##            filename = os.path.join(save_path, filename)
##            plt.savefig(filename)
            
        j+=1
    
    plt.show()
    
    return

p_percent(2)

def onepercentt():
    
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
      

def modelcheck(mag, zpicks, plot_var, firstderivs_key):
    
    t = plot_var.get('t')
    dl = plot_var.get('dl')
    a = plot_var.get('a')
    ombar_m = plot_var.get('ombar_m')
    ombar_r = plot_var.get('ombar_r')
    gamma = plot_var.get('gamma')
    zeta = plot_var.get('zeta')
    ombar_de = plot_var.get('ombar_de')
    
    t = -t
    
    if min(ombar_m) < 0:
        print('unphysical ombar_m', str(min(ombar_m)))
    elif min(ombar_de) < 0:
        print('unphysical ombar_de', str(min(ombar_de)))    
    
    # Scale factor vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel('a')
    plt.grid(True)
    plt.plot(zpicks, a, 'xkcd:crimson', lw=1)
    plt.title('Scale factor evolution, model = %s, $\gamma$ = %s, $\zeta$ = %s'
              %(firstderivs_key, gamma, zeta))
    
    Hz = []
    for i in range(len(ombar_m)):
        H = (ombar_m[i] + ombar_de[i])**(1/2)
        if np.isnan(H):
            print('plots.modelcheck got NaN value for H')
        Hz.append(H)

    # H vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel('H')
    plt.grid(True)
    plt.plot(zpicks, Hz, color='xkcd:blue', lw=1)
    plt.title(r'H evolution, model = %s, $\gamma$ = %s, $\zeta$ = %s'
              %(firstderivs_key, gamma, zeta))
    
    if ombar_r is None:
        # ombar_m, ombar_de vs redshift.
        plt.figure()
        plt.xlabel('redshift $z$')
        plt.ylabel(r'$\bar \Omega $')
        plt.grid(True)
        plt.plot(zpicks, ombar_m, label=r'$\bar \Omega_{m}$', 
                 color='xkcd:coral', linestyle=':')
        plt.plot(zpicks, ombar_de, label=r'$\bar \Omega_{DE}$', 
                 color='xkcd:aquamarine')
        plt.legend()
        plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, model = %s, $\gamma$ = %s'
              %(firstderivs_key, gamma))
        
        # ombar_m, ombar_de vs redshift log x axis.
        plt.figure()
        plt.xlabel('redshift $z$')
        plt.ylabel(r'$\bar \Omega $')
        plt.grid(True)
        plt.semilogx(zpicks, ombar_m, label=r'$\bar \Omega_{m}$', 
                     color='xkcd:coral', linestyle=':')
        plt.semilogx(zpicks, ombar_de, label=r'$\bar \Omega_{DE}$', 
                     color='xkcd:aquamarine')
        plt.legend()
        plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, model = %s, $\gamma$ = %s'
              %(firstderivs_key, gamma))
        
    else:
        # ombar_m, ombar_de, ombar_r evolution
        plt.figure()
        plt.xlabel('redshift $z$')
        plt.ylabel(r'$\bar \Omega $')
        plt.grid(True)
        plt.plot(zpicks, ombar_m, label=r'$\bar \Omega_{m}$', 
                 color='xkcd:coral', linestyle=':')
        plt.plot(zpicks, ombar_de, label=r'$\bar \Omega_{DE}$', 
                 color='xkcd:aquamarine')
        plt.plot(zpicks, ombar_r, label=r'$\bar \Omega_{r}$', 
                 color='xkcd:black')
        plt.legend()
        plt.title(r'$\bar \Omega_{r}$ evolution, model = %s, $\gamma$ = %s, $\zeta$ = %s'
              %(firstderivs_key, gamma, zeta)) 

    # Luminosity distance vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel('$d_L$*($H_0$/c)')
    plt.grid(True)     
    plt.plot(zpicks, dl, 'xkcd:lightgreen', lw=1)
    plt.title('$d_L$ evolution, model = %s, $\gamma$ = %s, $\zeta$ = %s'
              %(firstderivs_key, gamma, zeta))

    # Redshift vs -time.
    plt.figure()
    plt.xlabel('age')
    plt.ylabel('redshift $z$')
    plt.grid(True)
    plt.plot(t, zpicks, 'm', lw=1)
    plt.title('Redshift evolution, model = %s, $\gamma$ = %s, $\zeta$ = %s'
          %(firstderivs_key, gamma, zeta))

    # Scale factor vs -time.
    plt.figure()
    plt.xlabel('age')
    plt.ylabel('a')
    plt.grid(True)
    plt.plot(t, a, color='xkcd:crimson', lw=1)
    plt.title('Scale factor evolution, model = %s, $\gamma$ = %s, $\zeta$ = %s'
          %(firstderivs_key, gamma, zeta))
    
    # Magnitude vs redshift.
    plt.figure()
    plt.xlabel('redshift $z$')
    plt.ylabel('magnitude')
    plt.title('Magnitude evolution, model = %s, $\gamma$ = %s, $\zeta$ = %s'
          %(firstderivs_key, gamma, zeta))
    plt.scatter(zpicks, mag, marker='.', lw='1', c='xkcd:tomato')

    plt.show()
    return

def multi_modelcheck(zpicks, firstderivs_key, plot_var_list, label):
    
    mag = []
    t = []
    dl = []
    a = []
    ombar_m = []
    ombar_r = []
    ombar_de = []
    ombar_m0 = []
    ombar_r0 = []
    ombar_de0 = []
    gamma = []
    zeta = []
    
    for i in range(len(plot_var_list)):
        for plot_var in plot_var_list:
            mag.append(plot_var.get('mag'))
            t.append(plot_var.get('t'))
            dl.append(plot_var.get('dl'))
            a.append(plot_var.get('a'))
            ombar_m.append(plot_var.get('ombar_m'))
            ombar_r.append(plot_var.get('ombar_r'))
            ombar_de.append(plot_var.get('ombar_de'))
            ombar_m0.append(plot_var.get('ombar_m0'))
            ombar_r0.append(plot_var.get('ombar_r0'))
            ombar_de0.append(plot_var.get('ombar_de0'))
            gamma.append(plot_var.get('gamma'))
            zeta.append(plot_var.get('zeta'))
    
    # Changing time into age
    t[0], t[1], t[2] = -t[0], -t[1], -t[2]
    
    # Scale factor vs redshift.
    fig, ax = plt.subplots()
    plt.xlabel('redshift $z$')
    plt.ylabel('a')
    plt.grid(True)
    ax.plot(zpicks, a[0], 'r-', label=label[0])
    ax.plot(zpicks, a[1], 'g-.', label=label[1])
    ax.plot(zpicks, a[2], 'b:', label=label[2])
    plt.title('Scale factor evolution, model = %s'%(firstderivs_key))
    ax.legend()
    
    # H vs redshift.
    Hz_1 = (ombar_m[0] + ombar_de[0])**(1/2)
    Hz_2 = (ombar_m[1] + ombar_de[1])**(1/2)
    Hz_3 = (ombar_m[2] + ombar_de[2])**(1/2)
    
    fig, ax = plt.subplots()
    plt.xlabel('redshift $z$')
    plt.ylabel('H')
    plt.grid(True)
    ax.plot(zpicks, Hz_1, 'r-', label=label[0])
    ax.plot(zpicks, Hz_2, 'g-.', label=label[1])
    ax.plot(zpicks, Hz_3, 'b:', label=label[2])
    plt.title('H evolution, model = %s'%(firstderivs_key))
    ax.legend()

    if ombar_r[0] is None:
        
        # ombar_m, ombar_de vs redshift.
        fig, ax = plt.subplots()
        plt.xlabel('redshift $z$')
        plt.ylabel(r'$\bar \Omega $')
        plt.grid(True)
        ax.plot(zpicks, ombar_m[0], 'r:', label='m')
        ax.plot(zpicks, ombar_m[1], 'g:', label='m')
        ax.plot(zpicks, ombar_m[2], 'b:', label='m')
        ax.plot(zpicks, ombar_de[0], 'm-.', label='de')
        ax.plot(zpicks, ombar_de[1], 'k-.', label='de')
        ax.plot(zpicks, ombar_de[2], 'c-.', label='de')
        plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, model = %s'%(firstderivs_key))
        ax.legend()
    
        # ombar_m, ombar_de vs redshift log x axis.
        fig, ax = plt.subplots()
        plt.xlabel('redshift $z$')
        plt.ylabel(r'$\bar \Omega $')
        plt.grid(True)
        ax.semilogx(zpicks, ombar_m[0], 'r:', label='m')
        ax.semilogx(zpicks, ombar_m[1], 'g:', label='m')
        ax.semilogx(zpicks, ombar_m[2], 'b:', label='m')
        ax.semilogx(zpicks, ombar_de[0], 'm-.', label='de')
        ax.semilogx(zpicks, ombar_de[1], 'k-.', label='de')
        ax.semilogx(zpicks, ombar_de[2], 'c-.', label='de')
        plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, model = %s'%(firstderivs_key))
        ax.legend()

    else:        
        # ombar_r, ombar_m, ombar_de vs redshift.
        fig, ax = plt.subplots()
        plt.xlabel('redshift $z$')
        plt.ylabel(r'$\bar \Omega $')
        plt.grid(True)
        ax.plot(zpicks, ombar_r[0], 'r:', label='r')
        ax.plot(zpicks, ombar_r[1], 'g:', label='r')
        ax.plot(zpicks, ombar_r[2], 'b:', label='r')
        ax.plot(zpicks, ombar_m[0], 'r-.', label='m')
        ax.plot(zpicks, ombar_m[1], 'g-.', label='m')
        ax.plot(zpicks, ombar_m[2], 'b-.', label='m')
        ax.plot(zpicks, ombar_de[0], 'm-', label='de')
        ax.plot(zpicks, ombar_de[1], 'k-', label='de')
        ax.plot(zpicks, ombar_de[2], 'c-', label='de')
        plt.title(r'$\bar \Omega_{r}$, $\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, model = %s'%(firstderivs_key))
        ax.legend()
    
        # ombar_r, ombar_m, ombar_de vs redshift log x axis.
        fig, ax = plt.subplots()
        plt.xlabel('redshift $z$')
        plt.ylabel(r'$\bar \Omega $')
        plt.grid(True)
        ax.semilogx(zpicks, ombar_r[0], 'r:', label='r')
        ax.semilogx(zpicks, ombar_r[1], 'g:', label='r')
        ax.semilogx(zpicks, ombar_r[2], 'b:', label='r')        
        ax.semilogx(zpicks, ombar_m[0], 'r-.', label='m')
        ax.semilogx(zpicks, ombar_m[1], 'g-.', label='m')
        ax.semilogx(zpicks, ombar_m[2], 'b-.', label='m')
        ax.semilogx(zpicks, ombar_de[0], 'm-', label='de')
        ax.semilogx(zpicks, ombar_de[1], 'k-', label='de')
        ax.semilogx(zpicks, ombar_de[2], 'c-', label='de')
        plt.title(r'$\bar \Omega_{m}$, $\bar \Omega_{DE}$ evolution, model = %s'%(firstderivs_key))
        ax.legend()
    
        # ombar_r vs redshift.
        fig, ax = plt.subplots()
        plt.xlabel('redshift $z$')
        plt.ylabel('$\Omega_{m}$')
        plt.grid(True)
        ax.plot(zpicks, ombar_r[0], 'r-', label=label[0])
        ax.plot(zpicks, ombar_r[1], 'g-.', label=label[1])
        ax.plot(zpicks, ombar_r[2], 'b:', label=label[2])
        plt.title('$\Omega_{m}$ evolution, model = %s'%(firstderivs_key))
        ax.legend()
	
    # Luminosity distance vs redshift.
    fig, ax = plt.subplots()
    plt.xlabel('redshift $z$')
    plt.ylabel('$d_L$*($H_0$/c)')
    plt.grid(True)
    ax.plot(zpicks, dl[0], 'r-', label=label[0])
    ax.plot(zpicks, dl[1], 'g-.', label=label[1])
    ax.plot(zpicks, dl[2], 'b:', label=label[2])
    plt.title('$d_L$ evolution, model = %s'%(firstderivs_key))
    ax.legend()

    # Redshift vs age.
    fig, ax = plt.subplots()
    plt.xlabel('age')
    plt.ylabel('redshift $z$')
    plt.grid(True)
    ax.plot(t[0], zpicks, 'r-', label=label[0])
    ax.plot(t[1], zpicks, 'g-.', label=label[1])
    ax.plot(t[2], zpicks, 'b:', label=label[2])
    plt.title('Redshift evolution, model = %s'%(firstderivs_key))
    ax.legend()
    
    # Scale factor vs -time.
    fig, ax = plt.subplots()
    plt.xlabel('age')
    plt.ylabel('a')
    plt.grid(True)
    ax.plot(t[0], a[0], 'r-', label=label[0])
    ax.plot(t[1], a[1], 'g-.', label=label[1])
    ax.plot(t[2], a[2], 'b:', label=label[2])
    plt.title('Scale factor evolution, model = %s'%(firstderivs_key))
    ax.legend()

    # Magnitude vs redshift.
    fig, ax = plt.subplots()
    plt.xlabel('redshift $z$')
    plt.ylabel('magnitude')
    plt.scatter(zpicks, mag[0], c='r', marker='.', lw=0.1, label=label[0])
    plt.scatter(zpicks, mag[1], c='g', marker='.', lw=0.1,  label=label[1])
    plt.scatter(zpicks, mag[2], c='b', marker='.', lw=0.2, label=label[2])
    plt.title('Magnitude evolution, model = %s'%(firstderivs_key))
    ax.legend()
    
    plt.show()
    return