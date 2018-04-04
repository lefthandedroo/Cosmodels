#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 10:31:30 2018

@author: BallBlueMeercat
"""
import csv
#import plots

def rslt(mag, zpicks, z, dlpc, dl, gamma, e_dash0m, e_dash0de, t, a, a_dot, t_cut, a_cut, a_dotcut, e_dashm, e_dashde):
    
    with open('results.csv', 'wb') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=' ',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        spamwriter.writerow([mag, zpicks, z, dlpc, dl, gamma, e_dash0m, e_dash0de, t, a, a_dot, t_cut, a_cut, a_dotcut, e_dashm, e_dashde])


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
    
    # Plotting results.
    # plots.plots(mag, zpicks, z, dlpc, dl, gamma, e_dash0m, e_dash0de, t, a, a_dot, t_cut, a_cut, a_dotcut, e_dashm, e_dashde)
    return