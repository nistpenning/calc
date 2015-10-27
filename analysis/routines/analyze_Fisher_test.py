# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 14:43:24 2015

@author: rabi
"""

import os, shutil, importlib
import numpy as np
import numpy.random as rand
from numpy import pi, sqrt, sin
import matplotlib.pyplot as plt

import hfGUIdata as hf
import plot_style as ps
import plot_model_fit as pt
importlib.reload(pt)

samples = 1000
num_tip = 25
N = 106.0
std_dev=1/sqrt(N)
z_fluc = N/2.0*sin(std_dev)

#define scale
#bin_def = np.arange(-N/2.0,N/2.0,(N/num_bins))
bin_width_Strobel = 4.0/N  # width for z = 2*S_z/N
bin_def = np.arange(-N/2.0,N/2.0,bin_width_Strobel*N/2.0) 

#create reference histogram
data_ref = rand.normal(loc=0.0, scale=z_fluc, size=samples)
vals_ref, bins, patches = plt.hist(data_ref,bin_def,normed=True)
plt.show()
plt.close()

#create fake data
tipping_angle = np.linspace(0.1,2.1,num=num_tip)
data = np.array([rand.normal(loc= (N/2.0*sin(i/180*pi)), scale=z_fluc, size=samples) for i in tipping_angle])

data_hist = []
n_bin_fills = []
samps = []
for row in data:
    vals = plt.hist(row, bin_def,normed=True)[0]
    
    samp = np.size(row)
    n_bin_fill = np.size(vals[vals!=0])  #number of discrete bins for which P!=0
    
    data_hist.append(vals)
    n_bin_fills.append(n_bin_fill)
    samps.append(samp)
data_hist = np.array(data_hist)
diff_hist = np.array([(sqrt(vals_ref) - sqrt(row)) for i,row in enumerate(data_hist) ])
sq_hist = diff_hist**2
h_dist = np.array([np.sum(row) for row in sq_hist])
plt.close()

l=['Tip angle (rad)', 'dh^2','Extract Fisher Info']
h_err = sqrt(N)*(tipping_angle*pi/180.0)/sqrt(8*np.mean(samps))
fit_res = pt.plot_polyfit(tipping_angle*pi/180.0, h_dist,np.array([0,0,1]),
                        yerr=h_err,
                        hold=np.array([False,True,False]),
                        labels=l)
plt.show()
plt.close(0)

for row in diff_hist:
    plt.plot(bins[:-1],row)

k2 = fit_res[0][1]
n = np.mean(n_bin_fills) + 1
M = np.mean(samps)
NormF = k2/N/( (1/8.0) + (n/32/M) ) #note that for 100 ions and 1000 samples, the added term is neglible

print("F = 8*k2, and then F/N is {:.3g}".format(NormF))