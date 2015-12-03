# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 14:42:25 2015

@author: jgb
"""

import os
import numpy as np
from numpy import pi, sqrt, sin
import matplotlib.pyplot as plt

import hfGUIdata as hf
import plot_style as ps
import plot_model_fit as pt

#options
colors = ['k', ps.red, ps.blue, ps.orange, ps.pink]

files_to_use = [-1] 
Ncal = 1.70

base_path = os.getcwd()
fns = [os.listdir(base_path)[i] for i in files_to_use]
fns = [i for i in fns if not i.startswith('.DS_Store')]

#%%
#_____________________________________________________________________
# data processing here
fn = os.listdir(base_path)[files_to_use]
folder = os.path.join(base_path,fn)
os.chdir(folder)
files = os.listdir(os.getcwd())

#Load the data 
file_name, scandata, avg_pmt_counts, pmterr, trials, data = hf.get_raw_counts()
bm = hf.get_ionProp_value('det_brightMean')
dm = hf.get_ionProp_value('det_darkMean')
tpi = hf.get_ionProp_value('sf_fitParam_tpi')
det_t = hf.get_ionProp_value('det_t')
k = bm-dm  # phtns per N atoms
N = k/(det_t*1e-3)/Ncal

os.chdir(base_path)


# analysis here
def n_pi_std_fit(n_pi, y0, A):
    sig_tot_2 = y0**2 + (A*n_pi)**2
    return sqrt(sig_tot_2)

guess=np.array([90.0,1.0])
hold=np.array([False,False])
pl = ["# pi pulses","Std Dev. photon count",'uWave amplitude noise']
pout,perr = pt.plot_fit(scandata,pmterr,n_pi_std_fit,guess,
            hold=hold,
            labels=pl)

plt.show()
