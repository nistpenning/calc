# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 09:48:06 2015

@author: jgb
"""

import os, shutil, importlib
import numpy as np
from numpy import pi, sqrt
import matplotlib.pyplot as plt
import scipy.stats.mstats as mstats

import hfGUIdata as hf
importlib.reload(hf)
import plot_style as ps
importlib.reload(ps)
import plot_model_fit as pf
importlib.reload(pf)

#options
save = False
base_path = os.getcwd()
save_file_name = "9_30_counts_data"
data_path = '/Users/jgb/Data/20151002/squeeze_raw/165546'
data_path = '/Users/jgb/Data/20150930/Load331/squeeze_raw/182856'
data_path_name = data_path[-6:]
save_file_name  = save_file_name + data_path_name
os.chdir(data_path)

walsh_num = 4  # use 4 for a walsh seq, 2 for regular (mulitplies the arm time)

# containers for data sets
psis=[]
its=[]
z_datas=[]
Ns = []
names = []

fns = os.listdir(data_path)

for i,fn in enumerate(fns):
    folder = os.path.join(data_path,fn)
    os.chdir(folder)
    files = os.listdir(os.getcwd())

    max_c = hf.get_ionProp_value('detection%det_brightMean')
    min_c = hf.get_ionProp_value('detection%det_darkMean')
    pi_time = hf.get_ionProp_value('sf%fitParams%sf_fitParam_tpi')

    file_name, scandata, counts_data, data = hf.get_raw_counts_hist()

    rot_time = np.mean(data['final_t'])
    final_phase = np.mean(data['final_p'])
    psi_deg = final_phase*180/pi  # degrees
    tau = np.mean(data['arm_t']) * walsh_num * 1e-3  # ms  

    # Calculate derived quantities
    #detune = det_n*1e3/(arm_time)  # kHz

    z_data = 2*(((counts_data-min_c)/(max_c - min_c)) - 0.5)
    
    psis.append(psi_deg)
    z_datas.append(z_data)

os.chdir(data_path)
os.chdir('..')

name_mod = "_{:.0g}ms".format(tau)

if save is True:
    psis_names = ','.join(["{:.2f}".format(ele) for ele in psis])  # this formats the list of floats into string    
    pf.save_data_txt(save_file_name+name_mod+".csv", z_datas, col_names=psis_names)
        