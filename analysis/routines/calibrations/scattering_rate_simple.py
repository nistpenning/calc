# -*- coding: utf-8 -*-
"""
Created on Mon Jul 06 09:26:44 2015

@author: jgb
"""
import os
import numpy as np
from numpy import pi, sin, cos
import matplotlib.pyplot as plt

import hfGUIdata
import plot_model_fit as pt
import squeeze_func_time as squ

ats = []
cons = []
file_to_use = [-1]

fns = [hfGUIdata.get_immediate_subdirectories(os.getcwd())[i] for i in file_to_use]
# fit guess values
g = np.array([60.0])

def expdecay(t,G):
    return np.exp(-G*t*1e-3)

def gaussdecay(t,GG):
    return np.exp(-(GG*t*1e-3)**2)

for filename in fns:
    os.chdir(filename)
    fn, data = hfGUIdata.get_gen_csv('phase', skip_header=True)
    arm_time = np.array(data.T[0][1:],dtype='float')
    t0_ind = np.where(arm_time==10.0)
    int_time = arm_time*(2e-3)
    A = np.array(data.T[1][1:],dtype='float')
    B = np.array(data.T[2][1:],dtype='float')
    contrast = ((A - B)/(A[t0_ind]-B[t0_ind]))
    ats.append(arm_time)
    cons.append(contrast)
    l = ["Interaction time [ms]","2|S|/N",fn]
    res, res_err = pt.plot_fit(int_time,contrast,expdecay,fitguess=g,labels=l )
    loss = res[0]
    os.chdir('..')
    plt.show()
    plt.close()