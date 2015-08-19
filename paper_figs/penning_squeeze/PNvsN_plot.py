# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 11:24:06 2015

@author: jgb
"""

import os,importlib, csv
import numpy as np
from numpy import sqrt
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FormatStrFormatter

import hfGUIdata as hf
import plot_style as ps
importlib.reload(ps)

base = os.getcwd()

path = "/Volumes/688/Public/penning_britton/dailyLabBookFiles/2015/20150813"

os.chdir(path)
data = np.genfromtxt("PNvsN_corr.csv",delimiter=',',names=True)
print(data)

Ns = data['N']
Ns_err = data['N_err']
sig_subs = data['sig_sub_rad']

sfe = data['sig_full_err_rad']**2
N_pred = np.linspace(18,200)

fig, ax = plt.subplots()
plt.grid('off')
plt.errorbar([round(n) for n in Ns], np.array(sig_subs)**2, yerr=sfe, xerr=Ns_err, fmt='o')
plt.plot(N_pred,1/np.array(N_pred),'-')
plt.yscale('log')
plt.xscale('log')
plt.axis([0,220, 0.003,0.06])
plt.ylabel(r'Fractional spin noise $(\Delta S_z/|\vec{S}|)^2$')
plt.xlabel('Ion number N')

#spectroscopic enhancement data
SE = data['SE']
SE_N = data['SE_N']
squ_spin_var = SE/SE_N
plt.plot(SE_N, squ_spin_var,'s')

majorLocator = FixedLocator([20,50,100,200])
majorFormatter = FormatStrFormatter('%d')
majorLocatorY = FixedLocator([0.005,0.01,0.06])
majorFormatterY = FormatStrFormatter('%g')
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_major_locator(majorLocatorY)
ax.yaxis.set_major_formatter(majorFormatterY)

os.chdir(base)

plt.tight_layout()

plt.savefig("PNvsN_fig_alone.png",dpi=300,transparent=True)