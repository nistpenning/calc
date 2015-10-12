# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 20:36:13 2015

@author: jgb
"""

import os, shutil, importlib
import numpy as np
from numpy import pi, sqrt
import matplotlib.pyplot as plt

import hfGUIdata as hf
importlib.reload(hf)
import plot_style as ps
importlib.reload(ps)
import squeeze_func_time as squ

#options
verbose = True
save = True
img_name = "perp_hist_plot_8_11_wPrediction.pdf"
base_path = os.getcwd()
data_path = '/Volumes/688/Public/penning_britton/dailyLabBookFiles/2015/20150811/Load306/squ_hist_after_811'
os.chdir(data_path)

#theory calc info
N=21

# containers for data sets
psis=[]
its=[]
sig_obs = []
sig_ins = []
sig_pns = []
SE = []
Ns = []
names = []

fns = [os.listdir(data_path)[i] for i in [3,12]]
J1ks = (475.0*3.03)*np.ones(np.shape(fns))
Ncals = 1.3999 * np.ones(np.shape(fns))  # #photons per ion per ms
bs = np.arange(-1.1,1.1,0.1)

for i,fn in enumerate(fns):
    folder = os.path.join(data_path,fn)
    os.chdir(folder)
    files = os.listdir(os.getcwd())
    print(folder)

    max_c = hf.get_ionProp_value('detection%det_brightMean')
    min_c = hf.get_ionProp_value('sf%fitParams%sf_fitParam_darkMean')

    file_name, scandata, counts_data, data = hf.get_raw_counts_hist()

    rot_time = np.mean(data['final_t'])
    final_phase = np.mean(data['final_p'])
    tau = np.mean(data['arm_t']) * 2e-3  # ms

    # Calculate derived quantities
    #detune = det_n*1e3/(arm_time)  # kHz
    #phi = (rot_time/pi_time*180.0)  # degrees

    z_data = 2*(((counts_data-min_c)/(max_c - min_c)) - 0.5)
    
    lab = r"$\psi$={0:.2f} deg".format(final_phase/pi*180)
    plt.hist(z_data,bs,label=lab,alpha=0.6,align='left')
    
    os.chdir(base_path)
#plt.legend(fontsize=11)
plt.axis([-1.0,1.0,0,60])
#plt.title(r"Spin Delocalization -- $\tau$: {0:.0f} ms, N: {1:.0f} ".format(tau,N))
plt.xlabel(r"Transverse spin projection 2$S_\psi$/N")
plt.ylabel("Experiments")

print(base_path)

data170name = "/Users/jgb/data170.csv"
data170 = np.genfromtxt(data170name,dtype='float',delimiter=',')
data90name = "/Users/jgb/data90.csv"
data90 = np.genfromtxt(data90name,dtype='float',delimiter=',')

Sx170 = data170.T[0]
trials170 = data170.T[1]
Sx90 = data90.T[0]
trials90 = data90.T[1]

plt.plot(Sx170,trials170,color='k')
plt.plot(Sx90,trials90,color=ps.red)
plt.grid('off')

os.chdir(base_path)
if save is True:
    plt.savefig(img_name,bbox='tight',transparent=True)