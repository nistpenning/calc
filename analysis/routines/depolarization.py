# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 10:37:46 2015

@author: jgb
"""
import os
import numpy as np
from numpy import pi, sin, cos, sqrt
import matplotlib.pyplot as plt

import hfGUIdata as hf
import plot_tools_jgb as pt
import squeeze_func_time as squ

props = [hf.brightMean, hf.darkMean, hf.det_t]

save = False
name = "depol_8_07.png"

# containers for data sets
ats=[]
Cs = []
Cerrs = []
Ns = []
J1ks = []
Ncals = []
names = []

base_path = os.getcwd()
add_path = "depolarization"
fns = [os.listdir(os.path.join(base_path,add_path))[i] for i in [1]]
J1ks = (475.0*3.03)*np.ones(np.shape(fns))
Ncals = 1.3999 * np.ones(np.shape(fns))  # #photons per ion per ms

#_____________________________________________________________________
# data processing here
for i,fn in enumerate(fns):
    folder = os.path.join(base_path,add_path,fn)
    print folder
    os.chdir(folder)
    files = os.listdir(os.getcwd())

 #Load properties data
    prop_name = [x for x in files if "_props.csv" in x][0]
    file_name, data_p = hf.get_gen_csv(prop_name, skip_header=False)
    bm = data_p['det_brightMean']
    dm = data_p["det_darkMean"]
    det_t = data_p["det_t"]
    int_t = 2e-6*data_p["squeeze_arm_t"]  #total interaction time in secs
    k = bm-dm  # phtns per N atoms
    N = k/(det_t*1e-3)/Ncals[i]

    # load experiment data
    data_name = [x for x in files if "_data.csv" in x][0]
    file_name, data = hf.get_gen_csv(data_name, skip_header=True)
    arm_time = np.array(data.T[0][0:],dtype='float')
    avg_pmt_counts = np.array(data.T[1][0:],dtype='float')
    pmterr = np.array(data.T[2][0:],dtype='float')
    trials = np.array(data.T[3][0:],dtype='float')

    b_prob = (avg_pmt_counts - dm)/(float(bm - dm))
    pmterr_of_mean = pmterr/sqrt(trials)
    b_prob_err = pmterr_of_mean/(float(bm - dm))
    contrast_est = 1-(2*b_prob)
    contrast_est_err = 2 * b_prob_err

    ats.append(arm_time)
    Cs.append(contrast_est)
    Cerrs.append(contrast_est_err)
    Ns.append(N)
    names.append(hf.n_slice(fn))
    os.chdir(base_path)

#%%
#________________________________________________________________________
# visualizing the experimental data
for i,data in enumerate(ats[0:3]):
    l = "N: {:.0f}, J: {:.0f}".format(Ns[i],J1ks[i])
    plt.errorbar(2e-3*ats[i],Cs[i],yerr=Cerrs[i],fmt='o',label=l)
plt.legend(loc=3, fontsize=10)
plt.xlabel("Interaction time [ms]")
plt.ylabel("Ramsey fringe contrast")


#________________________________________________________________________
#add some theory curves
mult = 2
G_el = mult* 60.56
G_ud = mult* 9.075
G_du = mult* 6.413
G_tot = 70.0  # per s

ti = np.linspace(1e-6,4.0e-3,num=100)  # seconds
spem = np.exp(-G_tot*ti)
plt.plot(ti*1e3, spem,'--k',label='Spon. Emiss.')

colors = ['k', pt.red, pt.blue]
for j,Jbar1k in enumerate(J1ks[0:3]):
    Jbar = Jbar1k/(0.002/ti)
    out = squ.OAT_decoh(0.0, ti, Jbar, Ns[j], G_el, G_ud, G_du)
    C_coherent_pred = np.real(out[1])
    #plt.plot(ti*1e3,C_coherent_pred,c=colors[j])
    plt.plot(ti*1e3,C_coherent_pred,c=colors[j])


plt.axis([0,4.0,0.0,1.1])
plt.grid('off')

if save is True:
    plt.savefig(name,dpi=300,bbox='tight',transparent=True)
