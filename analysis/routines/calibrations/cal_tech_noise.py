# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 10:37:46 2015

@author: jgb
"""
import os, shutil, importlib, os.path
import numpy as np
from numpy import pi, sin, cos, sqrt
import matplotlib.pyplot as plt

import hfGUIdata as hf
import plot_style as ps
import plot_model_fit as pt
importlib.reload(pt)
import squeeze_func_time as squ

props = [hf.brightMean, hf.darkMean, hf.det_t]

#inputs
files_to_use = [-1]
Ncal = 3.672

verbose = False
save = False
name = "cal_tech_noise"
# make a copy of the analysis at the folder
if save is True:
    shutil.copy(__file__, os.path.normpath(os.getcwd()))

# containers for data sets
dts=[]
p_counts=[]
sig_pns = []
sig_obs = []
sig_tes = []
sig_ob_errs = []
sig_te_errs = []
Ns = []
Ncals = []
names = []

base_path = os.getcwd()
add_path = ""
fns = [os.listdir(os.path.join(base_path,add_path))[i] for i in files_to_use]
Ncals = Ncal * np.ones(np.shape(fns))  # #photons per ion per ms

#_____________________________________________________________________
# data processing here
for i,fn in enumerate(fns):
    folder = os.path.join(base_path,add_path,fn)
    print(folder)
    os.chdir(folder)
    files = os.listdir(os.getcwd())

 #Load properties data
    prop_name = [x for x in files if "_props.csv" in x][0]
    file_name, data_p = hf.get_gen_csv(prop_name, skip_header=False)
    bm = data_p['det_brightMean']
    dm = data_p["det_darkMean"]
    det_t = data_p["det_t"]
    k = bm-dm  # phtns per N atoms
    N = k/(det_t*1e-3)/Ncals[i]
    print("Number of ions: {:.4g}".format(N))

    # load experiment data
    data_name = [x for x in files if "_data.csv" in x][0]
    file_name, data = hf.get_gen_csv(data_name, skip_header=True)
    reps = np.mean(data.T[3])
    det_t = data.T[0]
    count_avg = data.T[1]
    sig_ob = data.T[2]
    sig_sn = sqrt(count_avg)

    sig_te = sqrt(sig_ob**2 - sig_sn**2)  # subtract poissonian shot noise
    sig_ob_err = sig_ob * 1/sqrt(2*reps)
    sig_te_err= sig_te * 1/sqrt(2*reps)
    
    sig_pn = sqrt(k**2/4.0/N)
    
    dts.append(det_t)
    p_counts.append(count_avg)
    sig_pns.append(sig_pn)
    sig_obs.append(sig_ob)
    sig_tes.append(sig_te)
    sig_ob_errs.append(sig_ob_err)
    sig_te_errs.append(sig_te_err)
    Ns.append(N)
    names.append((fn))

    os.chdir(base_path)

#%%
#________________________________________________________________________
# visualizing the experimental data
def added_noise(m, eps):
    sig_tot_2 = m + eps*m**2 
    return sig_tot_2
    

for i,data in enumerate(dts):
    #sort the data
    sort_ind = data.argsort()
    dts_s = data[sort_ind]
    sig_obs_s = sig_obs[i][sort_ind]
    sig_tes_s = sig_tes[i][sort_ind]
    sig_ob_errs_s = sig_ob_errs[i][sort_ind]
    p_counts_s = p_counts[i][sort_ind]
    guess=np.array([0.0000])
    hold=np.array([False])
    pl = ["Detected photons","Var of Photon number",'Extract Tech Noise']
    print(plt.xlabel)
    plt.plot(dts_s,sig_tes_s**2,'o')

    plt.plot(p_counts_s,sig_tes_s**2,'o') 
    plt.close()
    pout,perr = pt.plot_fit(p_counts_s,sig_obs_s**2,added_noise,guess,
                yerr=2*sig_obs_s*sig_ob_errs_s,
                hold=hold,
                labels=pl)

    plt.show()
    plt.close()

    p = np.array([0.0,1.0,0.0])
    pout,perr = pt.plot_polyfit(dts_s,p_counts_s,p,
                yerr=sig_obs_s,
                hold=np.array([True,False,True]),
                labels=["Det time [ms]","Photons","Check Depumping"])
    plt.show()
    plt.close()
    p = np.array([0.2,0.0,0.0])
    pout,perr = pt.plot_polyfit(dts_s,sig_tes_s**2/(p_counts_s),p,
                hold=np.array([False,False,True]),
                labels=["Det time [ms]",r"Tech Noise to PSN ratio $\sigma^2_{tec}$/m","Tech Noise ratio"])
    plt.show()
    plt.close()


    
#print("Added std dev (degrees): {0:.4g} t^2, {1:.4g} t^4)".format(pout_deg[0],pout_deg[1]))
