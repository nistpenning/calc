# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 10:37:46 2015

@author: jgb
"""
import os, shutil, os.path
import numpy as np
from numpy import pi, sin, cos, sqrt
import matplotlib.pyplot as plt

import hfGUIdata as hf
import plot_style as ps
import plot_model_fit as pt
import squeeze_func_time as squ

props = [hf.brightMean, hf.darkMean, hf.det_t]

#inputs
Ncal = 1.034
files_to_use = [-1]

verbose = False
save = False
name = "Jy_noise_batch"
# make a copy of the analysis at the folder
if save is True:
    shutil.copy(__file__, os.path.normpath(os.getcwd()))

# containers for data sets
ats=[]
p_counts=[]
sig_pns = []
sig_obs = []
sig_dephs = []
sig_ob_errs = []
sig_deph_errs = []
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

    # load experiment data
    data_name = [x for x in files if "_data.csv" in x][0]
    file_name, data = hf.get_gen_csv(data_name, skip_header=True)
    reps = np.mean(data.T[3])
    arm_t = data.T[0]
    count_avg = data.T[1]
    sig_sn = sqrt(count_avg)
    sig_ob = data.T[2]
    sig_in = sqrt(sig_ob**2 - sig_sn**2)  # subtract poissonian shot noise
    sig_ob_err = sig_ob * 1/sqrt(2*reps)
    sig_in_err= sig_in * 1/sqrt(2*reps)
    
    sig_pn = sqrt(k**2/4.0/N)
    
    ats.append(arm_t)
    p_counts.append(count_avg)
    sig_pns.append(sig_pn)
    sig_obs.append(sig_ob)
    sig_dephs.append(sig_in)
    sig_ob_errs.append(sig_ob_err)
    sig_deph_errs.append(sig_in_err)
    Ns.append(N)
    names.append((fn))

    os.chdir(base_path)

#%%
#________________________________________________________________________
# visualizing the experimental data

if verbose is True:
    for i,data in enumerate(ats[0:3]):
        l = "N: {:.0f}".format(Ns[i])
        plt.errorbar(2e-3*ats[i],sig_obs[i],yerr=sig_ob_errs[i],fmt='o',label=l)
        plt.plot(2e-3*ats[i],sig_pns[i]*np.ones(np.shape(ats[i])), label="Projection_noise")
    plt.legend(loc=0, fontsize=10)
    plt.xlabel("Interaction time [ms]",fontsize=14)
    #plt.ylabel(r"Avgerage spin 2$S_x$/N")
    plt.ylabel("Std Dev. photon count",fontsize=14)
    if len(names) is 1:
        plt.title(names[0])
    
    plt.show()
    plt.close()
    
    for i,data in enumerate(ats[0:3]):
        l = "N: {:.0f}, Proj. noise is {:.4g}".format(Ns[i],180/sqrt(Ns[i])/pi)
        sig_add = sqrt(sig_dephs[i]**2 - (sig_pns[i]*np.ones(np.shape(ats[i])))**2)
        plt.plot(2e-3*ats[i],(sig_add/p_counts[i])*180.0/pi,'o',label=l)
        plt.plot(2e-3*ats[i],180/sqrt(Ns[i])/pi*np.ones(np.shape(ats[i])), label="Projection noise")
    plt.legend(loc=0, fontsize=10)
    plt.xlabel("Interaction time [ms]",fontsize=14)
    #plt.ylabel(r"Avgerage spin 2$S_x$/N")
    plt.ylabel("Added noise, std dev [deg]",fontsize=14)
    if len(names) is 1:
        plt.title("Added noise alone")
        
    plt.show()
    plt.close()

def added_noise(tau, N, m, k, A, B):
    sig_tot_2 = m + (k**2/4.0/N) + A*tau**2 + B*tau**4
    return sqrt(sig_tot_2)

for i,data in enumerate(ats[0:3]):
    guess=np.array([Ns[i],np.mean(p_counts[i]),k,0.1,0.0])
    hold=np.array([True,True,True,False,False])
    pl = ["Interaction time [ms]","Std Dev. photon count",'Extract added']
    pout,perr = pt.plot_fit(2e-3*ats[i],(sig_obs[i]),added_noise,guess,
                yerr=sig_ob_errs[i],
                hold=hold,
                labels=pl)
    pout_deg = sqrt(pout)/np.mean(p_counts[i])*180/pi
    plt.show()
print("Added std dev (degrees): {0:.4g} t^2, {1:.4g} t^4)".format(pout_deg[0],pout_deg[1]))
                
"""
#________________________________________________________________________
#add some theory curves
ti = np.linspace(1e-6,4.0e-3,num=100)  # seconds
spem = np.exp(-G_tot*ti)
plt.plot(ti*1e3, spem,'--k',label='Spon. Emiss.')

colors = ['k', ps.red, ps.blue]
for j,Jbar1k in enumerate(J1ks[0:3]):
    Jbar = Jbar1k/(0.002/ti)
    out = squ.OAT_decoh(0.0, ti, Jbar, Ns[j], G_el, G_ud, G_du)
    C_coherent_pred = np.real(out[1])
    #plt.plot(ti*1e3,C_coherent_pred,c=colors[j])
    plt.plot(ti*1e3,C_coherent_pred,color=colors[j])


plt.axis([0,4.0,0.0,1.1])
plt.grid('off')

if save is True:
    plt.savefig(name,dpi=300,bbox='tight',transparent=True)
"""