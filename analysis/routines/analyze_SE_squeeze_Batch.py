# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 14:42:25 2015

@author: jgb
"""

import os, shutil, importlib
import numpy as np
from numpy import pi, sqrt, sin
import matplotlib.pyplot as plt

import hfGUIdata as hf
import plot_style as ps
importlib.reload(ps)
import squeeze_func_time as squ

#options
verbose = True
save = False
img_name = "spinNoise_9_23"


#theory calc info
G_el =  67.10
G_ud =  10.07
G_du =  7.10
G_tot = 42.1
#adjust for extra decohrence
G_tot = 0.5*(67.10+17.1+38.6)
print(G_tot)
G_el = 67.10 + 38.6

# containers for data sets
psis=[]
its=[]
sig_obs = []
sig_ins = []
sig_pns = []
SE = []
Ns = []
names = []

base_path = os.getcwd()
fns = [os.listdir(base_path)[i] for i in [1]]
J1ks = (1940.0)*np.ones(np.shape(fns))
Ncals = 1.2 * np.ones(np.shape(fns))  # #photons per ion per ms

#%%
#_____________________________________________________________________
# data processing here
for i,fn in enumerate(fns):
    folder = os.path.join(base_path,fn)
    os.chdir(folder)
    files = os.listdir(os.getcwd())

    #Load properties data
    prop_name = [x for x in files if "_props.csv" in x][0]
    file_name, data_p = hf.get_gen_csv(prop_name, skip_header=False)
    bm = data_p['det_brightMean']
    dm = data_p["det_darkMean"]
    det_t = data_p["det_t"]
    int_t = 2e-6*data_p["squeeze_arm_t"]  #total interaction time in secs
    #reps = data_p["reps"]
    k = bm-dm  # phtns per N atoms
    N = k/(det_t*1e-3)/Ncals[i]

    #load calibration data
    cal_name = [x for x in files if "_cal.csv" in x][0]
    file_name, data = hf.get_gen_csv(cal_name, skip_header=True)
    cal_counts_avg = np.mean(data.T[1])
    sig_cal_ob = np.mean(data.T[2])
    sig_cal_sn = sqrt(cal_counts_avg)

    # load experiment data
    data_name = [x for x in files if "_data.csv" in x][0]
    file_name, data = hf.get_gen_csv(data_name, skip_header=True)
    reps = np.mean(data.T[3])
    psi_deg = data.T[0]
    count_avg = data.T[1]
    sig_sn = sqrt(data.T[1])
    sig_ob = data.T[2]
    sig_in = sqrt(sig_ob**2 - sig_sn**2)  # subtract poissonian shot noise
    sig_ob_err = sig_ob * 1/sqrt(2*reps)
    sig_in_err= sig_in * 1/sqrt(2*reps)

    sig_pn = sqrt(k**2/4.0/N)  # calclated from the atom number

    sig_a = sqrt(sig_cal_ob**2 - sig_cal_sn**2 - sig_pn**2)
    sig_a_deg = (sig_a/cal_counts_avg)*180.0/pi

    dB_squ_in = 10*np.log10((sig_in**2)/(sig_pn**2))
    dB_squ_ob = 10*np.log10((sig_ob**2)/(sig_pn**2))
    
    # calc reduction in constrast (from model shown to represent data)
    Jbar = J1ks[i] /(0.002/int_t)
    out = squ.OAT_decoh(0.0, int_t, Jbar, N, G_el, G_ud, G_du)
    C_coherent_pred = np.real(out[1])
    csi_R2 = (sig_ob**2)/(sig_pn**2)/C_coherent_pred**2
    dB_csi_R2 = 10*np.log10(csi_R2)
    

    #Load data messages
    print( "______________Data set: {}__________________".format(hf.n_slice(file_name)))
    print( "Ion number from loading calibration: {:.0f}".format(N))
    if verbose:
        print( "Meas. noise sig_cal_ob [Stdv PMT counts]: {:.3f}".format(sig_cal_ob))
        print( "Infer. proj. noise sqrt(cal_ob^2 - cal_sn^2):  {:.3f}".format(sqrt(sig_cal_ob**2 - sig_cal_sn**2)) )
        print( "Pred. proj. noise sig_pn = sqrt(k**2/4.0/N): {:.3f}".format(sig_pn) )
        print( "Est. added noise sqrt((cal_ob^2 - cal_sn^2)-pn^2): {:.3f}, or {:.3f} deg".format(sig_a,sig_a_deg))
    print( "Minimum inferred spin variance: {:.3f} dB".format(np.min(dB_squ_in))) 
    print( "Minimum observed spin variance: {:.3f} dB".format(np.min(dB_squ_ob)))
    print( r"Minimum observed $\csi_R^2$: {:.3f} dB".format(np.min(dB_csi_R2)))


    psis.append(psi_deg)
    its.append(int_t)
    sig_obs.append(sig_ob)
    sig_ins.append(sig_in)
    sig_pns.append(sig_pn)
    SE.append(np.max(10*np.log10(csi_R2**(-1))))
    Ns.append(N)
    names.append(hf.n_slice(file_name))
    os.chdir(base_path)

#%%
#________________________________________________________________________
# visualizing the experimental data

for i,data in enumerate(sig_obs):
    l = r"$\tau=$ {:.1f} ms, N: {:.0f}".format(its[i]*1e3,Ns[i])
    spin_noise = (sig_ins[i]**2)/(sig_pns[i]**2)
    spin_noise_dB = 10*np.log10(spin_noise)
    spin_noise_err_dB = 10*np.log10(spin_noise) - 10*np.log10(spin_noise-2*spin_noise/sqrt(2*reps))
    plt.errorbar(psis[i],spin_noise_dB,yerr=spin_noise_err_dB, fmt='o',label=l)

#plt.yscale('log')
#plt.xscale('log')
plt.axis([-1,181,-10,15])
plt.xlabel(r"Tomography angle $\psi$ [deg]",fontsize=14)
plt.ylabel("Spin variance [dB]",fontsize=14)
plt.grid('off')
plt.legend(loc=0,fontsize=10)


#________________________________________________________________________
#add some theory curves

psi = np.linspace(0.001,pi,num=100) # radians
cs = ['b','g']
for i,name in enumerate(names):
    Jbar = J1ks[i]/(0.002/its[i])
    out = squ.OAT_decoh(-psi, its[i], Jbar, Ns[i], G_el, G_ud, G_du)
    out_u = squ.OAT_decoh(-psi, its[i], Jbar, Ns[i]+5, G_el, G_ud, G_du)
    out_l = squ.OAT_decoh(-psi, its[i], Jbar, Ns[i]-5, G_el, G_ud, G_du)
    R = np.real(out[0]/(sqrt(Ns[i])/(2.0)))**2
    R_add = R + 1.0*sin(psi)
    R_dB = 10*np.log10(R)
    R_add_dB = 10*np.log10(R_add)
    plt.plot(psi*180/pi,R_dB,color=cs[i])
    plt.plot(psi*180/pi,R_add_dB,color=cs[i],linestyle='--')
    #plt.fill_between(ti*1e3,C_l,C_u,facecolor=colors[j],alpha=0.5)

if save is True:
    os.chdir('..')
    plt.savefig(img_name+".png",dpi=300,bbox='tight',transparent=True)
    # make a copy of the analysis at the folder
    shutil.copy(__file__, os.getcwd())
    os.chdir(base_path)

plt.show()
plt.close()
int_times = np.array(its)*1e3
plt.plot(int_times,SE,'o')
plt.ylabel('Spectroscopic Enhancement [dB]')
plt.xlabel('Interaction time [ms]')
plt.show()