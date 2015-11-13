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
import resample_tools as re

#options
colors = ['k', ps.red, ps.blue, ps.orange, ps.pink]
raw = False
verbose = True
save = False

img_name = "spinNoise_11_10"
files_to_use = [-1]
J1k = 1500   
Ncal = 1.70

#theory calc info
G_el =  67.4
G_ud =  10.1
G_du =  7.1

#adjust for extra decohrence
G_add = 60.0
G_tot = 0.5*(G_el + (G_ud+G_du) + G_add)
print(G_tot)
G_el = G_el + G_add

#added noise from Jy noise fit
A = 0.00068 # rad^2/ms^2
B = 0.0002323 # rad^2/ms^4

# containers for data sets
psis=[]
its=[]
sig_obs = []
sig_ob_errs = []
sig_robs = []
sig_rob_errs = []
sig_ins = []
sig_2_in_err = []
sig_pns = []
sig_psns = []
SE = []
Ns = []
names = []

base_path = os.getcwd()
fns = [os.listdir(base_path)[i] for i in files_to_use]
J1ks = J1k*np.ones(np.shape(fns))
Ncals = Ncal * np.ones(np.shape(fns))  # #photons per ion per ms

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
    sig_sn_2_err = sig_sn**2 * sqrt(2/reps)
    sig_ob = data.T[2]
    sig_in = sqrt(sig_ob**2 - sig_sn**2)  # subtract poissonian shot noise
    sig_ob_2_err = sig_ob**2 * sqrt(2/reps)
    sig_in_2_err= sqrt( (sig_ob_2_err)**2  + (sig_sn_2_err)**2)

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
    
    if raw is True:
    # Load histgram
        files = os.listdir(os.getcwd())
        data_name = [x for x in files if "_raw.csv" in x][0]
        hdata = np.genfromtxt(data_name, delimiter=",", dtype='float')
        print(np.shape(hdata))
        
        jack_sig = np.zeros(np.size(sig_ob))
        jack_err_on_sig = np.zeros(np.size(sig_ob))
        for i,row in enumerate(hdata):
            counts = hf.parse_raw_counts(row)
            jack_sig[i], jack_err_on_sig[i] = re.jackknife_est(counts,np.std)
            
        sig_robs.append(jack_sig)
        sig_rob_errs.append(jack_err_on_sig)


    #Load data messages
    print( "______________Data set: {}__________________".format(hf.n_slice(file_name)))
    print( "Ion number from loading calibration: {:.0f}".format(N))
    if verbose:
        print( "Meas. noise sig_cal_ob [Stdv PMT counts]: {:.3f}".format(sig_cal_ob))
        print( "Infer. proj. noise sqrt(cal_ob^2 - cal_sn^2):  {:.3f}".format(sqrt(sig_cal_ob**2 - sig_cal_sn**2)) )
        print( "Pred. proj. noise sig_pn = sqrt(k**2/4.0/N): {:.3f}".format(sig_pn) )
        print( "Est. added noise sqrt((cal_ob^2 - cal_sn^2)-pn^2): {:.3f}, or {:.3f} deg".format(sig_a,sig_a_deg))
        print( "Est. added noise (ratio of var to photon shot noise: {:.3f}".format(sig_a**2/np.mean(sig_sn)**2))
    print( "Minimum inferred spin variance: {:.3f} dB".format(np.min(dB_squ_in))) 
    print( "Minimum observed spin variance: {:.3f} dB".format(np.min(dB_squ_ob)))
    print( r"Minimum observed $\csi_R^2$: {:.3f} dB".format(np.min(dB_csi_R2)))


    psis.append(psi_deg)
    its.append(int_t)
    sig_obs.append(sig_ob)
    sig_ob_errs.append(sqrt(sig_ob_2_err))
    sig_ins.append(sig_in)
    sig_2_in_err.append(sig_in_2_err)
    sig_pns.append(sig_pn)
    sig_psns.append(sig_sn)
    SE.append(np.max(10*np.log10(csi_R2**(-1))))
    Ns.append(N)
    names.append(hf.n_slice(file_name))
    os.chdir(base_path)

#%%
#________________________________________________________________________
# visualizing the experimental data

for i,data in enumerate(sig_obs):
    l = r"$\tau=$ {:.3g} ms, N: {:.0f}".format(its[i]*1e3,Ns[i])
    spin_noise = (sig_ins[i]**2)/(sig_pns[i]**2)
    spin_noise_err = sig_2_in_err[i]/(sig_pns[i]**2)
    spin_noise_err = sqrt( (sig_2_in_err[i]/sig_pns[i]**2)**2 + ((sig_ins[i]/sig_pns[i])**2 * 0.05)**2 ) # accounting for 5% uncertainty in PN
    spin_noise_dB = 10*np.log10(spin_noise)
    #spin_noise_err_dB = 10*np.log10(spin_noise) - 10*np.log10(spin_noise-2*spin_noise/sqrt(2*reps))
    spin_noise_err_dB = 10*np.log10(spin_noise) - 10*np.log10(spin_noise-spin_noise_err)
    plt.errorbar(psis[i],spin_noise_dB,yerr=spin_noise_err_dB, fmt='o',label=l, color=colors[i])

#plt.yscale('log')
#plt.xscale('log')
plt.axis([-1,181,-10,17])
plt.xlabel(r"Tomography angle $\psi$ [deg]",fontsize=14)
plt.ylabel("Spin variance [dB]",fontsize=14)
plt.grid('off')

#________________________________________________________________________
#add some theory curves

psi = np.linspace(0.001,pi,num=100) # radians
for i,name in enumerate(names):
    Jbar = J1ks[i]/(0.002/its[i])
    out = squ.OAT_decoh(-psi, its[i], Jbar, Ns[i], G_el, G_ud, G_du)
    out_u = squ.OAT_decoh(-psi, its[i], Jbar, Ns[i]+5, G_el, G_ud, G_du)
    out_l = squ.OAT_decoh(-psi, its[i], Jbar, Ns[i]-5, G_el, G_ud, G_du)
    C = np.exp(-G_tot*its[i])  # reduction in contrast from spontaneous emission
    R = np.real(out[0]/(sqrt(Ns[i])/(2.0)))**2
    R_add = R + (A*(its[i]*1e3)**2)*(N*C**2) * sin(psi)**2 + (B*(its[i]*1e3)**4)*(N*C**2) * sin(psi)**2
    R_dB = 10*np.log10(R) 
    R_add_dB = 10*np.log10(R_add)
    plt.plot(psi*180/pi,R_dB,color=colors[i])
    plt.plot(psi*180/pi,R_add_dB,color=colors[i],linestyle='--')
    
    #where is the limit just due to technical noise?
    tech_limit = 10*np.log10((0.3*sig_psns[i]**2)/(sig_pns[i]**2))
    plt.plot(psis[i],tech_limit,color=colors[i],linestyle='--')
    
    #plt.fill_between(ti*1e3,C_l,C_u,facecolor=colors[j],alpha=0.5)
    print("added dephasing: {:.3g} (ratio of var to proj noise)".format((A*(its[i]*1e3)**2)*N))

plt.legend(loc=0,fontsize=10)
if len(names) is 1:
    plt.title(names[0])

if save is True:
    os.chdir('..')
    plt.savefig(img_name+".png",dpi=300,bbox='tight',transparent=True)
    # make a copy of the analysis at the folder
    shutil.copy(__file__, os.getcwd())
    os.chdir(base_path)

plt.show()
plt.close()

if verbose is True:

    int_times = np.array(its)*1e3
    plt.plot(int_times,SE,'o')
    plt.ylabel('Spectroscopic Enhancement [dB]')
    plt.xlabel('Interaction time [ms]')
    plt.show()

    #comparing different error estimates
    if raw is True:
        plt.plot(sig_ob_errs[0], sig_rob_errs[0])

