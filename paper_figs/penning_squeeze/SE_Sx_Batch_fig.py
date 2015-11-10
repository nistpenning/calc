# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 10:37:46 2015

@author: jgb
"""
import os, importlib, shutil
import numpy as np
from numpy import pi, sin, cos, sqrt
import matplotlib.pyplot as plt

import resample_tools as re
import hfGUIdata as hf
importlib.reload(hf)
import plot_style as ps
import plot_model_fit as pt
importlib.reload(pt)
import squeeze_func_time as squ

props = [hf.brightMean, hf.darkMean, hf.det_t]

raw = True
save = False
save_txt = False
name = "SxVsN_fig_v2.pdf"
colors = ['k', ps.red, ps.blue, ps.purple]
shapes = ['o','s','D','^']

# containers for data sets
ats=[]
Cs = []
Craws = []
Cerrs = []
Crerrs = []
Ns = []
#J1ks = [1440,1776,2316,1785]
J1ks = [1440,1776,2316]

#Ncals = [1.7425,1.1,0.85,0.8]
Ncals = [1.7425,1.1,0.85]

#G_add = np.array([120.0,37.6,46.6,50.0])
G_add = np.array([120.0,37.6,46.6])
names = []
hist = []

base_path = os.getcwd()
#for the figure creation, point to the file folders that have the data sets we need
fns = ["/Users/jgb/Data/20150811/Load306/depolarization/2015-08-11--19.53.00.339",
       "/Users/jgb/Data/20151001/Load333/Sx/2015-10-01--11.42.37.735",
       "/Users/jgb/Data/20150918/depolarization/AtDecouplings/2015-09-18--17.48.36.467"]
       #"/Users/jgb/Data/20150929/Load330/depol/2015-09-29--17.34.44.654"]
#store the parameters for the N value in the props file
#J1ks = (475.24*3.03)*np.ones(np.shape(fns)) # per sec at 1 kHz detuning
#Ncals = 1.4924 * np.ones(np.shape(fns))  # #photons per ion per ms

#_____________________________________________________________________
# data processing here
for i,fn in enumerate(fns):
    os.chdir(fn)
    print(fn)
    files = os.listdir(os.getcwd())

 #Load properties data
    prop_name = [x for x in files if "_props.csv" in x][0]
    file_name, data_p = hf.get_gen_csv(prop_name, skip_header=False)
    bm = data_p['det_brightMean']
    dm = data_p["det_darkMean"]
    det_t = data_p["det_t"]
    int_t = 2e-6*data_p["squeeze_arm_t"]  #total interaction time in secs
    Ncal = Ncals[i] #data_p['Ncal']
    J1k = J1ks[i] # data_p['J1k']
    
    k = bm-dm  # phtns per N atoms
    N = k/(det_t*1e-3)/Ncal

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
    
    # Load histgram
    if raw is True:
        data_name = [x for x in files if "_raw.csv" in x][0]
        hdata = np.genfromtxt(data_name, delimiter=",", dtype='float')
        print(np.shape(hdata))
        
        raw_mean = np.zeros(np.size(hdata, axis=0))
        raw_mean_err = np.zeros(np.size(hdata, axis=0))
        for i,row in enumerate(hdata):
            det_array = np.copy(row)
            counts_data = hf.parse_raw_counts(det_array)
       
            Sz_data = (1-2*((counts_data-dm)/k))
            raw_mean[i],raw_mean_err[i] = re.jackknife_est(Sz_data, np.mean)
            

    ats.append(arm_time)
    Cs.append(contrast_est)
    Craws.append(raw_mean)
    Cerrs.append(contrast_est_err)
    Crerrs.append(raw_mean_err)
    Ns.append(N)
    names.append(hf.n_slice(fn))

    os.chdir(base_path)

#%%
#________________________________________________________________________
# visualizing the experimental data
#fig = plt.figure(figsize=(6.0,4.5))
data_for_save = []
for i,data in enumerate(ats):
    l = "N: {:.0f}".format(Ns[i])
    Jbar = J1ks[i]/(0.002/(data*2e-6))
    Jt_opt = (24**(1/6.)*(Ns[i]/2)**(-2/3.))/4.*Ns[i]
    Jt = (Jbar*data*2e-6) #/ Jt_opt
    plt.errorbar(data*2e-3,Cs[i],yerr=Cerrs[i],fmt='o',label=l,marker=shapes[i],color=colors[i])
    if save_txt is True:
        data_for_save.append(2e-3*ats[i])
        data_for_save.append(Cs[i])
        data_for_save.append(Cerrs[i])
if save_txt is True:
    names = "t_ms_21, Sx_21, Sx_err_21, t_ms_66, Sx_66, Sx_err_66, t_ms_100, Sx_100, Sx_err_100"
    pt.save_data_txt("Sx_data.txt",data_for_save, col_names=names)
        
#plt.legend(loc=3, fontsize=10)
plt.xlabel(r"Interaction time $\tau$ (ms)")
plt.ylabel(r"Contrast  2$|\left \langle \vec{S} \right \rangle|$/N")


#________________________________________________________________________
#add some theory curves
#theory calc info
G_el =  61.6
G_ud =  9.24
G_du =  6.52
G_tot = 38.7
#adjust for extra decohrence

G_tot = 0.5*(G_el + (G_ud+G_du) + G_add)
print(G_tot)
G_els = G_el + G_add


ti = np.linspace(1e-6,4.0e-3,num=100)  # seconds
spem = np.exp(-60.0*ti)
Jbar_theory = np.mean(J1ks)/(0.002/ti)
plt.plot(ti*1e3, spem,'--k',label='Spon. Emiss.')

for j,Jbar1k in enumerate(J1ks):
    print("Total scattering rate: {0:0.3g}".format(0.5*(G_els[j] + (G_ud+G_du) )))
    Jbar = Jbar1k/(0.002/ti)
    Jt_opt = (24**(1/6.)*(Ns[j]/2)**(-2/3.))/4.*Ns[j]
    Jt = Jbar*ti # / Jt_opt
    out = squ.OAT_decoh(0.0, ti, Jbar, Ns[j], G_els[j], G_ud, G_du)
    C_coherent_pred = np.real(out[1])
    #plt.plot(ti*1e3,C_coherent_pred,c=colors[j])
    plt.plot(ti*1e3,C_coherent_pred,color=colors[j])


plt.axis([0.,4.,0.0,1.05])
plt.grid('off')
#plt.xscale('log')

if True:
    # get decohernce data to use as comparion to depolarization
    os.chdir("/Users/jgb/Data/20151001/Load333/decoh/2015-10-01--11.45.55.754")
    fn, data = hf.get_gen_csv('phase', skip_header=True)
    arm_time = np.array(data.T[0][1:],dtype='float')
    t0_ind = np.where(arm_time==10.0)
    int_time = arm_time*(2e-3)
    A = np.array(data.T[1][1:],dtype='float')
    B = np.array(data.T[2][1:],dtype='float')
    contrast = ((A - B)/(A[t0_ind]-B[t0_ind]))
    #use Jbar from that data set
    Jbar = 1776.0/(0.002/(int_time*1e-3))
    Jt = Jbar*int_time*1e-3
    plt.plot(int_time,contrast,'s',color='k',zorder = 1)
    os.chdir(base_path)

if save is True:
    plt.savefig(name,dpi=300,bbox='tight',transparent=True)

plt.show()

