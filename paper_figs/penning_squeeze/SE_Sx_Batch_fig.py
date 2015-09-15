# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 10:37:46 2015

@author: jgb
"""
import os, importlib, shutil
import numpy as np
from numpy import pi, sin, cos, sqrt
import matplotlib.pyplot as plt

import hfGUIdata as hf
importlib.reload(hf)
import plot_style as ps
import squeeze_func_time as squ

props = [hf.brightMean, hf.darkMean, hf.det_t]

raw = False
save = True
name = "depolarizationVsN_noJ_fig.pdf"

# containers for data sets
ats=[]
Cs = []
Cerrs = []
Ns = []
J1ks = []
names = []
hist = []

base_path = os.getcwd()
#for the figure creation, point to the file folders that have the data sets we need
fns = ["/Volumes/688/Public/penning_britton/dailyLabBookFiles/2015/20150811/Load306/depolarization/2015-08-11--19.53.00.339",
       "/Volumes/688/Public/penning_britton/dailyLabBookFiles/2015/20150807/Load304/depolarization/2015-08-07--15.10.19.547",
       "/Volumes/688/Public/penning_britton/dailyLabBookFiles/2015/20150820/Load308/depolarization/2015-08-20--19.04.50.675"]
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
    Ncal = data_p['Ncal']
    J1k = data_p['J1k']
    
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
        hist.append(hdata)

    J1ks.append(J1k)
    ats.append(arm_time)
    Cs.append(contrast_est)
    Cerrs.append(contrast_est_err)
    Ns.append(N)
    names.append(hf.n_slice(fn))

    os.chdir(base_path)

#Have to get a data set by hand from this day, since analysis is different
"""
fn_8_6 = "/Volumes/688/Public/penning_britton/dailyLabBookFiles/2015/20150806/depolarization/phaseflop_datasets_8_06_save.csv"
data_8_6 = np.genfromtxt(fn_8_6, delimiter=",", names=True, dtype=None)
J1ks.append(472.0*3.03)
ats.append(data_8_6['tau_ms']/2.0*1e3)
Cs.append(data_8_6['C'])
Cerrs.append(data_8_6['C_err'])
Ns.append(97)
names.append("phaseflop_datasets_8_06")
"""

"""
fn_7_22 = "/Volumes/688/Public/penning_britton/dailyLabBookFiles/2015/20150722/Load296/depolarization/phaseflop_datasets_7_22_L296.csv"
data_7_22 = np.genfromtxt(fn_7_22, delimiter=",", names=True, dtype=None)
J1ks.append(420*3.03)
ats.append(data_7_22['tau_ms']/2.0*1e3)
Cs.append(data_7_22['C'])
Cerrs.append(data_7_22['C_err'])
Ns.append(130)
names.append("phaseflop_datasets_7_22_L296")
"""

#%%
#________________________________________________________________________
# visualizing the experimental data
#fig = plt.figure(figsize=(6.0,4.5))
for i,data in enumerate(ats):
    l = "N: {:.0f}".format(Ns[i])
    plt.errorbar(2e-3*ats[i],Cs[i],yerr=Cerrs[i],fmt='o',label=l)
plt.legend(loc=3, fontsize=10)
plt.xlabel(r"Interaction time  $\tau$ [ms]")
plt.ylabel(r"Spin Coherence  $ 2\left \langle S_x \right \rangle /N$")


#________________________________________________________________________
#add some theory curves
G_el =  57.1 * 1.0
G_ud =  8.56157395
G_du =  6.04795106
G_tot = 35.9  # per s

ti = np.linspace(1e-6,4.0e-3,num=100)  # seconds
spem = np.exp(-G_tot*ti)
plt.plot(ti*1e3, spem,'--k',label='Spon. Emiss.')

colors = ['k', ps.red, ps.blue, ps.purple]
for j,Jbar1k in enumerate(J1ks):
    Jbar = Jbar1k/(0.002/ti)
    out = squ.OAT_decoh(0.0, ti, Jbar, Ns[j], G_el, G_ud, G_du)
    C_coherent_pred = np.real(out[1])
    #plt.plot(ti*1e3,C_coherent_pred,c=colors[j])
    plt.plot(ti*1e3,C_coherent_pred,color=colors[j])


plt.axis([0,4.0,0.0,1.1])
plt.grid('off')

if save is True:
    plt.savefig(name,dpi=300,bbox='tight',transparent=True)
    
#%%
#play with histdata
#plt.close()
#
#def parse_raw_counts(array):
#    bad = 0
#    for x in np.nditer(array, op_flags=['readwrite']):
#        if x == -1:
#            print('Found bad data point')
#            bad += 1
#            x[...] = -1
#        else:
#            x[...] = int(x) & 0x1fff
#    if bad > 0:
#        print("# of bad points: {}".format(bad))
#
#det_array = np.copy(hist[0])
#parse_raw_counts(det_array)
#histim = []
#bins = range(0,450,15)
##for row in [det_array[i] for i in [17]]:
#for row in det_array:
#    h, bin_edges = np.histogram(row, density=False, bins=bins)
#    histim.append(h)
#    #plt.hist(row,bins)
#plt.imshow(np.transpose(histim), aspect='auto', origin='lower')
#
#pt.save_data_txt('histdata.txt', histim)


