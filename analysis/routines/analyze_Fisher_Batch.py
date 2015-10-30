# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 14:43:24 2015

@author: rabi
"""

import os, shutil, importlib
import numpy as np
from numpy import pi, sqrt, sin
import matplotlib.pyplot as plt

import hfGUIdata as hf
import plot_style as ps
import plot_model_fit as pt

#run script in the with "Fisher_batch# as the working directory
#inputs for loading data and histograms

files_to_use = [2]
Ncal = 1.32

bin_width = 2  # found from Strobel this was optimum
h = 20  #block size for resampling
base_path = os.getcwd()
data_path = base_path
os.chdir(data_path)

# containers for data sets
tipping_angles=[]
its = []
Ns = []
names = []
datas=[]
fns = [os.listdir(data_path)[i] for i in files_to_use]
Ncals = Ncal * np.ones(np.shape(fns))  # #photons per ion per ms

for i,fn in enumerate(fns):
    print("_________________________________________")
    folder = os.path.join(data_path,fn)
    os.chdir(folder)
    files = os.listdir(os.getcwd())
    print(folder)
   
    # Load histgram data
    files = os.listdir(os.getcwd())
    data_name = [x for x in files if "_raw.csv" in x][0]
    hdata = np.genfromtxt(data_name, delimiter=",", dtype='float')
    print(np.shape(hdata))
    
    #Load properties data
    prop_name = [x for x in files if "_props.csv" in x][0]
    file_name, data_p = hf.get_gen_csv(prop_name, skip_header=False)
    bm = data_p['det_brightMean']
    dm = data_p["det_darkMean"]
    det_t = data_p["det_t"]
    int_time = 2*data_p["squeeze_arm_t"]
    tpi = data_p["sf_fitParam_tpi"]
    k = bm-dm  # phtns per N atoms
    N = k/(det_t*1e-3)/Ncals
    print(N)
    
    #load batch data
    data_name = [x for x in files if "_data.csv" in x][0]
    file_name, data = hf.get_gen_csv(data_name, skip_header=True)
    scan_data = np.array(data.T[0][0:],dtype='float')
    avg_pmt_counts = np.array(data.T[1][0:],dtype='float')
    tip_angle_rad = scan_data*np.pi/180.0
    
    print("_________________________________________")
    print("File name: "+file_name)
    lab = r"Squeeze time = {0:.4g}".format(int_time)
    print(lab)
 
    # parse histogram raw data
    for i,row in enumerate(hdata):
        det_array = np.copy(row)
        counts_data = hf.parse_raw_counts(det_array)
   
        Sz_data = 2*(((counts_data-dm)/(bm - dm)) - 0.5) * (N/2.)
        datas.append(Sz_data)

    its.append(int_time)
    Ns.append(N)
    tipping_angles.append(tip_angle_rad)
    os.chdir(data_path)

#define scale
bin_def = np.arange(-N/2.0,N/2.0,bin_width)  

#calculate reference histogram
data_ref = datas[2]
vals_ref, bins, patches = plt.hist(data_ref,bin_def,normed=True)

"""
#remove the reference data from the rest
datas = datas[1:]
"""

#containers
data_hist = []
data_hist_jacks = []
data_hist_jack_errs = []
n_bin_fills = []
samps = []
for row in datas:
    vals = plt.hist(row, bin_def,normed=True)[0]
    h_dist = np.sum((sqrt(vals_ref) - sqrt(vals))**2)
    samp = np.size(row)
    n_bin_fill = np.size(vals[vals!=0])  #number of discrete bins for which P!=0
    
    #resampling    
    g = int(samp/h)
    jack_trials = np.reshape(row, (g,h))
    re_hds = []
    for i in range(0,g):
        re_row = np.reshape(np.concatenate((jack_trials[:i],jack_trials[i+1:])), (samp-h))
        #vals = plt.hist(re_row, bin_def,normed=True)[0]
        re_vals = plt.hist(re_row, bin_def,normed=True)[0]
        re_h_dist = np.sum((sqrt(vals_ref) - sqrt(re_vals))**2)
        re_hds.append(re_h_dist)
    h_dist_jack = g*h_dist - ((g-1)/float(g)* np.sum(re_hds))
    h_dist_jack_err = np.sqrt(np.var(re_hds))

    
    data_hist.append(h_dist)
    data_hist_jacks.append(h_dist_jack)
    data_hist_jack_errs.append(h_dist_jack_err)
    n_bin_fills.append(n_bin_fill)
    samps.append(samp)

l=['Tip angle (rad)', 'dh^2','Extract Fisher Info']

to_use = 12

#note: tipping angle starts at 1 because 0 corresponds to the reference
fit_res = pt.plot_polyfit(tipping_angles[0],data_hist_jacks,np.array([0,0,1]),
                          yerr = h_dist_jack_err,
                          hold=np.array([False,True,False]),
                          labels=l)
k2 = fit_res[0][1]
NormF = k2/N[0]/( (1/8.0) ) #dont' need the 1/M term from teh estimate
print("F = 8*k2, and then F/N is {:.3g}".format(NormF))
print("number of samples is {:.3g}".format(np.mean(samps)))


nonCx = np.linspace(0,0.8,num=200)

nonCy = nonCx**2 * N[0]/8.0
plt.show()
plt.close(0)
"""
for row in diff_hist:
    plt.plot(bins[:-1],row)

plt.show()
plt.close()
"""
plt.plot(nonCx,nonCy,'b-')
plt.errorbar(tipping_angles[0], data_hist_jacks,yerr=data_hist_jack_errs,fmt='o')
plt.title("Data wrt entanglement witness")