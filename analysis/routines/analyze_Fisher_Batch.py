# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 14:43:24 2015

@author: rabi
"""

import os, shutil, importlib
import numpy as np
import numpy.random as rand
from numpy import pi, sqrt, sin
import matplotlib.pyplot as plt

import hfGUIdata as hf
import plot_style as ps
import plot_model_fit as pt
importlib.reload(pt)

files_to_use = [3]
h_to_use = [0,1,2,3,4,5,6,7]
Ncal = 1.22
num_bins = 53.0 #sqrt(len(z_data))
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
   
    # Load histgram
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
    lab = r"Squeeze time = {0:.4g}".format(int_time)
    print(lab)
 
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
bin_def = np.arange(-N/2.0,N/2.0,(N/num_bins))

data_ref = datas[0]
vals_ref, bins, patches = plt.hist(data_ref,bin_def,normed=True)

datas = datas[1:]

data_hist = []
for row in datas:
    vals = plt.hist(row, bin_def,normed=True)[0]
    data_hist.append(vals)
data_hist = np.array(data_hist)
diff_hist = np.array([(sqrt(vals_ref) - sqrt(row)) for i,row in enumerate(data_hist) ])
sq_hist = diff_hist**2
h_dist = np.array([np.sum(row) for row in sq_hist])
plt.close()

l=['Tip angle (rad)', 'dh^2','Extract Fisher Info']
to_use = 8
fit_res = pt.plot_polyfit(tipping_angles[0][1:to_use], h_dist[0:(to_use-1)],np.array([0,0,1]),
                          hold=np.array([True,True,False]),
                          labels=l)
nonCx = np.linspace(0,0.02,num=100)
nonCy = nonCx**2 * N[0]/8.0

plt.show()
plt.close(0)

for row in diff_hist:
    plt.plot(bins[:-1],row)
print("F = 8*k2, and then F/N is {:.3g}".format(fit_res[0][0]*8/N[0]) )
plt.show()
plt.close()

plt.plot(nonCx,nonCy,'b-')
plt.plot(tipping_angles[0][1:to_use], h_dist[0:(to_use-1)],'o')
plt
