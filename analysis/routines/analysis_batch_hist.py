# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 20:36:13 2015

@author: jgb
"""

import os, shutil, importlib
import numpy as np
from numpy import pi, sqrt
import matplotlib.pyplot as plt
import scipy.stats.mstats as mstats

import hfGUIdata as hf
importlib.reload(hf)
import plot_style as ps
importlib.reload(ps)

#options
Ncal = 1.22
verbose = True
save = True
ymax = 350
files_to_use = [-4]
hist_to_use = [0,1,2,3,4,5,6]
text_name = "batch_hist_1016_wODF_tau3000.pdf"
img_name = "batch_hist_img_1016"
num_bins = 37#sqrt(len(z_data))
base_path = os.getcwd()
data_path = base_path
os.chdir(data_path)

# containers for data sets
psis=[]
its=[]
sig_obs = []
sig_ins = []
sig_pns = []
SE = []
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
    k = bm-dm  # phtns per N atoms
    N = k/(det_t*1e-3)/Ncals
    print(N)
    
    #load batch data
    data_name = [x for x in files if "_data.csv" in x][0]
    file_name, data = hf.get_gen_csv(data_name, skip_header=True)
    scan_data = np.array(data.T[0][0:],dtype='float')
    avg_pmt_counts = np.array(data.T[1][0:],dtype='float')

    
    hdata_to_use = np.array([hdata[i] for i in hist_to_use])    
    for i,row in enumerate(hdata_to_use):
        print("_________________________________________")
        det_array = np.copy(row)
        counts_data = hf.parse_raw_counts(det_array)
    
        Sz_data = 2*(((counts_data-dm)/(bm - dm)) - 0.5)
        datas.append(Sz_data)
    
        bs = np.arange(-1.01,1.01,(2.02/num_bins))

        lab = r"Scan data = {0:.4g}, Squeeze time = {1:.4g}".format(scan_data[i], int_time)
        lab = r"Scan value = {0:.4g}, Squeeze time = {1:.4g}".format(scan_data[hist_to_use[i]], int_time)

        plt.hist(datas[i],bs,label=lab,alpha=0.6)#, align='right')
        plt.axis([-1.1,1.1,0,ymax])
        plt.xlabel(r"Spin projection 2$S_\psi$/N")
        plt.ylabel("Trials")
        
        plt.show()
        plt.close()
        
        k2,pval = mstats.normaltest(datas[i])
        print(lab)
        print("# of trials: {}".format(len(datas[i])))
        print("Mean: {0:.3g}, Median: {1:.3g}, Min: {2:.3g}, Max {3:.3g}".format(np.mean(datas[i]),np.median(datas[i]),np.min(datas[i]),np.max(datas[i])))
        print("Std Dev: {0:3g}, Variance: {1:.3g}".format(np.std(datas[i]),np.var(datas[i])))        
        #print("Normality tests: skew+kurtosis: {0:.4g}, pval: {1:.4g}".format(k2,pval))

    os.chdir(data_path)
if save is True:    
    ps.save_data_txt(text_name+".txt",datas)

"""
#plt.legend(fontsize=11)
plt.axis([-1.0,1.0,0,100])
plt.xlabel(r"Transverse spin projection 2$S_\psi$/N")
plt.ylabel("Experiments")
if len(files_to_use) == 1:
    plt.title(r"$\tau$: {0:.3g} ms, N: {1:.0f}, $\psi$: {2:.3g}".format(tau,N,final_phase*180/pi))
        
    plt.show()
    plt.close()
    plt.plot(z_data,'o')
    #plt.axis([0,len(z_data),-1.1,1.1])
    plt.xlabel("Trial number")
    plt.ylabel(r"Transverse spin projection 2$S_\psi$/N" )


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
"""