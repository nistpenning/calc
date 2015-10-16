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

def parse_raw_counts(array):
    bad = 0
    for x in np.nditer(array, op_flags=['readwrite']):
        if x == -1:
            print('Found bad data point')
            bad += 1
            x[...] = -1
        elif np.isnan(x) == True:
            print('Found bad data point Nan')
            bad += 1
            x[...] = -1
        else:
            x[...] = int(x) & 0x1fff
    if bad > 0:
        print("# of bad points: {}".format(bad))
        print("removing all bad points in return value")
        array = array[array!=-1]
    return array
    

#options
verbose = True
save = False
files_to_use = [-1]
img_name = "transverse_hist.pdf"
base_path = os.getcwd()
data_path = '/Users/jgb/Data/20151002/FisherW/2015-10-02--17.15.08.249'
os.chdir(data_path)

walsh_num = 4  # use 4 for a walsh seq, 2 for regular (mulitplies the arm time)
num_bins = 40#sqrt(len(z_data))

# containers for data sets
psis=[]
its=[]
datas=[]
Ncals = [1.1]
names = []

#fns = [os.listdir(data_path)[i] for i in files_to_use]

#for i,fn in enumerate(fns):
#folder = os.path.join(data_path,fn)

os.chdir(data_path)
files = os.listdir(os.getcwd())
#Load properties data
prop_name = [x for x in files if "_props.csv" in x][0]
file_name, data_p = hf.get_gen_csv(prop_name, skip_header=False)
bm = data_p['det_brightMean']
dm = data_p["det_darkMean"]
det_t = data_p["det_t"]
int_t = walsh_num * 1e-6*data_p["squeeze_arm_t"]  #total interaction time in secs
k = bm-dm  # phtns per N atoms
N = k/(det_t*1e-3)/Ncals[0]

# load experiment data
data_name = [x for x in files if "_data.csv" in x][0]
file_name, data = hf.get_gen_csv(data_name, skip_header=True)
tip_angle_deg = np.array(data.T[0][0:],dtype='float')
avg_pmt_counts = np.array(data.T[1][0:],dtype='float')
pmterr = np.array(data.T[2][0:],dtype='float')
trials = np.array(data.T[3][0:],dtype='float')

b_prob = (avg_pmt_counts - dm)/(float(bm - dm))
pmterr_of_mean = pmterr/sqrt(trials)
b_prob_err = pmterr_of_mean/(float(bm - dm))

# Load histgram
data_name = [x for x in files if "_raw.csv" in x][0]
hdata = np.genfromtxt(data_name, delimiter=",", dtype='float')
print(np.shape(hdata))

for i,row in enumerate(hdata):
    det_array = np.copy(row)
    counts_data = parse_raw_counts(det_array)

    z_data = 2*(((counts_data-dm)/(bm - dm)) - 0.5)
    datas.append(z_data)

    bs = np.arange(-1.1,1.1,(2.2/num_bins))

    lab = r"Tip angle (deg)={0:.2f} ms".format(tip_angle_deg[i])
    plt.hist(datas[i],bs,label=lab,alpha=0.6)#, align='right')
    
    k2,pval = mstats.normaltest(datas[i])
    print(lab)
    print("# of trials: {}".format(len(datas[i])))
    print("Normality tests: skew+kurtosis: {0:.4g}, pval: {1:.4g}".format(k2,pval))

os.chdir(base_path)

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


os.chdir(base_path)
if save is True:
    plt.savefig(img_name,bbox='tight',transparent=True)
"""
