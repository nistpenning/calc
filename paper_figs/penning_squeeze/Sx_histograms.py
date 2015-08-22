# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 12:13:03 2015

@author: jgb
"""

import os, shutil
import numpy as np
from numpy import pi, sin, cos, sqrt
import matplotlib.pyplot as plt

import hfGUIdata as hf
import plot_style as ps
import squeeze_func_time as squ

base = os.getcwd()

os.chdir("/Volumes/688/Public/penning_britton/dailyLabBookFiles/2015/20150811/Load306/depolarization/2015-08-11--19.53.00.339")

# Load histgram
files = os.listdir(os.getcwd())
data_name = [x for x in files if "_raw.csv" in x][0]
hdata = np.genfromtxt(data_name, delimiter=",", dtype='float')
print(np.shape(hdata))

#Load properties data
Ncals = 1.7425
prop_name = [x for x in files if "_props.csv" in x][0]
file_name, data_p = hf.get_gen_csv(prop_name, skip_header=False)
bm = data_p['det_brightMean']
dm = data_p["det_darkMean"]
det_t = data_p["det_t"]
int_t = 2e-6*data_p["squeeze_arm_t"]  #total interaction time in secs
k = bm-dm  # phtns per N atoms
N = k/(det_t*1e-3)/Ncals
print(N)

data_name = [x for x in files if "_data.csv" in x][0]
file_name, data = hf.get_gen_csv(data_name, skip_header=True)
arm_time = np.array(data.T[0][0:],dtype='float')

def parse_raw_counts(array):
    bad = 0
    for x in np.nditer(array, op_flags=['readwrite']):
        if x == -1:
            print('Found bad data point')
            bad += 1
            x[...] = -1
        else:
            x[...] = int(x) & 0x1fff
    if bad > 0:
        print("# of bad points: {}".format(bad))

det_array = np.copy(hdata)
parse_raw_counts(det_array)

Sz_array = 2*((det_array-dm)/k) - 1.0

histim = []
bins = [np.arange(-1.1,1.1,n) for n in [0.15,0.1,0.06]]
#for row in [det_array[i] for i in [17]]:
data_sets = [6,11,19]

taus = []
hists = []
bin_list  = []

for j,row in enumerate([Sz_array[i] for i in data_sets]):
    #h, bin_edges = np.histogram(row, density=False, bins=bins[j])
    #histim.append(h)
    l = r"$\tau$ = {0:.2f}".format(2e-3*arm_time[data_sets[j]])
    plt.hist(row,bins[j],histtype='stepfilled',
             alpha=0.35, normed=False,label=l)
    taus.append(2e-3*arm_time[data_sets[j]])
    hists.append(row)
    bin_list.append(bins[j])
             
#plt.imshow(np.transpose(histim), aspect='auto', origin='lower')

plt.hist()

plt.ylabel("Number of trials")
plt.xlabel(r"Spin projection $2 S_x/N$" )
#plt.title("2015-08-11 -- N: {0:.0f}".format(N))
plt.axis([-1.1,1.1,0,90])
plt.legend(fontsize=12)
plt.grid('off')

os.chdir(os.path.dirname(os.path.realpath(__file__)))
plt.savefig("Sx_hist_8_11.pdf",bbox='tight',transparent=True)

ps.save_data_txt('histdata.txt', hists)

os.chdir(base)
#ps.save_data_txt('histdata.txt', histim)
