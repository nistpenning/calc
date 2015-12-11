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
save = False
fig_name = "SqHelDist_10_28.pdf"
files_to_use = [-2]
Ncal = 3.3
hist_to_use = 'all'
bin_width = 2.
h = 50  #block size for resampling
first_points = 1 # how many of the first data sets are at the same angle
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

tipping_angles = tip_angle_rad
#define scale
data_ref = np.hstack(np.array(datas[:first_points]))
bin_def = np.arange(-N/2.0,N/2.0,bin_width)
w = np.zeros_like(data_ref) + (1/np.size(data_ref)) #for weighting the normalized hist  
#calculate reference histogram

vals_ref, bins, patches = plt.hist(data_ref,bin_def,normed=False,weights=w)
plt.close()


#remove the reference data from the rest
datas = datas[first_points:]

#containers
data_hist = []
data_hist_jacks = []
data_hist_jack_errs = []
n_bin_fills = []
samps = []

#first setup a figure for displaying the histograms
num_of_plots = np.size(hdata,axis=0)
num_rows = int(np.ceil(num_of_plots/2.))
plt.figure(1,figsize=(5.0,8.0))
plt.subplot(num_rows,2,1)
l = r"$\theta$: {:.3g} deg".format(tip_angle_rad[0]*180/pi)
gauss = (1/sqrt(2.*pi)/(sqrt(N)/2.))*np.exp(-((bin_def)**2)/2.0/(sqrt(N)/2.)**2)
plt.hist(data_ref, bin_def,normed=True)
plt.plot(bin_def,gauss,color=ps.red,label=l)
plt.locator_params(axis='y',nbins=2)
plt.title(l,fontsize=9)
    
for i,row in enumerate(datas):
    plt.figure(1,figsize=(5.0,8.0))
    l = r"$\theta$: {:.3g} deg".format(tip_angle_rad[i+1]*180/pi)
    plt.subplot(num_rows,2,i+2)
    gauss = (1/sqrt(2.*pi)/(sqrt(N)/2.))*np.exp(-((bin_def)**2)/2.0/(sqrt(N)/2.)**2)
    plt.hist(row, bin_def,normed=True)
    plt.plot(bin_def,gauss,color=ps.red,label=l)
    plt.locator_params(axis='y',nbins=2)
    plt.title(l,fontsize=9)
    
    plt.figure(2)
    w = (1/np.size(row))+np.zeros_like(row) #for weighting the normalized hist      
    vals = plt.hist(row, bin_def,normed=False,weights=w)[0]
    h_dist = 0.5*np.sum((sqrt(vals_ref) - sqrt(vals))**2)
    samp = np.size(row)
    n_bin_fill = np.size(vals[vals!=0])  #number of discrete bins for which P!=0
    
    #resampling    
    g = int(samp/h)
    jack_trials = np.reshape(row, (g,h))
    re_hds = []
    for i in range(0,g):
        re_row = np.reshape(np.concatenate((jack_trials[:i],jack_trials[i+1:])), (samp-h))
        #vals = plt.hist(re_row, bin_def,normed=True)[0]
        w = np.zeros_like(re_row) + (1/np.size(re_row)) #for weighting the normalized hist  
        re_vals = plt.hist(re_row, bin_def,normed=False,weights=w)[0]
        re_h_dist = 0.5*np.sum((sqrt(vals_ref) - sqrt(re_vals))**2)
        re_hds.append(re_h_dist)
    h_dist_jack = g*h_dist - ((g-1)/float(g)* np.sum(re_hds))
    h_dist_jack_err = np.sqrt(np.var(re_hds))

    
    data_hist.append(h_dist)
    data_hist_jacks.append(h_dist_jack)
    data_hist_jack_errs.append(h_dist_jack_err)
    n_bin_fills.append(n_bin_fill)
    samps.append(samp)

l=['Tipping angle (rad)', r'dh$^2$','Extract Fisher Info']
plt.close()

plt.figure(1)
plt.tight_layout()
plt.show()
plt.close()

#note: tipping angle starts at 1 because 0 corresponds to the reference
if hist_to_use == 'all':
    fit_res = pt.plot_polyfit(tipping_angles[first_points:],data_hist_jacks[0:],
                          np.array([0,0,1]),
                          yerr = data_hist_jack_errs[0:],
                          hold=np.array([False,True,False]),
                          labels=l)

else:
    if hist_to_use>len(tipping_angles[1:]):
        print("Error, hist_to_use must be <= {}".format(len(tipping_angles[1:])))
    else: 
        fit_res = pt.plot_polyfit(tipping_angles[first_points:hist_to_use],data_hist_jacks[0:hist_to_use-1],
                          np.array([0,0,1]),
                          yerr = data_hist_jack_errs[0:hist_to_use-1],
                          hold=np.array([False,True,False]),
                          labels=l)


k2 = fit_res[0][1]
frac_err_k2 = fit_res[1][1]/k2
NormF = 8*k2/N[0] #dont' need the 1/M term from teh estimate
print("F = 8*k2, and then F/N is {:.3g} +- {:.3g}".format(NormF, NormF*frac_err_k2))
print("number of samples is {:.3g}".format(np.mean(samps)))

plt.show()
plt.close(0)

#%%

fig, ax = plt.subplots(1,figsize=[4.0,4.0])
xmax = np.max(tipping_angles)*1.1
ymax = np.max(data_hist_jacks)*1.1
plt.locator_params(axis='y',nbins=5)
plt.locator_params(axis='x',nbins=5)
nonCx = np.linspace(0,xmax,num=200)
nonCy = nonCx**2 * N[0]/8.0
fitCy = fit_res[0][0] + fit_res[0][1]*nonCx**2
#ax.plot(nonCx,nonCy,'b-',fillstyle='top')
ax.fill_between(nonCx, 10*np.ones_like(nonCx), nonCy, facecolor='grey', alpha=0.5)
if hist_to_use == 'all':
    ax.errorbar(tipping_angles[first_points:],data_hist_jacks[0:],
             yerr=data_hist_jack_errs[0:],fmt='ko')
else:
    ax.errorbar(tipping_angles[first_points:hist_to_use],data_hist_jacks[0:hist_to_use-1],
             yerr=data_hist_jack_errs[0:hist_to_use-1],fmt='ko')
ax.plot(nonCx,fitCy,'-',color=ps.red)
ax.set_xlabel(r"Angle $\theta$ (rad)")
ax.set_ylabel("Squared Hellinger Distance")
ax.set_ylim([0,ymax])
ax.set_xlim([0,xmax])

ax.grid()
if save is True:
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    plt.tight_layout()
    plt.savefig(fig_name,dpi=300,bbox='tight',transparent=True)
    os.chdir(base_path)
