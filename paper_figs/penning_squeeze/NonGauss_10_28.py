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
Ncal = 1.1598
verbose = True
save = True
files_to_use = [3]
hist_to_use = [0,1,2,4]
axis_list = [-1.1,1.1,0.0,0.30]
text_name = "batch_hist_1028"
img_name = "batch_hist_img_1028"
num_bins = 27#sqrt(len(z_data))
base_path = os.path.normpath("/Users/jgb/Data/20151028/Load 346/squeeze")
data_path = base_path
os.chdir(data_path)

#load theory data
show_theory = False
if show_theory is True:
    theory_file = "ConvolvedTheoryPDF_N127.txt"
    theory_path = os.path.normpath("/Users/jgb/Data/20151016/")
    os.chdir(theory_path)
    
    tdata = np.genfromtxt(theory_file, dtype='float', names=True)

#%%

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
    print("__________________________________________________________________")
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
    bm = float(data_p['det_brightMean'])
    dm = float(data_p["det_darkMean"])-4.0
    det_t = float(data_p["det_t"])
    int_time = 2*data_p["squeeze_arm_t"]
    k = bm-dm  # phtns per N atoms
    N = k/(det_t*1e-3)/Ncals[i]
    print("Number of ions: {:.4g}".format(N))
    
    #load batch data
    data_name = [x for x in files if "_data.csv" in x][0]
    file_name, data = hf.get_gen_csv(data_name, skip_header=True)
    scan_data = np.array(data.T[0][0:],dtype='float')
    avg_pmt_counts = np.array(data.T[1][0:],dtype='float')
    
    sig_S_psn = np.mean(sqrt(avg_pmt_counts)*N/k)
    print("Std Dev of S_psi from PSN (norm to N/2): {:5g}".format(sig_S_psn / (N/2.)))
    print("Ratio SNvar/PNvar: {:.5g}".format((sig_S_psn/(N/2.) / (1/sqrt(N)))**2))

    hdata_to_use = np.array([hdata[i] for i in hist_to_use])
    
    num_plots = len(hist_to_use)
    fill_c = [ps.purple,ps.yellow]

    for i,row in enumerate(hdata_to_use):
        plt.figure(figsize=(6.0, 3.0))        
        
        print("_________________________________________")
        det_array = np.copy(row)
        counts_data = hf.parse_raw_counts(det_array)
    
        Sz_data = 2*(((counts_data-dm)/(bm - dm)) - 0.5)
        datas.append(Sz_data)
    
        bs = np.arange(-1.0,1.01,(2.0/num_bins))
        trials = len(datas[i])

        lab = r"Scan data = {0:.4g}, Squeeze time = {1:.4g}".format(scan_data[i], int_time)
        lab = r"Scan value: {0:.4g}, Squeeze time: {1:.4g}, Detect time: {2:.4g} ms".format(scan_data[hist_to_use[i]], int_time, det_t*1e-3)
        
        
        w = np.zeros_like(datas[i]) + (1/trials) #for weighting the normalized hist
        vals, bins, patchs = plt.hist(datas[i],bs,label=lab,alpha=1,normed=False,weights=w,
                 histtype='stepfilled',color='lightgray',edgecolor='k')  # rel freq
        """
        vals,bins,patchs = plt.hist(datas[i],bs,label=lab,alpha=1,normed=True,
                 histtype='stepfilled',color='lightgray',edgecolor='k')  #PDF
        """
        #gauss = (2.0/num_bins)*trials*(sqrt(N)/sqrt(2*pi))*np.exp(-((bs*sqrt(N))**2)/2.0) #frequency
        gauss = (2.0/num_bins)*(sqrt(N)/sqrt(2*pi))*np.exp(-(((bs-np.mean(datas[i]))*sqrt(N))**2)/2.0) #rel. freq
        #gauss = (sqrt(N)/sqrt(2*pi))*np.exp(-(((bs-np.mean(datas[i]))*sqrt(N))**2)/2.0)  # PDF
        plt.plot(bs,gauss,color=ps.red)
        plt.xlabel(r"Spin projection 2$S_\psi$/N")
        plt.ylabel("Probability")
        plt.grid('off')
        if save is True:
            plt.axis(axis_list)
            plt.locator_params(axis='y',nbins=3)
            plt.tight_layout()
            os.chdir(os.path.dirname(os.path.realpath(__file__)))
            plt.savefig(text_name+str(scan_data[hist_to_use[i]])+"_"+str(int_time)+".pdf",dpi=300,bbox='tight',transparent=True)
            os.chdir(base_path)
        
        """        
        #add theory curves
        if hist_to_use[i] == 0:
            yt = np.zeros_like(tdata['pdf88'])
            plt.plot(bs,gauss,color=ps.red)
            axis_list = [-1.1,1.1,0.0,5.5]
            plt.axis(axis_list)
            plt.locator_params(axis='y',nbins=3)
            os.chdir(os.path.dirname(os.path.realpath(__file__)))
            plt.tight_layout()
            #plt.savefig('127ions_pdf180'+".pdf",dpi=300,bbox='tight',transparent=True)
            os.chdir(base_path)
        elif hist_to_use[i] == 5:
            yt = tdata['pdf88']
            plt.plot(tdata['S_psi'],yt,'-',color='k')
            ymax = 3.0
            axis_list = [-1.1,1.1,0.0,ymax]
            plt.axis(axis_list)
            plt.locator_params(axis='y',nbins=3)
            os.chdir(os.path.dirname(os.path.realpath(__file__)))
            plt.tight_layout()
            #plt.savefig('127ions_pdf88'+".pdf",dpi=300,bbox='tight',transparent=True)
            os.chdir(base_path)
        elif hist_to_use[i] == 2:
            bs_to_use = bins[1:]
            yt = tdata['pdf174p6']
            plt.plot(bs,gauss,color=ps.red)
            plt.plot(tdata['S_psi'],yt,'-',color='k')
            scaling = (yt/((sqrt(N)/sqrt(2*pi))*np.exp(-(((tdata['S_psi']*sqrt(N))**2)/2.0))))
            yt_save = yt
            scaling = (vals/((sqrt(N)/sqrt(2*pi))*np.exp(-(((bs_to_use-np.mean(datas[i]))*sqrt(N))**2)/2.0)))            
            yt_save = vals
            ymax = 8.5
            axis_list = [-1.1,1.1,0.0,ymax]        
            plt.axis(axis_list)
            plt.locator_params(axis='y',nbins=4)
            os.chdir(os.path.dirname(os.path.realpath(__file__)))
            plt.tight_layout()
            #plt.savefig('127ions_pdf174p6'+".pdf",dpi=300,bbox='tight',transparent=True)
            os.chdir(base_path)
        """
        
        plt.show()
        plt.close()
        
        k2,pval = mstats.normaltest(datas[i])
        if type(k2) == np.ma.core.MaskedConstant:
            k2 = 0.0
            pval = 0.0
            print("Normality test had divide by zero error")
        print(lab)
        print("# of trials: {}".format(trials))
        print("Mean: {0:.3g}, Median: {1:.3g}, Min: {2:.3g}, Max {3:.3g}".format(np.mean(datas[i]),np.median(datas[i]),np.min(datas[i]),np.max(datas[i])))
        print("Std Dev: {0:3g}, Variance: {1:.3g}".format(np.std(datas[i]),np.var(datas[i])))        
        print("Normality tests: skew+kurtosis: {0:.4g}, pval: {1:.4g}".format(k2,pval))

    os.chdir(data_path)
    
if save is True:
    os.chdir(os.path.dirname(os.path.realpath(__file__)))    
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