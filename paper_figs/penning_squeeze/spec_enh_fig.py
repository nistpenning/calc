# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 11:24:06 2015

@author: jgb
"""

import os,importlib, shutil
import numpy as np
from numpy import sqrt
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FormatStrFormatter

import hfGUIdata as hf
import plot_style as ps
importlib.reload(ps)
import squeeze_func_time as squ

save = True
img_name = "spec_enh"

base_path = os.getcwd()

#first look into the parameter database needed for each data set
db_path = "/Users/jgb/calc/paper_figs/penning_squeeze/sq_data_paramdb.csv"
datadb = np.genfromtxt(db_path, delimiter=",", names=True, dtype=None)

Ns = []
data_names = []
dataSE_names = []
N_SEs =[]
sig_PN_subs = []
sig_PN_fulls = []
sig_PN_fulls_errs = []
sig_SE_fulls = []
sig_SE_fulls_errs = []

#first, just get data sets that were just measuring projection noise of CSS
data_path = "/Users/jgb/Data/20150813/Collect_justPN_data"
folders = os.listdir(data_path)

for i,fn in enumerate(folders[1:]):
    folder = os.path.join(data_path,fn)
    os.chdir(folder)
    files = os.listdir(os.getcwd())

#Load properties data
    prop_name = [x for x in files if "_props.csv" in x][0]
    file_name, data_p = hf.get_gen_csv(prop_name, skip_header=False)
    bm = data_p['det_brightMean']
    dm = data_p["det_darkMean"]
    det_t = data_p["det_t"]
    reps_flag = 0
    try: 
        reps = data_p["reps"]
    except ValueError:
        print(file_name+" no reps in props, look in data")
        reps_flag = 1
    Ncal = data_p['pho_ion_ms']
    k = bm-dm  # phtns per N atoms
    N = k/(det_t*1e-3)/Ncal

    # load experiment data
    data_name = [x for x in files if "_data.csv" in x][0]
    file_name, data = hf.get_gen_csv(data_name, skip_header=True)
    if reps_flag == 1:
        reps = np.mean(data.T[3])
    avg_pmt_counts = np.array(data.T[1][0:],dtype='float')
    pmterr = np.array(data.T[2][0:],dtype='float')
    pmterr_err = pmterr/np.sqrt(2*reps)
    m_len = avg_pmt_counts-dm
    
    sig_sub = np.mean(sqrt(pmterr**2 - avg_pmt_counts)/m_len)
    sig_full = np.mean(pmterr/m_len)
    sig_full_err = np.std((pmterr/m_len))/sqrt(2*np.size(avg_pmt_counts))
    
    Ns.append(N)  
    sig_PN_subs.append(sig_sub)
    sig_PN_fulls.append(sig_full)  
    sig_PN_fulls_errs.append(sig_full_err)
    data_names.append(file_name)
    os.chdir(base_path)

# now get CSS projection noise from calibration data
data_path = "/Users/jgb/Data/20150813/Collect_PNvsN_data"

folders = os.listdir(data_path)

for i,fn in enumerate(folders[1:]):
    folder = os.path.join(data_path,fn)
    os.chdir(folder)
    files = os.listdir(os.getcwd())

#Load properties data
    prop_name = [x for x in files if "_props.csv" in x][0]
    file_name, data_p = hf.get_gen_csv(prop_name, skip_header=False)
    bm = data_p['det_brightMean']
    dm = data_p["det_darkMean"]
    det_t = data_p["det_t"]
    reps_flag = 0
    try: 
        reps = data_p["reps"]
    except ValueError:
        print(file_name+" no reps in props, look in data")
        reps_flag = 1
    #get data base entry
    date = file_name[:10]
    match = np.array([i.decode('ascii')==date for i in datadb['data_set']])
    row = datadb[match]
    Ncal = float(row['Ncal'])
    k = bm-dm  # phtns per N atoms
    N = k/(det_t*1e-3)/Ncal
    
#load calibration data
    cal_name = [x for x in files if "_cal.csv" in x][0]
    file_name, data = hf.get_gen_csv(cal_name, skip_header=True)
    if reps_flag == 1:
        reps = np.mean(data.T[3])
    cal_counts_avg = np.mean(data.T[1])
    m_len = (cal_counts_avg-dm)
    sig_cal_ob = np.mean(data.T[2])/m_len
    sig_cal_ob_err = sqrt(sig_cal_ob**2/sqrt(2*reps))
    sig_cal_pn = sqrt((np.mean(data.T[2])**2 - cal_counts_avg))/m_len
    
    Ns.append(N)  
    sig_PN_subs.append(sig_cal_pn)
    sig_PN_fulls.append(sig_cal_ob)  
    sig_PN_fulls_errs.append(sig_cal_ob_err)
    data_names.append(file_name)
    os.chdir(base_path)

fig, ax = plt.subplots(figsize=(5.0,3.7)) 
Nround = np.array([round(n) for n in Ns])
plt.errorbar(Nround, np.array(sig_PN_subs)**2, fmt='o')
plt.errorbar(Nround, np.array(sig_PN_fulls)**2, fmt='d')
N_pred = np.linspace(18,250)
plt.plot(N_pred,1/np.array(N_pred),'-')

#sfe = 2* np.array(sig_fulls) * np.array(sig_fulls_errs)


#_________________________________________________
#Here get actual spectroscopic enhancement

data_path = "/Users/jgb/Data/20150813/Collect_Spec_Enh_Data"

folders = os.listdir(data_path)

for i,fn in enumerate(folders[1:]):
    folder = os.path.join(data_path,fn)
    os.chdir(folder)
    files = os.listdir(os.getcwd())

#Load properties data
    prop_name = [x for x in files if "_props.csv" in x][0]
    file_name, data_p = hf.get_gen_csv(prop_name, skip_header=False)
    try:    
        arm_t = data_p['squeeze_arm_t']
    except ValueError:
        print(file_name+" No arm_t found. Use 500us")
        arm_t = 500.0
    bm = data_p['det_brightMean']
    dm = data_p["det_darkMean"]
    k = bm-dm  # phtns per N atoms
    det_t = data_p["det_t"]
    

    reps_flag = 0
    try: 
        reps = data_p["reps"]
    except ValueError:
        print(file_name+" no reps in props, look in data")
        reps_flag = 1
    #get data base entry
    date = file_name[:10]
    match = np.array([i.decode('ascii')==date for i in datadb['data_set']])
    row = datadb[match]
    Ncal = float(row['Ncal'])
    N = k/(det_t*1e-3)/Ncal
    J1k = row['J1k']
    walsh_number = row['walsh_num']
    G_el = row['Gel'] + row['Gadd']
    G_ud = row['Gud']
    G_du = row['Gdu']
    
    int_t = walsh_number * 1e-6 * arm_t
    Jbar = J1k/(1000.0/arm_t)
    out = squ.OAT_decoh(0.0, int_t, Jbar, N, G_el, G_ud, G_du)
    C_coherent_pred = np.real(out[1])*(k/2.0)
    
# load experiment data
    data_name = [x for x in files if "_data.csv" in x][0]
    file_name, data = hf.get_gen_csv(data_name, skip_header=True)
    if reps_flag == 1:
        reps = np.mean(data.T[3])
    arm_time = np.array(data.T[0][0:],dtype='float')
    avg_pmt_counts = np.array(data.T[1][0:],dtype='float')
    pmterr = np.array(data.T[2][0:],dtype='float')
    pmterr_err = pmterr/np.sqrt(2*reps)
    m_len = avg_pmt_counts-dm
    
    sig_sub = np.min(pmterr**2 - avg_pmt_counts)**2/C_coherent_pred**2
    sig_full = np.min(pmterr)**2/C_coherent_pred**2
    sig_full_err = sig_full/sqrt(2*reps)
    
    dataSE_names.append(file_name)
    N_SEs.append(N)  
    sig_SE_fulls.append(sig_full)  
    sig_SE_fulls_errs.append(sig_full_err)
    os.chdir(base_path)

NSEround = np.array([round(n) for n in N_SEs])
plt.errorbar(NSEround, np.array(sig_SE_fulls), fmt='s',color=ps.blue)
plt.yscale('log')
plt.xscale('log')
plt.axis([18,250, 0.001,0.06])
plt.ylabel(r'Spin variance $(\Delta S_z/\langle|\vec{S}|\rangle)^2$')
plt.xlabel('Ion number N')
plt.grid('off')

majorLocator = FixedLocator([20,50,100,200])
majorFormatter = FormatStrFormatter('%d')
majorLocatorY = FixedLocator([0.001,0.01,0.06])
majorFormatterY = FormatStrFormatter('%g')
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_major_locator(majorLocatorY)
ax.yaxis.set_major_formatter(majorFormatterY)

dephase_est = 0.002143  # for interaction time of 1ms, from 9/29
x = np.linspace(10,300,num=10)
heisenberg = 1/x**2
plt.plot(x,dephase_est*np.ones(np.shape(x)),'r--' )
plt.plot(x,heisenberg,'k')

if save is True:
    os.chdir('..')
    # make a copy of the analysis at the folder
    shutil.copy(__file__, os.path.normpath("/Users/jgb/Data/20150813"))
    os.chdir(base_path)
    #save figure in the dir with the script, since it is a figure maker
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    plt.tight_layout()
    plt.savefig(img_name+".pdf",dpi=300,bbox='tight',transparent=True)
    os.chdir(base_path)


plt.show()

#labels = ['data_set','N', 'N_err','sig_full_rad', 'sig_sub_rad', 'sig_full_err_rad']
#data = [folders,Ns,Ns_err, sig_fulls, sig_subs,sqrt(sfe)]
#ps.save_data_txt("PNvsN_autodata.csv",data,labels)