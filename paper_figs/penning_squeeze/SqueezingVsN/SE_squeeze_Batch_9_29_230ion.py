# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 14:42:25 2015

@author: jgb
"""

import os, importlib
import numpy as np
from numpy import pi, sqrt, sin
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FormatStrFormatter
import scipy.interpolate as inter

import hfGUIdata as hf
import plot_style as ps
importlib.reload(ps)
import squeeze_func_time as squ
import resample_tools as re
from matplotlib.ticker import FixedLocator, FormatStrFormatter

#options
colors = ['k', ps.red, ps.blue, ps.orange, ps.pink, ps.aqua, ps.navy]
verbose = True
save = True
show_dynamics = False
save_dynamics = False
img_name = "spinNoise_9_29_N230"
legend = True
raw = False

files_to_use = [2,4,3,5]
files_to_use = [2,4,3]
J1k = 1785.0    
Ncal = 0.85
ODF_seq = 2 # use 2 for a simple sequence, 4 for Walsh, set up to get the detunings correct
axis_set = [4,190,-16,17]

#theory calc info
G_el =  62.3
G_ud =  9.35
G_du =  6.59

#adjust for extra decohrence
G_add = 50.0
G_tot = 0.5*(G_el + (G_ud+G_du) + G_add)
print(G_tot)
G_el = G_el + G_add

#added noise from Jy noise fit
A =  0.002143  # rad^2/ms^2
B = 0#.0001607 # rad^2/ms^4


# containers for data sets
psis=[]
its=[]
sig_obs = []
sig_ob_errs = []
sig_robs = []
sig_rob_errs = []
sig_ins = []
sig_2_in_err = []
sig_pns = []
sig_psns = []
SE = []
SN = []
SN_max = []
xi_R2 = []
Ns = []
R_cals = []
names = []

base_path = os.getcwd()
base_path = os.path.normpath("/Users/jgb/Data/20150929/Load330/squeeze")
fns = [os.listdir(base_path)[i] for i in files_to_use]
fns = [i for i in fns if not i.startswith('.DS_Store')]
J1ks = J1k*np.ones(np.shape(fns))
Ncals = Ncal * np.ones(np.shape(fns))  # #photons per ion per ms

#%%
#_____________________________________________________________________
# data processing here
for i,fn in enumerate(fns):
    folder = os.path.join(base_path,fn)
    os.chdir(folder)
    files = os.listdir(os.getcwd())

    #Load properties data
    prop_name = [x for x in files if "_props.csv" in x][0]
    file_name, data_p = hf.get_gen_csv(prop_name, skip_header=False)
    bm = data_p['det_brightMean']
    dm = data_p["det_darkMean"]
    det_t = data_p["det_t"]
    int_t = 2e-6*data_p["squeeze_arm_t"]  #total interaction time in secs
    #reps = data_p["reps"]
    k = bm-dm  # phtns per N atoms
    N = k/(det_t*1e-3)/Ncals[i]

    #load calibration data
    cal_name = [x for x in files if "_cal.csv" in x][0]
    file_name, data = hf.get_gen_csv(cal_name, skip_header=True)
    cal_counts_avg = np.mean(data.T[1])
    sig_cal_ob = np.mean(data.T[2])
    sig_cal_sn = sqrt(cal_counts_avg)
    sig_cal_pn = sqrt(k**2/4.0/N)
    R_cal = (sig_cal_ob**2 - sig_cal_sn**2) / (sig_cal_pn**2)
    
    # load experiment data
    data_name = [x for x in files if "_data.csv" in x][0]
    file_name, data = hf.get_gen_csv(data_name, skip_header=True)
    reps = np.mean(data.T[3])
    psi_deg = data.T[0]
    count_avg = data.T[1]
    sig_sn = sqrt(data.T[1])
    sig_sn_2_err = sig_sn**2 * sqrt(2/reps)
    sig_ob = data.T[2]
    sig_in = sqrt(sig_ob**2 - sig_sn**2)  # subtract poissonian shot noise
    sig_ob_2_err = sig_ob**2 * sqrt(2/reps)
    sig_in_2_err= sqrt( (sig_ob_2_err)**2  + (sig_sn_2_err)**2)

    sig_pn = sqrt(k**2/4.0/N)  # calclated from the atom number

    sig_a = sqrt(sig_cal_ob**2 - sig_cal_sn**2 - sig_pn**2)
    sig_a_deg = (sig_a/cal_counts_avg)*180.0/pi

    dB_squ_in = 10*np.log10((sig_in**2)/(sig_pn**2))
    dB_squ_ob = 10*np.log10((sig_ob**2)/(sig_pn**2))
    
    # calc reduction in constrast (from model shown to represent data)
    Jbar = J1ks[i] /(0.002/int_t)
    out = squ.OAT_decoh(0.0, int_t, Jbar, N, G_el, G_ud, G_du)
    C_coherent_pred = np.real(out[1])
    csi_R2 = (sig_ob**2)/(sig_pn**2)/C_coherent_pred**2
    csi_R2_err = sqrt( (sig_ob_2_err/sig_pn**2/C_coherent_pred**2)**2 + ((sig_ob/sig_pn**2/C_coherent_pred**2)**2 * 0.05)**2 + (csi_R2*0.05)**2 ) # accounting for 5% uncertainty in PN
    dB_csi_R2 = 10*np.log10(csi_R2)
    
    csi_R2_inf = (sig_in**2)/(sig_pn**2)/C_coherent_pred**2
    csi_R2_inf_err = sqrt( (sig_in_2_err/sig_pn**2/C_coherent_pred**2)**2 + ((sig_in/sig_pn**2/C_coherent_pred**2)**2 * 0.05)**2 + (csi_R2_inf*0.05)**2 )
    
    if raw is True:
    # Load histgram
        files = os.listdir(os.getcwd())
        data_name = [x for x in files if "_raw.csv" in x][0]
        hdata = np.genfromtxt(data_name, delimiter=",", dtype='float')
        print(np.shape(hdata))
        
        jack_sig = np.zeros(np.size(sig_ob))
        jack_err_on_sig = np.zeros(np.size(sig_ob))
        for i,row in enumerate(hdata):
            counts = hf.parse_raw_counts(row)
            jack_sig[i], jack_err_on_sig[i] = re.jackknife_est(counts,np.std)
            
        sig_robs.append(jack_sig)
        sig_rob_errs.append(jack_err_on_sig)


    #Load data messages
    print( "______________Data set: {}__________________".format(hf.n_slice(file_name)))
    print( "Ion number from loading calibration: {:.0f}".format(N))
    if verbose:
        print( "Meas. noise sig_cal_ob [Stdv PMT counts]: {:.3f}".format(sig_cal_ob))
        print( "Infer. proj. noise sqrt(cal_ob^2 - cal_sn^2):  {:.3f}".format(sqrt(sig_cal_ob**2 - sig_cal_sn**2)) )
        print( "Pred. proj. noise sig_pn = sqrt(k**2/4.0/N): {:.3f}".format(sig_pn) )
        print( "Est. added noise sqrt((cal_ob^2 - cal_sn^2)-pn^2): {:.3f}, or {:.3f} deg".format(sig_a,sig_a_deg))
        print( "Est. added noise (ratio of var to photon shot noise: {:.3f}".format(sig_a**2/np.mean(sig_sn)**2))
    print( "Minimum inferred spin variance: {:.3f} dB".format(np.min(dB_squ_in))) 
    print( "Minimum observed spin variance: {:.3f} dB".format(np.min(dB_squ_ob)))
    print( r"Minimum observed $\csi_R^2$: {:.3f} dB".format(np.min(dB_csi_R2)))
    
    print( r"Minimum observed $\csi_R^2$: {:.5g} +- {:3g}".format(np.min(csi_R2),csi_R2_err[np.argmin(csi_R2)] ))
    print( r"Minimum inferred $\csi_R^2$: {:.5g} +- {:3g}".format(np.min(csi_R2_inf),csi_R2_inf_err[np.argmin(csi_R2_inf)] ))


    psis.append(psi_deg)
    its.append(int_t)
    sig_obs.append(sig_ob)
    sig_ob_errs.append(sqrt(sig_ob_2_err))
    sig_ins.append(sig_in)
    sig_2_in_err.append(sig_in_2_err)
    sig_pns.append(sig_pn)
    sig_psns.append(sig_sn)
    SE.append(np.max(10*np.log10(csi_R2**(-1))))
    if int_t == 0.0026:
        temp = (sig_ins[i]**2)/(sig_pns[i]**2)
        SN.append(temp[np.argmin(temp)])
    else:
        SN.append(np.min((sig_ins[i]**2)/(sig_pns[i]**2)))
    SN_max.append(np.max((sig_ins[i]**2)/(sig_pns[i]**2))) 
    xi_R2.append(np.min((csi_R2)))
    Ns.append(N)
    R_cals.append(R_cal)
    names.append(hf.n_slice(file_name))
    os.chdir(base_path)

#%%
#________________________________________________________________________
# visualizing the experimental data
fig, ax = plt.subplots() 
SN_err = []
SNmax_err = []
for i,data in enumerate(sig_obs):
    l = r"$\tau=$ {:.3g} ms, N: {:.0f}".format(its[i]*1e3,Ns[i])
    l = r"$\tau=$ {:.3g} ms".format(its[i]*1e3)
    spin_noise = (sig_ins[i]**2)/(sig_pns[i]**2)
    spin_noise_err = sig_2_in_err[i]/(sig_pns[i]**2)
    spin_noise_err = sqrt( (sig_2_in_err[i]/sig_pns[i]**2)**2 + ((sig_ins[i]/sig_pns[i])**2 * 0.05)**2 ) # accounting for 5% uncertainty in PN
    SN_err.append(np.min(spin_noise_err))
    SNmax_err.append(np.max(spin_noise_err))        
    spin_noise_dB = 10*np.log10(spin_noise)
    #spin_noise_err_dB = 10*np.log10(spin_noise) - 10*np.log10(spin_noise-2*spin_noise/sqrt(2*reps))
    spin_noise_err_dB = 10*np.log10(spin_noise) - 10*np.log10(spin_noise-spin_noise_err)
    plt.errorbar(np.abs(psis[i]-180),spin_noise_dB,yerr=spin_noise_err_dB, fmt='o',label=l, color=colors[i])

#plt.yscale('log')
plt.xscale('log')
plt.axis(axis_set)
plt.xlabel(r"Tomography angle $\psi$ (deg)",fontsize=14)
plt.ylabel("Spin variance $(\Delta S_\psi')^2$/N/4 (dB)",fontsize=14)
plt.grid('off')
plt.locator_params(axis='y',nbins=8)

majorLocator = FixedLocator([1,10,100,180])
majorFormatter = FormatStrFormatter('%d')
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)

#________________________________________________________________________
#add some theory curves

psi = np.linspace(0.001,pi,num=500) # radians
for i,name in enumerate(names):
    Jbar = J1ks[i]/(0.002/its[i])
    out = squ.OAT_decoh(-psi, its[i], Jbar, Ns[i], G_el, G_ud, G_du)
    out_perf = squ.OAT_decoh(-psi, its[i], Jbar, Ns[i], 0.0, 0.0, 0.0)
    out_l = squ.OAT_decoh(-psi, its[i], Jbar, Ns[i]-5, G_el, G_ud, G_du)
    C = np.exp(-G_tot*its[i])  # reduction in contrast from spontaneous emission
    R = np.real(out[0]/(sqrt(Ns[i])/(2.0)))**2
    R_perf = np.real(out_perf[0]/(sqrt(Ns[i])/(2.0)))**2
    R_add = R + (A*(its[i]*1e3)**2)*(N*C**2) * sin(psi)**2 + (B*(its[i]*1e3)**4)*(N*C**2) * sin(psi)**2
    R_dB = 10*np.log10(R)
    R_perf_dB = 10*np.log10(R_perf)
    R_add_dB = 10*np.log10(R_add)
    plt.plot(psi*180/pi,R_dB,color=colors[i])
    plt.plot(psi*180/pi,R_perf_dB,'--',color=colors[i])
    plt.plot(psi*180/pi,np.zeros_like(psi),color='gray', zorder=1)
    #plt.plot(psi*180/pi,R_add_dB,color=colors[i],linestyle='-.')
    
    #where is the limit just due to technical noise?
    tech_limit = 10*np.log10((0.3*sig_psns[i]**2)/(sig_pns[i]**2))
    #plt.plot(psis[i],tech_limit,color=colors[i],linestyle='--')
    
    #plt.fill_between(ti*1e3,C_l,C_u,facecolor=colors[j],alpha=0.5)
    print("added dephasing: {:.3g} (ratio of var to proj noise)".format((A*(its[i]*1e3)**2)*N))

if legend is True: plt.legend(loc=0,fontsize=10)

if save is True:
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    plt.savefig(img_name+".pdf",dpi=300,bbox='tight',transparent=True)
    # make a copy of the analysis at the folder
    #shutil.copy(__file__, os.getcwd())
    os.chdir(base_path)

plt.show()
plt.close()

#________________________________________________________________________
# This section shows quantities vs interaction time on y axis
if show_dynamics is True:

    int_times = np.array(its)*1e3
   
    fig, ax = plt.subplots()   
    int_times= np.insert(int_times,0,[0.])
    Ns = np.insert(Ns, 0, [85.0])
    SN = np.insert(SN,0,[np.mean(R_cals[:])])
    SN_max = np.insert(SN_max,0,[np.mean(R_cals[:])])
    SN_err = np.insert(SN_err,0,[0.])
    SNmax_err = np.insert(SNmax_err,0,[0.])
    SN_err_dB = 10*np.log10(np.array(SN) + np.array(SN_err)) - 10*np.log10(SN)
    SNmax_err_dB = 10*np.log10(np.array(SN_max) + np.array(SNmax_err)) - 10*np.log10(SN_max)
    #plot data points
    plt.errorbar(int_times[:],10*np.log10(SN[:]),yerr=SN_err_dB,fmt='o',color=ps.aqua, label=r"Min[$\psi$]")
    plt.errorbar(int_times[:],10*np.log10(SN_max[:]),yerr=SNmax_err_dB,fmt='ko',label="Max[$\psi$]")
    
    #make theory predictions
    it_theory = 1e-3*np.linspace(np.min(int_times),np.max(int_times),num=100)
    SN_the = np.zeros_like(it_theory)
    SN_the_p = np.zeros_like(it_theory)
    SN_max_the = np.zeros_like(it_theory)
    
    N_func = inter.interp1d(1e-3*int_times, Ns,kind='linear')

    for i,t in enumerate(it_theory):
            Jbar = J1k/(0.002/t)
            plusa = 5*pi/180.0
            out_p = squ.OAT_decoh(0.0, t, Jbar, N_func(t), G_el, G_ud, G_du)
            out = squ.OAT_decoh(np.real(out_p[2]), t, Jbar,  N_func(t), G_el, G_ud, G_du)
            out_u = squ.OAT_decoh(np.real(out_p[2])-plusa, t, Jbar,  N_func(t), G_el, G_ud, G_du)
            out_max = squ.OAT_decoh(np.real(out_p[2])+pi/2., t, Jbar,  N_func(t), G_el, G_ud, G_du)
            SN_the[i] = (np.real(out[0])**2) / ( N_func(t)/4.0)
            SN_the_p[i] = (np.real(out_u[0])**2) / ( N_func(t)/4.0)
            SN_max_the[i] = (np.real(out_max[0])**2) / ( N_func(t)/4.0)

    plt.fill_between(it_theory*1e3, 10*np.log10(SN_the), y2=10*np.log10(SN_the_p), color=ps.aqua, alpha=0.3)
    plt.plot(it_theory*1e3, 10*np.log10(SN_the), '-', color=ps.aqua)
    plt.plot(it_theory*1e3, 10*np.log10(SN_max_the),'k-')
    plt.plot(np.linspace(0.0,3.0),np.zeros_like(np.linspace(0.0,3.0)),color='k', zorder=1)
   
    plt.ylabel("Spin variance $(\Delta S_\psi')^2$/N/4 (dB)")
    plt.xlabel('Interaction time (ms)')
    plt.axis([np.min(int_times),np.max(int_times),axis_set[2],axis_set[3]])
    plt.grid('off')
    
    if save_dynamics is True:
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        plt.savefig(img_name+"SN"+".pdf",dpi=300,bbox='tight',transparent=True)
        os.chdir(base_path)
    plt.show()
    plt.close()
    
    #________________________________________________________________________
    # plot of the Ramsey squeezing parameter dyanmics
    xi_R2 = np.insert(xi_R2,0,[1.])
    plt.plot(int_times[:],10*np.log10(xi_R2[:]),'-o')
    plt.ylabel("Squeezing Parameter $\\xi_R^2$ (dB)")
    plt.xlabel("Interaction time $\\tau$ [ms]")

    plt.fill_between(np.linspace(0,3.0), np.zeros_like(np.linspace(0,3.0)), y2=-10*np.ones_like(np.linspace(0,3.0)),color='k',alpha=0.2 )
    plt.axis([0,2.65,-10,20])

    plt.grid('off')
    if save_dynamics is True:
        os.chdir('..')
        #plt.savefig(img_name+"xi_R2"+".pdf",dpi=300,bbox='tight',transparent=False)
        os.chdir(base_path)

    plt.show()
    plt.close()
