# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 21:02:26 2015

@author: justinbohnet
"""

import os, importlib, shutil
import numpy as np

import matplotlib.pyplot as plt

import hfGUIdata
importlib.reload(hfGUIdata)
import plot_model_fit as pt
importlib.reload(pt)

# make a copy of the analysis at the folder
if True:
    shutil.copy(__file__, os.getcwd())

def get_immediate_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]

folders = get_immediate_subdirectories(os.getcwd())

#Do the analysis on all the folders
#define intial parameters
Gamma = 80e-6
w_a = 1570.2
Jbar_1kHz = 2000
Jbar_guess = 1.0*Jbar_1kHz*1e-6

fitguess = np.array([Jbar_guess, Gamma])
hold = np.array([False,True])

#define analysis functions

def analysis():

    file_name, scandata, avg_pmt_counts, pmterr, trials, data = hfGUIdata.get_raw_counts()

    b_prob, b_prob_err = hfGUIdata.bright_fraction(avg_pmt_counts, err=pmterr)
    b_prob_err = b_prob_err/np.sqrt(trials-1)

#define fit model
    arm_time = np.mean(data['arm_t'])
    pi_time = np.mean(data['middle_t'])
    n = np.mean(data['det_n'])
    det = n*1e3/(arm_time)
    theta = scandata/pi_time*np.pi*180./np.pi #degrees

    def meanfieldfit(theta, J, Gamma):
        theta = theta * np.pi/180.
        return 0.5*(1 - np.exp(-Gamma*2*arm_time)*np.sin(theta)*np.sin(2*J*np.cos(theta)*2*arm_time))



    title = "Data set: %s, $t_{a}$: %d us, $\delta$: %.3g kHz"%(folder[12:-4],arm_time,det)

    label = ['Inital angle [degrees]', 'Spin up fraction', title]

    popt, perr = pt.plot_fit(theta,b_prob,meanfieldfit, fitguess, yerr=b_prob_err,
                hold=hold, labels=label, axis=[0.0,np.max(theta),0,1.0],save=True)

    return popt, perr, det, file_name

#bins for data
data_set = np.zeros(0)
Jbar = np.zeros(0)
Jbar_err = np.zeros(0)
delta = np.zeros(0)

for folder in folders:
    os.chdir(folder)

    plt.close()
    res, res_err, det, name = analysis()
    plt.show()

    Jbar = np.append(Jbar, res[0]*1e6)
    Jbar_err = np.append(Jbar_err, res_err[0]*1e6)
    data_set = np.append(data_set,folder)
    delta = np.append(delta, det)

    print("Fitted Jbar: {0:.3g} 1/s".format(res[0]*1e6))
    print("Ratio Gamma/Jbar: {0:.3g}".format(Gamma/(res[0])))


    os.chdir('..')
'''
wz = 2*pi*1570e3
def fJbar(delta, A): return A*F0**2/(4*hbar*m_Be*wz*(2*pi*delta*1e3))

mf.plot_fit(delta,Jbar, fJbar, np.array([1]), hold=np.array([True]),
            labels=['$\delta$ from COM [kHz]','Jbar [1/s]','Average interaction strength'],
            axis=[0,11,0,2000])

out = np.array([data_set, Jbar, Jbar_err, delta]).transpose()
np.savetxt('results.csv',out,delimiter=',',fmt="%s")
y = Gamma/Jbar*1e6
plot(delta,y,'o')
xlabel('$\delta$ from COM [kHz]')
ylabel('Gamma/Jbar')
'''
#%%
def get_avg(x, k0):
    return k0 * np.ones(np.shape(x))

y =Jbar*delta
l = ["detuning from COM [kHz]", "Jbar*delta", "Axial Freq. {0} kHz, Jbar Pred: {1}".format(w_a,Jbar_1kHz)]
popt, perr = pt.plot_fit(delta, y, get_avg, np.array([1000.0]),
                         yerr=delta*(Jbar_err),
                            labels=l)
print(popt/Jbar_1kHz)