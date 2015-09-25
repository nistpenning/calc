# -*- coding: utf-8 -*-
"""
Created on Sat May  2 12:44:14 2015

@author: justinbohnet
"""
import os, importlib
import numpy as np
from numpy import sin, cos, pi, sqrt
import matplotlib.pylab as plt

import plot_model_fit as pt
import spont_emission as spe
import ODF
importlib.reload(ODF)
from scicons import pi, hbar, m_Be, k_b
import hfGUIdata
importlib.reload(hfGUIdata)
       
folders = [hfGUIdata.get_immediate_subdirectories(os.getcwd())[i] for i in [-1]]

#known inputs
# 9 shells
F12 = 2*pi*185.5e3  # per sec
lACSS = 2*pi*20.3e3  # per sec
uACSS = 2*pi*20.3e3  # per sec

Gamma, Jbar_1kHz, F0 = ODF.IsingCalc(lACSS, uACSS, 2*pi*1573.0*1e3)
Gamma = Gamma
wa_hold = False

Nion = 150

print("F0: {}".format(F0))
print('Gamma: {}'.format(Gamma))  

#%%#define analysis function

def analysis():
    #get data
    file_name, scandata, avg_pmt_counts, pmterr, trials, data = hfGUIdata.get_raw_counts()

    b_prob = hfGUIdata.bright_fraction(avg_pmt_counts)
    
    #get values from ion properties
    fz_str = 'raman%raman_fz'
    w_a = hfGUIdata.get_ionProp_value(fz_str)
    w_a = 2*pi*w_a*1e3
    w_a = 2*pi*1571.0*1e3
    
    #Fit guesses
    K0 = (F0/ sqrt(2.0))**2 / (hbar*2*m_Be*Nion*w_a) #sqrt(2) tries to account for DW factor
    n = 10.0 

    arm_time = np.mean(data['arm_t'])*1e-6
    pi_time = np.mean(data['middle_t'])*1e-6
    pre_time = np.mean(data['pre_t'])*1e-3

    mu_R = scandata
    tau = arm_time
    
    #define fit model
    guess = np.array([n, K0, Gamma, w_a])
    hold = np.array([False, True, True, wa_hold])

    #plot labels
    extent = [1564.0, 1580.0, 0.0, 0.60]
    l = ['Raman detuning [kHz]', 'Bright Fraction', '%s, $t_{a}$: %d us, $t_{pre}$: %d ms'%(file_name[-14:-8],arm_time*1e6,pre_time)]

    # print out relevant parameters
    print("file name: {}".format(file_name))
    print("pre_time: {} ms".format(pre_time*1e-3))

    
    def COMfit(mu_R, nm, K0, G, w_a):
    
        mu_R = 2*pi*mu_R*1e3  # convert to 1/s

        res = np.array([])
        for m in mu_R:
            delta = m - w_a
            phi = (tau + pi_time)*delta
            
            asq = 2*sin(phi/2.)**2/(m**2-w_a**2)**2 * ((m**2-w_a**2)*(cos(phi)+cos(m*tau+phi)) - 
                            2*(-m**2 -w_a**2 +
                            cos(m*tau)*cos(w_a*tau)*(m**2*(1+cos(phi))+w_a**2*(1-cos(phi))) +
                            (w_a**2 - m**2)*cos(w_a*tau)*sin(m*tau)*sin(phi) +
                            2*m*w_a*sin(m*tau)*sin(w_a*tau)))
            
            alpha =  w_a*(1-cos(phi)) + 1j*m*sin(phi) - \
                        (np.exp(1j*w_a*tau) * 
                        (w_a*(cos(m*tau)-cos(m*tau+phi)) -  1j*(sin(m*tau)-sin(m*tau+phi))))
            #asq = alpha*alpha.conj()/((m**2-w_a**2)**2)
            res = np.append(res, asq)
        res = np.abs(res)
        return 0.5*(1 - np.exp(-G*2*tau)*np.exp(-2*K0*res*(2*nm+1)))        
    popt, perr = res, res_err = pt.plot_fit(mu_R, b_prob, COMfit, guess,
                           hold=hold, axis=extent, labels=l, save=True)
    
    return popt, perr, det, pre_time, file_name

#%%

#bins for data
data_set = np.zeros(0)
nbar = np.zeros(0)
nbar_err = np.zeros(0)
temp = np.zeros(0)
temp_err = np.zeros(0)
w_a = np.zeros(0)
pre_t = np.zeros(0)

for folder in folders:
    os.chdir(folder)
    
    plt.close()
    res, res_err, det, pre_time, name = analysis()
    plt.show()

    nbar = np.append(nbar, res[0])
    nbar_err = np.append(nbar_err, res_err[0])
    if wa_hold is False:
        w_a = np.append(w_a, res[1])
        t = res[0]*hbar*res[1]/k_b
        t_err = res_err[0] * hbar*res[1]/k_b  # assume most error in n, not w_a
    else:
        w = 2*pi*1487.6*1e3
        w_a = np.append(w_a, w)
        t = res[0]*hbar*w/k_b
        t_err = res_err[0] * hbar*w/k_b  # assume most error in n, not w_a
    temp = np.append(temp,t)
    temp_err = np.append(temp_err, t_err)
    data_set = np.append(data_set,folder)
    pre_t = np.append(pre_t,pre_time)
    
    print("Fitted nbar: {0:.3g} 1/s".format(res[0]))
    print("nbar uncertianty: {0:.3g} 1/s".format(res_err[0]))
    print("Temperature Estimate: {:.4g} mK".format(t*1e3))
    
    os.chdir('..')

#%%
plt.errorbar(pre_t, temp, yerr=temp_err, fmt='o')
plt.ylabel('Temperature [mK]')
plt.xlabel('Heating time')
plt.title('Endcaps grounded')
plt.show()

plt.close()

def lin(x,y0,m): 
    return y0+m*x

l = ['Heating time [ms]', 'nbar', 'Extract heating rate']
pt.plot_fit(pre_t[pre_t<4], nbar[pre_t<4], lin, np.array([5.0,500.0]),
            hold=np.array([False, False]),
            yerr=nbar_err[pre_t<4],
            labels=l)               