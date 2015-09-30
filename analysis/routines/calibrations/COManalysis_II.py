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
       
folders = [hfGUIdata.get_immediate_subdirectories(os.getcwd())[i] for i in [0]]

#known inputs
lACSS = 2*pi*18.8e3  # per sec
uACSS = 2*pi*18.8e3  # per sec
w_z = 2*pi*1573.0*1e3
wa_hold = False

Gamma, Jbar_1kHz, F0 = ODF.IsingCalc(lACSS, uACSS, w_z)

Nion = 145.0
Jbar_1kHz = 1753.0  # per s, measured for this data set
F0 = sqrt(Jbar_1kHz * (4*hbar*m_Be*w_z*(2*pi*1.0*1e3)))

print("F0: {}".format(F0))
print('Gamma: {}'.format(Gamma))  

#%%
#fit func def
def COMfit(scandata, nm, K0, G, w_a, tau):
    """
    scandata: ODF beatnote, kHz
    nm: motional quanta
    K0:
    G: total Ramsey contrast decay rate, 1/s
    w_a: axial mode frequency, 1/s
    """
    mu_R = 2*pi*scandata*1e3
    res = np.array([])
    for m in mu_R:
        delta = m - w_a
        phi = (tau + pi_time)*delta
        
        asq = 2*sin(phi/2.)**2/(m**2-w_a**2)**2 * ((m**2-w_a**2)*(cos(phi)+cos(m*tau+phi)) - 
                        2*(-m**2 -w_a**2 +
                        cos(m*tau)*cos(w_a*tau)*(m**2*(1+cos(phi))+w_a**2*(1-cos(phi))) +
                        (w_a**2 - m**2)*cos(w_a*tau)*sin(m*tau)*sin(phi) +
                        2*m*w_a*sin(m*tau)*sin(w_a*tau)))

        res = np.append(res, asq)
    res = np.abs(res)
    return 0.5*(1 - np.exp(-G*2*tau)*np.exp(-2*K0*res*(2*nm+1)))        

#bins for data
data_set = np.zeros(0)
G_fit = np.zeros(0)
G_fit_err = np.zeros(0)
w_rs = np.zeros(0)
w_as = np.zeros(0)

for folder in folders:
    os.chdir(folder)
    #get data
    file_name, scandata, avg_pmt_counts, pmterr, trials, data = hfGUIdata.get_raw_counts()

    b_prob = hfGUIdata.bright_fraction(avg_pmt_counts)
    wr = hfGUIdata.get_ionProp_value('wall%wall_f')
    t = np.mean(data['arm_t'])*1e-6
    pi_time = np.mean(data['middle_t'])*1e-6

    #Fit guesses
    K0 = (F0)**2 / (hbar*2*m_Be*Nion*w_z) 
    n = 7.0
    w_a_guess = np.mean(scandata)
    print("__________________________________________")
    print("Axial Freq guess: {0:.6g}".format(w_a_guess))
   
    #define fit model
    guess = np.array([n, K0, 1.6*Gamma, 2*pi*w_a_guess*1e3,t])
    hold = np.array([False, True, False, wa_hold,True])

    #plot labels
    extent = [np.min(scandata), np.max(scandata), 0.0, 0.60]
    l = [r'Raman detuning [kHz]', 'Bright Fraction', '%s, $t_{a}$: %d us, $\omega_{r}$: %d kHz'%(file_name[-14:-8],t*1e6,wr)]
    
    plt.close()
    res, res_err = pt.plot_fit(scandata, b_prob, COMfit, guess,
                           hold=hold, axis=extent, labels=l, save=False)
    plt.show()
    
    
    G_fit = np.append(G_fit, res[0])
    G_fit_err = np.append(G_fit_err, res_err[0])
    w_rs = np.append(w_rs, wr)
    w_as = np.append(w_as, res[1])
    data_set = np.append(data_set,folder)
    
    print("Fitted Gamma: {0:.3g} 1/s".format(res[0]))
    print("Fitted Gamma uncert: {0:.3g} 1/s".format(res_err[0]))
    
    os.chdir('..')
