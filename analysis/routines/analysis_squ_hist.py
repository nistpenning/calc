# -*- coding: utf-8 -*-
"""
Created on Sun May 10 11:46:38 2015

@author: justinbohnet
"""

import os, csv, sys
import numpy as np
from numpy import sin, cos, pi, sqrt
from scicons import hbar, m_Be, k_b
import matplotlib.pyplot as plt

sys.path.append("..")
import squeeze_func_time as squ
reload(squ)

import hfGUIdata
reload(hfGUIdata)
import plot_model_fit as pt
reload(pt)
import scipy
import scipy.stats

#inputs to the analysis
dataname = "1100us"
csvname = dataname + "_summary.csv"
N = 74.0
N_err = 4.0
Gamma = 78.3-6
confine_param = 0.8**2  # measured from coherent rotation
Jbar_1kHz = 2700 * confine_param # 1/s, predicted from F0
#w_a = hfGUIdata.get_ionProp_value('raman%raman_fz')
#calibration mask range
mask_range = 0  # [8,9]
fitguess = np.array([N, 0.01, 0.02,])
hold = np.array([True, False, False])

'''
Use ACSS to predict Ising interaction and decoherence
ACSS_u_pi = 2*pi*14.4e3
ACSS_L_pi = 2*pi*15.1e3
wz = 2*pi*1556.0e3

Gamma, Jbar_1kHz, F0 = ODF.IsingCalc(ACSS_u_pi, ACSS_L_pi, wz)
'''
folders = hfGUIdata.get_immediate_subdirectories(os.getcwd())

'''
# Use data set of just a rotated Bloch vector for noise calibration
os.chdir(folders[1])
max_c, min_c, N, sigA, k0 = squ.cal_analysis(mask_range, hold=hold, Nguess=N)
print('Ion number (fixed): {}'.format(N))
# get data for contrast compare plot
file_name, scandata, m_full, pmterr_full, trials, data = hfGUIdata.get_raw_counts()
brightMean = hfGUIdata.get_ionProp_value('detection%det_brightMean')
darkMean = hfGUIdata.get_ionProp_value('sf%fitParams%sf_fitParam_darkMean')
pmterr_full = pmterr_full/sqrt(trials-1)
phi_full = ((scandata/pi) * 180.0)  # degrees
os.chdir('..')
'''

# Squeezing data set
os.chdir(folders[0])
max_c = hfGUIdata.get_ionProp_value('detection%det_brightMean')
min_c = hfGUIdata.get_ionProp_value('sf%fitParams%sf_fitParam_darkMean')
#int_time, detuning, m, pmterr_min, theta, j_z_err, psi, stddevjz, RO = squ.sq_analysis(max_c, min_c, N, N_err, sigA ,k0, Jbar_1kHz)
#squ.sq_figures(max_c, min_c, N, N_err, sigA, k0, Jbar_1kHz, save=True, extent=[-5,350,-6.0,12.5])
#counts_data, x_data = squ.data_point_histogram(10, max_c, min_c, N, N_err, sigA, k0, Jbar_1kHz)
squ.hist_data_browser(max_c, min_c, 10)
os.chdir('..')

#%%
#______________________________________________________________________
os.chdir(folders[0])
pn=0
file_name, scandata, counts_data, data = hfGUIdata.get_raw_counts_hist()
arm_time = np.mean(data['arm_t'])
pi_time = np.mean(data['middle_t'])
det_n = np.mean(data['det_n'])
rot_time = np.mean(data['final_t']) * np.ones(np.size(scandata))

# Calculate derived quantities
detune = det_n*1e3/(arm_time)  # kHz
phi = (rot_time/pi_time*180.0)  # degrees

z_data = (counts_data-min_c)/(max_c - min_c)
z_data = z_data.flatten()

plt.close()
bin_width = 0.01
bin_range = np.arange(min(z_data), max(z_data) + bin_width, bin_width)
plt.hist(z_data,bins=bin_range)

param = scipy.stats.norm.fit(z_data)
pdf_x = np.linspace(0,1)
pdf_fitted = (np.size(z_data)*bin_width) * scipy.stats.norm.pdf(pdf_x, loc=param[0], scale=param[1])
plt.plot(pdf_x, pdf_fitted)

plt.xlabel('Counts')
plt.ylabel('Number of Instances')
title = 'Angle:{0:.1f} deg, t_a: {1} us, #{2}'.format(phi[pn],arm_time,pn)
plt.title(title)
plt.axis([0,1.0,0,400.0])

#___________________________________________________________________________

os.chdir('..')
