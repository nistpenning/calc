# -*- coding: utf-8 -*-

import os
import numpy as np
from numpy import pi, sin, cos, sqrt
import scipy.optimize
import matplotlib.pyplot as plt

import hfGUIdata as hf
import plot_tools_jgb as pt
import squeeze_func_time as squ

props = [hf.brightMean, hf.darkMean, hf.det_t]

save = False
name = "phaseflop_datasets_8_06.txt"

# containers for data sets
ats=[]
C_est = []
C_est_err = []
Cs = []
Cerrs = []
Ns=[]
phis = []
names = []

base_path = os.getcwd()
fns = os.listdir(base_path)[:-1]
J1ks = (475.0*3.03)*np.ones(np.shape(fns))
Ncals = 1.3999 * np.ones(np.shape(fns))  # #photons per ion per ms

for i,fn in enumerate(fns):
    folder = os.path.join(base_path,fn)
    os.chdir(folder)

    files = os.listdir(os.getcwd())
    ip = hf.get_ionProp_dict(props)
    bm = ip[hf.brightMean]
    dm = ip[hf.darkMean]
    k = bm-dm  # phtns per N atoms
    det_t = ip[hf.det_t]
    N = k/(det_t*1e-3)/Ncals[i]

    file_name, scandata, avg_pmt_counts, pmterr, trials, data = hf.get_raw_counts()
    phi = scandata
    b_prob = (avg_pmt_counts - dm)/(float(bm - dm))
    pmterr_of_mean = pmterr/sqrt(trials)
    b_prob_err = pmterr_of_mean/(float(bm - dm))
    arm_time = np.mean(data['arm_t'])  # us

    def confit(x,C,phi):
        return 0.5 - 0.5* C*(cos(x+phi))

    g = np.array([0.5, 0.0])
    p,pcov = scipy.optimize.curve_fit(confit,phi,b_prob,p0=g,sigma=b_prob_err,
                                           absolute_sigma=True)

    plt.close()
    plt.plot(phi,b_prob,'o')
    plt.plot(phi,confit(phi,p[0],p[1]))

    contrast = p[0]
    contrast_err = np.sqrt(np.diag(pcov))[0]
    contrast_est = 1-(2*b_prob[0])
    contrast_est_err = 2 * b_prob_err[0]

    C_est.append(contrast_est)
    C_est_err.append(contrast_est_err)
    ats.append(arm_time)
    Cs.append(contrast)
    Cerrs.append(contrast_err)
    Ns.append(N)
    names.append(hf.n_slice(fn))

    os.chdir(base_path)

plt.show()

plt.close()

#%%

plt.close()
tau = 2*np.array(ats)/1e3
#plt.errorbar(tau[8:],Cs[8:],yerr=Cerrs[8:],fmt='o')
plt.errorbar(tau[8:],C_est[8:],yerr=C_est_err[8:],fmt='o',c=pt.red)
plt.xlabel("Interaction time[us]")
plt.ylabel("Collective Spin Sx/(N/2)")

N = np.mean(Ns)

#________________________________________________________________________
#add some theory curves
mult = 2
G_el = mult* 60.56
G_ud = mult* 9.075
G_du = mult* 6.413
G_tot = 70.0  # per s

ti = np.linspace(1e-6,4.0e-3,num=100)  # seconds
spem = np.exp(-G_tot*ti)
plt.plot(ti*1e3, spem,'--k',label='Spon. Emiss.')

colors = ['k', pt.red, pt.blue]
for j,Jbar1k in enumerate(J1ks[0:3]):
    Jbar = Jbar1k/(0.002/ti)
    out = squ.OAT_decoh(0.0, ti, Jbar, Ns[j], G_el, G_ud, G_du)
    C_coherent_pred = np.real(out[1])
    #plt.plot(ti*1e3,C_coherent_pred,c=colors[j])
    plt.plot(ti*1e3,C_coherent_pred,c=colors[j])


plt.axis([0,4.1,0.0,1.1])
plt.grid('off')

if save is True:
    plt.savefig(name,dpi=300,bbox='tight',transparent=True)