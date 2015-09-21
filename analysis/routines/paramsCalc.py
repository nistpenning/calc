# -*- coding: utf-8 -*-
"""
A script to take in common parameters for the penning trap and return useful
number for analysis.
"""

import numpy as np
from scicons import pi, c, hbar, m_Be, mI, mu_b, ge, gIp, mu_n, AhfS

import spont_emission as spe
import ODF

#inputs
shells = 5
power_u = 4.85e-3  # W
power_L = 5.12e-3  # W
wpu = 60.0
wpl = 36.0
ACSS_u_pi = 2*pi*20.25e3
ACSS_L_pi = 2*pi*20.25e3
wz = 2*pi*1570.0*1e3  # 1/s

#known constants
wpu_offset = 3
wpl_offset = 4

Lcooling = 626.2252e-9  # wavemeter reading of red lasers
LODF = 626.2670e-9
Lrepumping = 626.1963e-9
#  get the laser wavelength, handling the mismatch between Plevels calc and wavemeter
wL, Delta_ODF = ODF.calc_wL(B=4.4588, L_c=0.5*Lcooling, L_ODF=0.5*LODF)

##################### calculated values ##################
Nion = 1 + 6* np.sum(range(1,shells+1))

I_fit_u = ODF.ACSS_I_cal(ACSS_u_pi, wL=wL)
I_fit_L = ODF.ACSS_I_cal(ACSS_L_pi, wL=wL)
I_fit_avg = (I_fit_u+I_fit_L)/2.

#  calc angles from exp offset
phip_u = (wpu - wpu_offset)*2  # degress, the offset is found from the best linear polarization
phip_L = (wpl - wpl_offset)*2  # degrees

#  spont emission
G_u = spe.spont_emission(I=I_fit_u, omegaL=wL, phip=phip_u*pi/180.)
G_L = spe.spont_emission(I=I_fit_L, omegaL=wL, phip=phip_L*pi/180.)

F0 = ODF.ODF(I=I_fit_avg, thetaR=10.0, omegaL=wL)

def Jbar(delta): return F0**2/(4*hbar*m_Be*wz*(2*pi*delta*1e3))

#print out
print("Number of ions: {}".format(Nion))
print("ODF freqency-- Absolute: {0:.4g}, Detuning from cooling: {1:.4g} GHz".format(wL,Delta_ODF*1e-9/2./pi) )
print("Spont. emission rates-- Upper: {0:.3g} 1/s, Lower: {1:.3g} 1/s, Sum: {2:.3g} 1/s".format(G_u[2],G_L[2],G_u[2]+G_L[2]))
#print("Intensity Calc-- Upper: {:.3f}, Lower: {:.3f}".format(I_u,I_L))
print("Intensity Fit ACSS-- Upper: {:.3f}, Lower: {:.3f}".format(I_fit_u,I_fit_L))
print("F0_avg: {0:.4g} yN".format(F0*1e24) )
print("Jbar pred, [1,2] kHz: [{:.0f},{:.0f} 1/s]".format(Jbar(1),Jbar(2)))
print("All decay rates (sum upper and lower) -- GammaEl, GammaRam, Gamma, phip, Gammaud, Gammadu:{}".format((np.array(G_u)+np.array(G_L))))
ratios = (np.array(G_u)+np.array(G_L))/(0.64*1500.0/1.11)
print("All as ratios -- {}".format(ratios))