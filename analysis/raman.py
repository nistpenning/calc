# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 10:25:34 2015

@author: bsawyer
"""
from __future__ import division

import numpy as np
import sys
import scipy
import scipy.constants
import matplotlib.pyplot as plt

sys.path.append("c:\\Users\\bsawyer\\Documents\\Python Scripts\\subroutines")
import wigner
import ODF
import SpontEmission
reload(wigner)
reload(ODF)
reload(SpontEmission)

pi = np.pi    
c = scipy.constants.c
mu_b = scipy.constants.physical_constants['Bohr magneton'][0]
mu_n = mu_b/scipy.constants.physical_constants['proton-electron mass ratio'][0]
hbar = scipy.constants.hbar
epsilon0 = scipy.constants.epsilon_0
ge = scipy.constants.physical_constants['electron g factor'][0] #electron g-factor
gIp = -0.784955 #Be nuclear g-factor
AhfS = -625.008837e6 * 2*pi #Hyperfine coefficient 
AhfP = -118.6e6 * 2*pi #Hyperfine coefficient for P state
mI = +3./2. #nuclear spin projection of Be+ used for experiments

def raman_res(I1 = 1e4, #Beam 1 intensity [W/m^2]
          I2 = 1e4, #Beam 2 intensity [W/m^2]
          phi1 = pi/2, #polarization angle of beam 1 [radians]
          phi2 = 0, #polarization angle of beam 2 [radians]
          Delta = 2*pi*-5e9, #Raman detuning from given intermediate level [s^-1]
          int_state_J = 3/2., # intermediate J state referencing Delta
          int_state_mJ = +1/2., # intermediate mJ state referencing Delta
          B = 4.4588): #magnetic field in [T]

    #calculate Paschen-Bach shift of each qubit level
    w_a = -1*(-0.5*ge*mu_b + mI*gIp*mu_n)*B/hbar + AhfS*mI*(-0.5) #down state
    w_b = -1*(0.5*ge*mu_b + mI*gIp*mu_n)*B/hbar + AhfS*mI*0.5 #up state
    
    #compute spin flip frequency from a and b state shifts
    delta_sf = w_b-w_a
    
    #calculate P level shifts for Be+
    Plevels = Plevelcalc(B)
    Plength = int(np.shape(Plevels)[1])
    #find index of intermediate state
    Pindex = int(3*(1.5 - int_state_J) + (1.5-int_state_mJ))
    w_1 = Delta + Plevels[2][Pindex]-w_a
    
    #compute beam 2 frequency assuming two-photon resonance with spin flip
    w_2 = w_1 - delta_sf
    
    #coefficient that multiplies final Rabi rate
    coeff =3*pi*epsilon0*hbar*c**3 * 2*np.sqrt(I1*I2)/(epsilon0*c) / (2*hbar**2)
    
    Ja = 0.5
    Jb = 0.5
    ma = -0.5
    mb = 0.5
    
    #S is the summand
    S = np.zeros(Plength)

    for ii in range(0,Plength):     
        # w_fa is frequency difference between level f and a
        w_fa = Plevels[2][ii] - w_a
        w_fb = Plevels[2][ii] - w_b
        
        Delta_1 = w_1 - w_fa
        Delta_2 = w_2 - w_fa
        
        Jf = Plevels[0][ii]
        mf = Plevels[1][ii]
        A_T = Plevels[3][ii]
        
        S1 = (-1)**(Jb-mb)*(np.cos(phi2)*wigner.W3J(Jb,1,Jf,mb,0,-1*mf)+\
        0.5*np.sin(phi2)*(wigner.W3J(Jb,1,Jf,mb,1,-1*mf)+wigner.W3J(Jb,1,Jf,mb,-1,-1*mf)))/\
        (w_fb**1.5)
        
        S2 = (-1)**(Jf-mf)*(np.cos(phi1)*wigner.W3J(Jf,1,Ja,mf,0,-1*ma)+\
        0.5*np.sin(phi1)*(wigner.W3J(Jf,1,Ja,mf,1,-1*ma)+wigner.W3J(Jf,1,Ja,mf,-1,-1*ma)))/\
        (w_fa**1.5)
        
        S3 = (-1)**(Jb-mb)*(np.cos(phi1)*wigner.W3J(Jb,1,Jf,mb,0,-1*mf)+\
        0.5*np.sin(phi1)*(wigner.W3J(Jb,1,Jf,mb,1,-1*mf)+wigner.W3J(Jb,1,Jf,mb,-1,-1*mf)))/\
        (w_fb**1.5)
        
        S4 = (-1)**(Jf-mf)*(np.cos(phi2)*wigner.W3J(Jf,1,Ja,mf,0,-1*ma)+\
        0.5*np.sin(phi2)*(wigner.W3J(Jf,1,Ja,mf,1,-1*ma)+wigner.W3J(Jf,1,Ja,mf,-1,-1*ma)))/\
        (w_fa**1.5)
        
        
        S[ii] = A_T*(2*Jf+1)*(S1*S2/Delta_1 + S3*S4/Delta_2) 

    # calculate Rabi rate in units of [s^-1]
    Omega_ab = abs(coeff*sum(S))
    Omega_single = abs(coeff*S[Pindex])
    
    # calculate total spontaneous emission decoherence rate [s^-1]
    Gamma1 = SpontEmission.SpontEmission(I1,B,w_1,phi1)
    Gamma2 = SpontEmission.SpontEmission(I2,B,w_2,phi2)
    
    Gamma_tot = Gamma1[2] + Gamma2[2]
    
    OmegaHz = Omega_ab/(2*pi)
    OmegasingleHz = Omega_single/(2*pi)
    t_pi = pi/(2*Omega_ab)
    t_pi_us = t_pi*1e6
    print('Rabi Rate (single level): %.2f Hz' %OmegasingleHz)
    print('Rabi Rate (all levels): %.2f Hz' % OmegaHz)
    print('Pi Time: %.2f us' % t_pi_us)
    print('Decoherence Rate: %.2f s^-1' % Gamma_tot)
    
    t = np.linspace(0,t_pi*5,500)
    t_us = t*1e6
    P_up = 0.5*(1+np.exp(-1*Gamma_tot*t)*np.cos(2*Omega_ab*t))
    
    plt.plot(t_us,P_up)
    plt.xlabel("Flopping Time ($\mu$s)")
    plt.ylabel("P$_{\uparrow}$")
    plt.axis([t_us[0],t_us[-1],0,1])
    
    return np.array([Omega_ab, Gamma_tot, Omega_ab/Gamma_tot])
                     
               
