# -*- coding: utf-8 -*-
"""
Created on Thu Dec 04 09:22:22 2014

@author: bsawyer
"""

from __future__ import division

import numpy as np
import sys
from scicons import *

#sys.path.append("c:\\Users\\bsawyer\\Documents\\Python Scripts\\subroutines")
import wigner
import ODF

AhfS = -625.008837e6 * 2*pi #Hyperfine coefficient 
AhfP = -118.6e6 * 2*pi #Hyperfine coefficient for P state
#ge is electron g-factor
gIp = -0.784955 #Be nuclear g-factor
mI = 3./2. #nuclear spin projection used for experiments

omega12 = 31928.744 * (2*pi*c*100) #Center of P1/2 state in cm**-1
omega32 = 31935.320 * (2*pi*c*100) #Center of P3/2 state in cm**-1
A32 = 1.1292e8  #A-coefficient of J=3/2 P state



def Aij(lamb,J,mi,mj,omegaL=6015492587866242.0,B=4.4588):
    
    mk = lamb + mi
        
    if J < 1:
        index = int(4.5-mk) #find index of (J,mJ) level
    else:
        index = int(1.5-mk)  #find index of (J,mJ) level

#Plevelcalc outputs a 4-column matrix [J mJ frequency lifetime]
    Plevels = ODF.Plevelcalc(B)

#calculate Paschen-Bach shift of each level
    omega0 = (-1*mi*ge*mu_b - mI*gIp*mu_n)*B/hbar + AhfS*mI*mi #S1/2,


# Calculate detuning of laser from the ik transition

    if index > (np.size(Plevels[2])-1):
        deltaik = 1e30 
    else:
        deltaik = Plevels[2][index] - omega0 - omegaL       

    return ((2*J+1)/deltaik)*wigner.W3J(1/2,1,J,mj,lamb+(mi-mj),-(lamb+mi))*\
    wigner.W3J(1/2,1,J,mi,lamb,-(lamb+mi))



def spont_emission(I=1e4, #laser intensity [W/m^2]
                  B=4.4588, #magnetic field [T]
                  omegaL=6015492587866242.0, #laser frequency [s^-1]
                  phip=1e6): #laser linear polarization [radians]
    
    Aacss = ODF.ACstarkshiftTRANS(0,I,omegaL,B)
    Bacss = 0.5*(ODF.ACstarkshiftTRANS(1,I,omegaL,B)+ODF.ACstarkshiftTRANS(-1,I,omegaL,B)) 
    #print Aacss, Bacss

    if phip > 1e5:
        phip = np.arcsin(np.sqrt(Aacss/(Aacss-Bacss)))
        #print phip
        
    J = np.array([0.5,1.5])
    
    blamb = np.array([np.sin(phip)/np.sqrt(2), np.cos(phip), -1.0*np.sin(phip)/np.sqrt(2)])
    #print blamb
    
    g = 1 + (1.5*2.5 + .5*1.5 - 2) / (2*1.5*2.5) # g for J=3/2 state
    
    omega00 = (-0.5*ge*mu_b - mI*gIp*mu_n)*B/hbar + AhfS*mI*0.5 #Zeeman shift of up state
    omegaik = ((g*(3/2)*mu_b - gIp*mI*mu_n)*B/hbar + omega32) - omega00    
    
    coeff = A32**2*(3*pi*c**2*I)/(2*hbar*omegaik**3)

    Eltot = 0.0
    Ramud = 0.0
    Ramdu = 0.0
    for ii in range(0,3):
        lamb = float(ii-1) #lamb runs from -1 to 1
    
        El1 = 0.0
        El2 = 0.0
        Ram1 = 0.0
        Ram2 = 0.0
        
        for jj in range(0,2):
        
            El1 = El1 + Aij(lamb,J[jj],1/2,1/2,omegaL,B) #up -> up scattering
            El2 = El2 + Aij(lamb,J[jj],-1/2,-1/2,omegaL,B) #down -> down scattering
                
            Ram1 = Ram1 + Aij(lamb,J[jj],1/2,-1/2,omegaL,B) #up -> down scattering
            Ram2 = Ram2 + Aij(lamb,J[jj],-1/2,1/2,omegaL,B) #down -> up scattering

    
        Eltot = Eltot + blamb[ii]**2*(El1 - El2)**2
        #print Eltot
        Ramud = Ramud + blamb[ii]**2*Ram1**2
        Ramdu = Ramdu + blamb[ii]**2*Ram2**2
    
    GammaEl = np.abs(coeff*Eltot)
    GammaRam = coeff*(np.abs(Ramud) + np.abs(Ramdu))
    Gamma = 0.5*(GammaEl + GammaRam)

    Gammaud = coeff*np.abs(Ramud)
    Gammadu = coeff*np.abs(Ramdu)
    
    return [GammaEl, GammaRam, Gamma, phip, Gammaud, Gammadu]