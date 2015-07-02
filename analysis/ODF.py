# -*- coding: utf-8 -*-


from __future__ import division

import numpy as np
from numpy import cos, sin, arctan2, pi
import sys
from scicons import c, m_Be, ge, epsilon0, hbar, mu_b, mu_n
import scipy.constants
import spont_emission as spe

sys.path.append("c:\\Users\\bsawyer\\Documents\\Python Scripts\\subroutines")
import wigner
reload(wigner)

gIp = -0.784955 #Be nuclear g-factor
AhfS = -625.008837e6 * 2*pi #Hyperfine coefficient
AhfP = -118.6e6 * 2*pi #Hyperfine coefficient for P state
mI = +3./2. #nuclear spin projection of Be+ used for experiments


def Plevelcalc(B=4.4588):

    omega12 = 31928.744 * (2*pi*c*100) #Center of P1/2 state in cm^-1
    omega32 = 31935.320 * (2*pi*c*100) #Center of P3/2 state in cm^-1
    Aik1 = 1.1285e8 #A-coefficient of J=1/2 P state
    Aik2 = 1.1292e8 #A-coefficient of J=3/2 P state

    J = np.array([3/2, 3/2, 3/2, 3/2, 1/2, 1/2])
    mJ = np.array([3/2, 1/2, -1/2, -3/2, 1/2, -1/2])
    omega = [omega32,omega32,omega32,omega32,omega12,omega12]
    Aik = np.array([Aik2, Aik2, Aik2, Aik2, Aik1, Aik1])

    E = np.zeros(6)
    for ii in range(0,6):
        g = 1 + (J[ii]*(J[ii]+1) + .5*1.5 - 2) / (2*J[ii]*(J[ii]+1))

        if ii>4:
            E[ii] = (g*mJ[ii]*mu_b - gIp*mI*mu_n)*B/hbar + AhfP*mI*mJ[ii] + omega[ii]
        else:
            E[ii] = (g*mJ[ii]*mu_b - gIp*mI*mu_n)*B/hbar + omega[ii]

    E = np.diag(E)
    melem = np.sqrt(2)*mu_b*B/(3*hbar)*(1+ge)

    E[1][4] = melem
    E[4][1] = melem
    E[2][5] = melem
    E[5][2] = melem

    levels = sorted(np.linalg.eig(E)[0],reverse=1)

    return J, mJ, levels, Aik

def calc_wL(B=4.4588, L_c=0.5*626.2252e-9, L_ODF=0.5*626.2670e-9):

    """
    A simple function that can correctly calcuate the laser frequency and detuning
    from the cooling transition given the wavelength of the red light on the
    angstrom wavementer
    """
    Plevels = Plevelcalc(B=B)
    mi=0.5
    omega0 = -1*(mi*ge*mu_b + mI*gIp*mu_n)*B/hbar + AhfS*mI*mi #S1/2,

    '''
    okay found that the wavelengths we use in the lab don't match up to the
    Plevel calcuation for B = 4.4588 T.

    For example, L_cool_lab = 626.2252 nm, and Plevels is 626.1841 nm

    so to get the laser frequency, I will calculate the detuning in the lab, then
    use that detuning from Plevels for the theory calculation
    '''
    L_c = 0.5*626.2252e-9
    L_ODF = 0.5*626.2670e-9
    Delta_ODF = c*(1/L_ODF - 1/L_c) #Hz

    w_cool = Plevels[2][0] - omega0#957.519919e12 #at B=4.4584 T
    w_L = w_cool + (2*pi*Delta_ODF)

    return w_L, 2*pi*Delta_ODF


def ACstarkshift(p=1, #polarization -1, 0 , or 1
                 I=1e4, #intensity W/m^2
                 mi=0.5, #initial m_j level
                 omegaL=6015492587866242.0, #laser frequency s^-1
                 B=4.4588): #Bfield tesla

    #calculate Paschen-Bach shift of each level
    omega0 = -1*(mi*ge*mu_b + mI*gIp*mu_n)*B/hbar + AhfS*mI*mi #S1/2,

    Plevels = Plevelcalc(B)

    #Remember that ACSS = (1/4)*alpha*|E|**2 (alpha is polarizability).  The 1/4
    #factor comes from time-averaging the electric field amplitude
    coeff = -3*pi*(c**2)*I

    #S is the summand
    S = 0

    for ii in range(0,6):


        #omegaik is absolute energy difference between level i and k
        omegaik = Plevels[2][ii] - omega0
        #print omegaik
        #S1 is numerator
        S1 = Plevels[3][ii]*(2*Plevels[0][ii]+1)*\
        wigner.W3J(0.5,1,Plevels[0][ii],mi,p,-1*Plevels[1][ii])**2

        #S2 is denominator
        S2 = omegaik**2 * (omegaik**2 - omegaL**2)
        #print S2

        S = S + S1/S2

        #print S
    #This spits out shift in angular frequency
    return coeff*S/hbar

def ACstarkshiftTRANS(p=1,I=1e4,omegaL=6015492587866242.0,B=4.4588):

    return ACstarkshift(p,I,1/2,omegaL,B) - ACstarkshift(p,I,-1/2,omegaL,B)


def ODF(I=1e4, thetaR=10.0, omegaL=6015492587866242.0):

    Aup = ACstarkshift(0, I, 1/2, omegaL)
    Bup = 0.5*(ACstarkshift(1, I, 1/2, omegaL)+ACstarkshift(-1, I, 1/2, omegaL))

    Adown = ACstarkshift(0,I,-1/2,omegaL)
    Bdown = 0.5*(ACstarkshift(1,I,-1/2,omegaL)+ACstarkshift(-1, I, -1/2, omegaL))

    phiup = np.arcsin(np.sqrt(Aup/(Aup-Bup))) #calculate null angle
    phidown = np.arcsin(np.sqrt(Adown/(Adown-Bdown)))

    k = omegaL/c

    Fup = 4*hbar*k*np.sin(thetaR*pi/180)*(Aup*np.cos(phiup)**2-Bup*np.sin(phiup)**2)
    Fdown = 4*hbar*k*np.sin(thetaR*pi/180)*(Adown*np.cos(phidown)**2-Bdown*np.sin(phidown)**2)

    return (Fup-Fdown)/2.0

#def spinmotion(N = 100, # number of Be+
#               I = 1e4, # ODF intensity in W/m^2
#               thetaR = 10, # ODF beam half-angle
#               fz = 1570e3, # COM frequency (Hz)
#               omegaL = 6015492587866242.0, # ODF laser frequency (s^-1)
#               T = 0.5e-3, # COM temperature (K)
#               Tarm = 100e-6, # spin echo arm time (s)
#               t_pi = 82e-6): # pi time (s)
#
#    omegaz = 2*pi*fz
#    delta = 2*pi*np.linspace(-3./Tarm,3./Tarm,100)
#
#    dk = 2*omegaL/c * np.sin(thetaR*pi/180.)
#    z0 = np.sqrt(hbar/(N*2*))
#    F0 = ODF(I=I,thetaR=thetaR,omegaL=omegaL)

def ACSS_I_cal(diffACSSpi, wL=6015492587866242.0,B=4.4588):
    Plevels = Plevelcalc(B)

    #calculate Paschen-Bach shift of each level
    mi = 0.5
    omegau = -1*(mi*ge*mu_b + mI*gIp*mu_n)*B/hbar + AhfS*mI*mi #S1/2,
    omegauu = Plevels[2][1] - omegau
    #lower
    mi = -0.5
    omegaL = -1*(mi*ge*mu_b + mI*gIp*mu_n)*B/hbar + AhfS*mI*mi #S1/2,
    omegaLL = Plevels[2][2] - omegaL

    up_term = Plevels[3][1]/(omegauu**2*(omegauu**2-wL**2))
    down_term = Plevels[3][2]/(omegaLL**2*(omegaLL**2-wL**2))

    I = diffACSSpi/(up_term-down_term)*(-hbar/(2*pi*c**2))

    return I

def waist_from_piACSS(ACSS, Pow, wL=6015492587866242.0, B=4.4588):
    I_fit = ACSS_I_cal(ACSS, wL=wL, B=B)
    return np.sqrt(2*Pow/(pi*I_fit))

def OATjgb(N, Jbar_at_1kHz, psi, int_time):
    '''    
    calclates the var of Jz as a func of the angle for a 1 axis twisting 
    Hamiltonian.  
    assumes that we are driving with the Penning trap ODF, and making complete
    loops in phase space, so the detuning from the COM mode is set by the
    interaction time, thus the squeezing strength chi is calculated in the 
    function
    '''
    t=int_time
    delta = 2e-3/t
    chi = 2*Jbar_at_1kHz/N/(delta)
    Acoef = 1-cos(2*chi*t)**(N-2)
    Bcoef = 4*sin(chi*t)*cos(chi*t)**(N-2)
    delt = 0.5*np.arctan2(Bcoef,Acoef)
    varJz = N/4.0*(1+(N/4.0-0.25)*(Acoef - np.sqrt(Acoef**2+Bcoef**2)*cos(2*(psi+delt))))

    return varJz, delt

def IsingCalc(ACSS_l, ACSS_u, wz):
    '''
    gives parameters needed for Ising model calcuations from measured ACSS
    Units:
        ACSS_l -- per sec
        
        ACSS_u -- per sec
        
        wz -- per sec
        
    '''

    Lcooling = 626.2252e-9  # wavemeter reading of red lasers
    LODF = 626.2670e-9
    #  get the laser wavelength, handling the mismatch between Plevels calc and wavemeter
    wL, Delta_ODF = calc_wL(B=4.4588, L_c=0.5*Lcooling, L_ODF=0.5*LODF)

    I_fit_u = ACSS_I_cal(ACSS_u, wL=wL)
    I_fit_L = ACSS_I_cal(ACSS_l, wL=wL)
    I_fit_avg = (I_fit_u+I_fit_L)/2.

    #  calc angles from exp offset
    phip_u = 57.0*2  # degress, the offset is found from the best linear polarization
    phip_L = 32.0*2  # degrees

    #  spont emission
    G_u = spe.spont_emission(I=I_fit_u, omegaL=wL, phip=phip_u*pi/180.)
    G_L = spe.spont_emission(I=I_fit_L, omegaL=wL, phip=phip_L*pi/180.)

    F0 = ODF(I=I_fit_avg, thetaR=10.0, omegaL=wL)

    def Jbar(delta): return F0**2/(4*hbar*m_Be*wz*(2*pi*delta*1e3))

    Jbar_1kHz = Jbar(1)

    return G_u[2]+G_L[2], Jbar_1kHz, F0
