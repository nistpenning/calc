# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 09:08:26 2015

@author: jgb
"""

import numpy as np
from numpy import pi, exp, sqrt, sin, cos
import matplotlib.pyplot as plt

import hfGUIdata
import plot_tools_jgb as pt
import squeeze_func_time as squ

# setup calculation
J = 1.0
N = 100
G_el = 0.0
G_ud = 0.1
G_du = 0.1
G_r = G_ud + G_du

g = (G_ud-G_du)/4.
G = (G_el + G_r)/2.

tsq = (3**(1/6.0))/(2*N**(2/3.0))
psi = np.linspace(0,pi,num=100)

# Key functions
def Phi(t, J, G_ud, G_du):
    theta = pi/2.
    g = (G_ud-G_du)/4. 
    G_r = G_ud + G_du
    s = 2*J + 2j*g
    r= G_ud*G_du
    
    a = t*sqrt(s**2-r)
    
    out = exp(-G_r*t/2.)*( cos(a) + 0.5*(G_r+4j*J*cos(theta))*t*np.sinc(a) )
    return out
    
def Psi(t, J, G_ud, G_du):
    theta = pi/2.
    g = (G_ud-G_du)/4. 
    G_r = G_ud + G_du
    s = 2*J + 2j*g
    r= G_ud*G_du
    
    a = t*sqrt(s**2-r)
    
    out = exp(-G_r*t/2.)*( cos(theta)*cos(a) + 0.5*(2j*s - 4*g - G_r*cos(theta))*t*np.sinc(a) )
    return out

# collective spin length
Sx = 0.5*exp(-G*tsq) * (Phi(tsq,J,G_ud, G_du))**(N-1)
Sx = np.real(Sx)

# pm correlations
Cmz = 0.5*exp(-G*tsq) * Psi(tsq, -J ,G_ud, G_du) * (Phi(tsq, -J, G_ud, G_du))**(N-2)
Cpz = 0.5*exp(-G*tsq) * Psi(tsq, J, G_ud, G_du) * (Phi(tsq, J, G_ud, G_du))**(N-2)

Cpp = 0.25*exp(-2*G*tsq) * (Phi(tsq, 2*J, G_ud, G_du))**(N-2)
Cmm = 0.25*exp(-2*G*tsq) * (Phi(tsq, -2*J, G_ud, G_du))**(N-2)
Cpm = 0.25*exp(-2*G*tsq) * (Phi(tsq, 0.0, G_ud, G_du))**(N-2)
Cmp = Cpm

# coorediate correlations
CyyA = -0.25*(Cpp+Cmm-Cpm-Cmp)
CyzA = -(0.25j)*(Cpz-Cmz)
CzyA = CyzA

SyyA = N*(N-1)*CyyA + N/4.
SyzA = N*(N-1)*CyzA + 0.5j*N*Sx
SzyA = N*(N-1)*CzyA - 0.5j*N*Sx
SzzA = N/4.

# Std Devs
DyyA = sqrt(SyyA)
DzzA = sqrt(SzzA)
DyzA = sqrt(SyzA)
DzyA = sqrt(SzyA)

# calc squeezing
z = (SzyA+SyzA)/(SzzA-SyyA)
opt_squ_angle = np.real(0.5*np.arctan(z))
Jz_std = sqrt(cos(-psi)**2*SzzA + sin(-psi)**2*SyyA + sin(-psi)*cos(-psi)*(SyzA+SzyA))
R = (Jz_std/(sqrt(N)/(2.0)))**2  # reduction in spin noise variance
C = 2*Sx  # fringe contrast

Xi2 = R * (C)**(-2)  # Ramsey squeezing parameter squared

# display
plt.plot(psi, sqrt(Xi2))
plt.plot(psi,sqrt(R))

