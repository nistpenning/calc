# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 12:30:50 2015

@author: jgb
"""

import numpy as np
from numpy import pi, sqrt

N = 100.0
var_PN = 1/N

var_coll_dephas = 0.0022  # at 1ms

phot_per_ion = 15.0  # for detection time 15 ms

var_SN = (1/sqrt(phot_per_ion*N/2))**2

var_TN = (0.002*15.+0.242)  #var relative to PSN, at 15ms detection time, from 7/21, low perp
var_TN = var_TN * var_SN

print("Class. detection noise: {:.3g} radian_sq, rel. PN 100 ions: {:.3g}".format(var_TN, var_TN/var_PN))
print("Photon shot noise: {:.3g} radian_sq, rel. PN 100 ions: {:.3g}".format(var_SN, var_SN/var_PN))