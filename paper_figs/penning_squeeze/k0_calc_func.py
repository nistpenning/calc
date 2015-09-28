# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 12:40:57 2015

@author: jgb
"""
import numpy as np

def k0_calc(m, std_dev_deg, N):
    pn_var = 1/float(N)
    k0 = (std_dev_deg*pi/180.0)**2/m
    return k0/pn_var
    