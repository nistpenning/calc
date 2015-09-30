# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 14:43:24 2015

@author: rabi
"""

import os, shutil, importlib
import numpy as np
import numpy.random as rand
from numpy import pi, sqrt, sin
import matplotlib.pyplot as plt

import hfGUIdata as hf
import plot_style as ps
import plot_model_fit as pt
importlib.reload(ps)

samples = 1000
num_bins = 30
num_tip = 4

#create reference histogram
data_ref = rand.normal(size=samples)
vals_ref, bins, patches = plt.hist(data_ref,bins=num_bins,normed=True)
plt.show()
plt.close()

#create fake data
tipping_angle = np.linspace(-1.0,1.0,num=num_tip)
data = np.array([rand.normal(loc=i, size=samples) for i in tipping_angle])
for row in data:
    plt.hist(row, bins=num_bins, normed=True)
plt.show()

#realized the problem in my code -- the bin edges are getting changed, need to have constant bin edges
data_hist = []
for row in data:
    vals = plt.hist(row, bins=num_bins,normed=True)[0]
    print(vals)
    data_hist.append(vals)
data_hist = np.array(data_hist)
diff_hist = np.array([(sqrt(vals_ref) - sqrt(row)) for row,i in enumerate(data_hist) ])
sq_hist = diff_hist**2
h_dist = np.array([np.sum(row) for row in sq_hist])
plt.close()

pt.plot_polyfit(tipping_angle, h_dist,np.array([0,0,1]),hold=np.array([True,True,False]))
plt.show()
plt.close(0)

for row in diff_hist:
    plt.plot(bins[:-1],row)