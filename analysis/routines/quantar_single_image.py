# -*- coding: utf-8 -*-
"""
Created on Sun May 17 08:03:10 2015

@author: justinbohnet
"""
import os, importlib
import numpy as np
from numpy import sin, cos, pi
import matplotlib.pyplot as plt

import quantar_image
importlib.reload(quantar_image)

from skimage.feature import peak_local_max

x0 = 54.75
y0 = -9.0
num_img = 8
wrot = 189e3

fdir = os.getcwd()

qi = quantar_image.QuantarImage(x0=x0,y0=y0,fwall=wrot)
xyt = qi.read_file_range(fdir, 462, num_img)
xytr = qi.rot_frame(xyt)
xyt_bg = qi.read_file_range(fdir, 540, num_img)
xytr_bg = qi.rot_frame(xyt_bg)

img = qi.make_image(xytr, bck=xytr_bg, gfilter=0.2,
              im_range=quantar_image.im_extent(100))
              
coord = qi.get_ion_positions(img, extent=[100,150,150,100], min_distance=2.0,threshold_rel=0.4)
print(np.shape(coord))