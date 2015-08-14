# -*- coding: utf-8 -*-
"""
Created on Sun May 17 08:03:10 2015

@author: justinbohnet
"""
import os, importlib
import numpy as np
from numpy import sin, cos, pi
import matplotlib.pyplot as plt
import matplotlib as mpl

import quantar_image
importlib.reload(quantar_image)

from skimage.feature import peak_local_max

#%%
x0 = 41.5
y0 = -1.1
num_img = 60
wrot = 190e3
data_num = 7000

fdir = os.getcwd()

qi = quantar_image.QuantarImage(x0=x0,y0=y0,fwall=wrot)
xyt = qi.read_file_range(fdir, data_num, num_img)
xytr = qi.rot_frame(xyt)
xyt_bg = qi.read_file_range(fdir, 5950, num_img)
xytr_bg = qi.rot_frame(xyt_bg)

#%%

img = qi.make_image(xytr, gfilter=0.2,
              im_range=quantar_image.im_extent(100), cmap=quantar_image.bluehot_cmap)
print("x0: {}, y0: {}, File#: {}".format(x0,y0,data_num))


              
coord = qi.get_ion_positions(img, extent=[90,160,90,160], min_distance=2.0,
                             threshold_rel=0.4)
print(np.shape(coord))