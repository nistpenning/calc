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

#%%
name = "10_14_124ions.png"
x0 = 48.25
y0 = -13.5
num_img = 20
num_bck_img = 60
wrot = 187e3
data_num = 1080
bck_num = 1200

intensity_range = 'auto'
intensity_range = [0,0.00006]

color_map = 'auto'
color_map = quantar_image.bluehot_cmap


base_path = os.getcwd()
fdir = os.path.normpath("/Users/jgb/Data/20151014/image")

qi = quantar_image.QuantarImage(x0=x0,y0=y0,fwall=wrot)

xyt = qi.read_file_range(fdir, data_num, num_img)
xytr = qi.rot_frame(xyt)
xyt_bg = qi.read_file_range(fdir, bck_num, num_bck_img)
xytr_bg = qi.rot_frame(xyt_bg)

#%%
if color_map == 'auto':
    img = qi.make_fig_image(xytr, bck=xytr_bg, gfilter=0.2, im_range=quantar_image.im_extent(110),
                        int_range=intensity_range)
else:
    img = qi.make_fig_image(xytr, bck=xytr_bg, gfilter=0.2, im_range=quantar_image.im_extent(110),
                        int_range=intensity_range, cmap=color_map, save_name=name)
    img = qi.make_fig_image(xytr, bck=xytr_bg, gfilter=0.2, im_range=quantar_image.im_extent(110),
                        int_range=intensity_range, cmap=color_map, fig_axis=False, save_name="noaxis_"+name)

print("x0: {}, y0: {}, File#: {}".format(x0,y0,data_num))

coord = qi.get_ion_positions(img, extent=[90,160,90,160], min_distance=1.5,
                             threshold_rel=0.14)

print(np.shape(coord))
