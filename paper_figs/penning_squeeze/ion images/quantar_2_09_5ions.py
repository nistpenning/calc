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
name = "2_09_5ions.png"
name = False
x0 = -20.5
y0 = -32.5
num_img = 20
num_bck_img = 0
wrot = 53e3
data_num = 500
bck_num = 530

intensity_range = [0,0.00003]
intensity_range = 'auto'
intensity_range = [1e-5,5e-5]

color_map = 'auto'
color_map = quantar_image.bluehot_cmap


base_path = os.getcwd()
fdir = os.path.normpath("/Volumes/688/Public/penning_britton/dailyLabBookFiles/2015/20150209/QuadWall_5brightIonCrystal_53kHz")

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
    img = qi.make_fig_image(xytr, gfilter=0.2, im_range=quantar_image.im_extent(110),
                        int_range=intensity_range, cmap=color_map, fig_axis=False, save_name=name)

coord = qi.get_ion_positions(img, extent=[90,160,90,160], min_distance=1.5,
                             threshold_rel=0.3)
                             
print("File for creation: "+str(fdir)+", total integration time: {} s".format(xyt[-1][2]- xyt[1][2]))
print("x0: {}, y0: {}, File#: {}, # of files: {}".format(x0,y0,data_num,num_img))
print("Ion number counting: {}".format(np.shape(coord)))