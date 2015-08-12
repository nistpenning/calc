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

x0 = 41.5
y0 = -0.5
num_img = 40
wrot = 190e3
int_time = 0.1  # sec

fdir = os.getcwd()

qi = quantar_image.QuantarImage(x0=x0,y0=y0,fwall=wrot)
xyt = qi.read_file_range(fdir, 910, 35)
xytr = qi.rot_frame(xyt)
xyt_bg = qi.read_file_range(fdir, 5900, 35)
xytr_bg = qi.rot_frame(xyt_bg)
qi.make_image(xytr,im_range=quantar_image.im_extent(100))

#%%
bck_first_file = "00005900.dat"
    
bck = qi.QuantarImage(x0, y0, num_img, int_time, wrot, bck_first_file)

bck.make_lab_image()

#%%
first_file = "00005720.dat"
a = qi.quantar_image(x0, y0, num_img, int_time, wrot, first_file)

a.set_background_hist(bck)

lab = a.make_lab_image()


#%%
first_file = "00006440.dat"
a = qi.quantar_image(x0, y0, num_img, int_time, wrot, first_file)

a.set_background_hist(bck)

lab = a.make_lab_image()

#%%
first_file = "00000910.dat"
a = qi.quantar_image(x0, y0, num_img, int_time, wrot, first_file)

a.set_background_hist(bck)

lab = a.make_lab_image()

#%%
img = a.make_rot_image(gfilter=0.35)
ex = [-100,100,-100,100]
a.show_rot_image(ex,low_threshold=2.0)

#plt.savefig("image_out.png",dpi=300)

"""
img = np.copy(a.rot_image)
img = a.crop_image(img, 0.3)

img = a.fix_contrast(img, debug=True)

ax = plt.axes()
plt.imshow(img, cmap=plt.cm.gray)
ax.set_axis_off()
plt.savefig("image_out.png",dpi=300)
"""

#%%
a.get_ion_positions(extent=[90,160,160,90],min_distance=2.2,threshold_rel=.35)
print np.shape(a.coordinates)