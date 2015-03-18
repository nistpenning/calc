# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 08:23:21 2015

@author: justinbohnet
"""

import os
import numpy as np
from numpy import pi
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from scipy import fftpack

from skimage import feature
from skimage.morphology import watershed

#setup the info for the data
file_time = 0.2
fwall = 54e3
first_file = "00001220.dat"
found_flag = 0
print_flag = 0
num_to_read = 5

x0 = -4.0
y0 = -17.0
bins = 120

xLab = np.empty(0)
yLab = np.empty(0)
tLab = np.empty(0)

file_list = os.listdir(os.getcwd())
num_read = 0

for name in file_list:
    if name == first_file:
        found_flag = 1
    if found_flag == 1 and num_read < num_to_read:
        with open(name, 'rb') as f:
            norf = f.read(16)
            if print_flag: print(norf)
            
            N_reps_per_scan, = np.fromfile(f, dtype=np.uint32, count=1)
            if print_flag: print(N_reps_per_scan)
            
            rep_boundary_message = f.read(14)
            if print_flag: print(rep_boundary_message)
            
            N_this_rep, = np.fromfile(f, dtype=np.uint32, count=1)
            if print_flag: print(N_this_rep)
            
            if N_this_rep > 0:
                #x = np.fromfile(f, dtype=np.float64, count=N_this_rep, sep="")
                x = np.fromfile(f, dtype=np.float64, count=N_this_rep, sep="")
                y = np.fromfile(f, dtype=np.float64, count=N_this_rep, sep="")
                t = np.fromfile(f, dtype=np.float64, count=N_this_rep, sep="")
            
        xLab = np.append(xLab,x + x0)
        yLab = np.append(yLab,y + y0)
        tLab = np.append(tLab,t + file_time * num_read)
        
        num_read += 1
  
#counts, xedges, yedges, image = plt.hist2d(xLab,yLab,bins=120,
#                                           cmap = mpl.cm.Blues, norm=mpl.colors.LogNorm())
counts, xedges, yedges, LabFrame = plt.hist2d(xLab,yLab,bins=bins, cmap = mpl.cm.Blues,
                                           normed=True)
plt.show(LabFrame)

phaseOfWall = 2*pi*fwall* tLab

xRot = xLab * np.cos(phaseOfWall) + yLab * np.sin(phaseOfWall)
yRot = yLab * np.cos(phaseOfWall) - xLab * np.sin(phaseOfWall)

#Make Rotating Frame Image
counts, xedges, yedges, RotFrame = plt.hist2d(xRot,yRot,bins=120,
                                           cmap = mpl.cm.Blues,normed=True)
extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
counts_filter = ndi.gaussian_filter(counts,0.7)
RotFrame = plt.imshow(counts_filter,extent=extent)
plt.show(RotFrame)
#RotFrame = plt.imshow(exposure.equalize_adapthist(counts_filter),extent=extent)
#plt.show(RotFrame)

#make binary figure for ion counting
#mask = np.where(counts_filter> np.max(counts_filter)/2.35, 1, 0)
#mask_fig = plt.imshow(mask, extent=extent, cmap = mpl.cm.Blues)
#plt.show(mask_fig)
#use watershed to seperate?
coordinates = feature.peak_local_max(counts_filter, min_distance=1.75,threshold_rel=0.4)
plt.autoscale(True)
plt.imshow(counts_filter, cmap = mpl.cm.Blues, origin='lower')
plt.plot(coordinates[:, 1], coordinates[:, 0], 'r.')
plt.axis([40,80,40,80])
plt.show()


#all the unused stuff
if False:
    sx = ndi.sobel(counts_filter, axis=0, mode='constant')
    sy = ndi.sobel(counts_filter, axis=1, mode='constant')
    sobel = np.hypot(sx,sy)
    derv_fig = plt.imshow(sobel, extent=extent, cmap = mpl.cm.Blues)
    plt.show(derv_fig)
    
    smask = np.where(sobel> np.max(counts)/1.5, 1, 0)
    morphology = ndi.binary_fill_holes(smask)
    smask_fig = plt.imshow(morphology, extent=extent, cmap = mpl.cm.Blues)
    plt.show(smask_fig)


if False:
    F1 = fftpack.fft2(counts)
    F2 = fftpack.fftshift(F1)
    PSD2d = np.log10(np.abs(F2)**2)
    
    plt.imshow(PSD2d[25:75,25:75], origin='lower')
    
    
def count_ion(image, extent):
    sx = ndi.sobel(image, axis=0, mode='constant')
    sy = ndi.sobel(image, axis=1, mode='constant')
    sobel = np.hypot(sx,sy)
    
    smask = np.where(sobel > np.max(counts)/3.5, 1, 0)
    smask_fig = plt.imshow(ndi.binary_fill_holes(smask), extent=extent, 
                           cmap = mpl.cm.Blues)
    plt.show(smask_fig)
    
    #use watershed to seperate?
    #local_maxi = peak_local_max(counts_filter, indices=False, labels=counts_filter)
    
    label_im, nb_labels = ndi.label(smask)
    
    return count_ion


#find ions using blob determinatant of hessian
#but I couldn't get this to find anything
#blobs = feature.blob_doh(sobel, threshold=0.01, min_sigma=0.1)
#fig,ax = plt.subplots(1,1)
#ax.imshow(sobel, cmap = mpl.cm.Blues, origin='lower')
#ax.set_title("Blob method")
#for blob in blobs:
#    y,x,r = blob
#    c = plt.Circle((x,y),r,fill=False)
#    ax.add_patch(c)
#plt.show()