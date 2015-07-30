# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 11:09:48 2015

@author: jgb
"""
from __future__ import division
from scipy.constants import pi
import numpy as np
import scipy.ndimage as ndi
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

from skimage.feature import peak_local_max

class quantar_image:
    """A class for images from quantar .dat files
    **parameters**, **types**, **return** and **return types**::

    :param x0:
    :param y0: 
    :param num_to_read: 
    :param file_time: time of each file
    :param fwall: rotating wall freq kHz
    :param first_file: string of first file to read
    :return: return description
    :rtype: the return type description    
    
    """
    bins = 250

    #the constructor stores the data from the .dat files in the class
    def __init__(self, x0, y0, num_to_read, file_time, fwall, first_file):
        #data for creating images
        self.x0 = x0
        self.y0 = y0
        self.fw = fwall
        
        self.rot_image = None
        self.extent = [0,0,0,0]
        
        #data for getting back to the raw data
        self.file_time = file_time
        self.first_file = first_file #the number of the first .dat file for data
        self.num_to_read = num_to_read #the number of files used to create raw data
        
        #initialize background histogram
        self.bckgnd = np.zeros((self.bins,self.bins))
        
        found_flag = 0
        print_flag = 0
        
        self.x = np.empty(0)
        self.y = np.empty(0)
        self.t = np.empty(0)
        
        file_list = os.listdir(os.getcwd())
        num_read = 0   
    
        for name in file_list:
            if name == self.first_file:
                found_flag = 1
            if found_flag == 1 and num_read < self.num_to_read:
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
                    
                self.x = np.append(self.x,x + x0)
                self.y = np.append(self.y,y + y0)
                self.t = np.append(self.t,t + self.file_time * num_read)
                
                num_read += 1
        '''attempt to scale the data to correct xy distances
        conversion from quantar data to um is (53mm/60)*(97um/38.5mm)
        based on the size of the cloud on 032420125 and the known size of 
        the qtxyt3 centering cicle
        '''
        conversion = (53/60.0)*(97.0/38.5)
        self.x = self.x*conversion
        self.y = self.y*conversion
    
    def crop_image(self, image, c):
        lx, ly = image.shape
        return image[lx/2*(1-c) : lx/2*(1+c), ly/2*(1-c): ly/2*(1+c)]
    
    def set_background_hist(self,image):
        """takes an instance of the class, stores it as a background"""
        xLab = image.x
        yLab = image.y
        phaseOfWall = 2*pi*image.fw* image.t

        xRot = xLab * np.cos(phaseOfWall) + yLab * np.sin(phaseOfWall)
        yRot = yLab * np.cos(phaseOfWall) - xLab * np.sin(phaseOfWall)
        
        #Make Rotating Frame Image
        counts_background, xedges, yedges, RotFrame = plt.hist2d(xRot,yRot,bins=self.bins,
                                                   cmap = mpl.cm.Blues,normed=False)
        self.bckgnd = counts_background*self.num_to_read/image.num_to_read
    
    def make_lab_image(self,  im_range=None, gfilter=0.0):
        xLab = self.x
        yLab = self.y
        
        counts, xedges, yedges, LabFrame = plt.hist2d(xLab,yLab,bins=self.bins, cmap = mpl.cm.Blues,
                                           normed=False)
        extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
        counts_filter = ndi.gaussian_filter(counts-self.bckgnd, gfilter)
        #LabFrame = plt.imshow(counts_filter-self.bckgnd, cmap = mpl.cm.Blues,
        #                      vmin = 0, vmax = np.max(counts_filter))
        if im_range!=None: plt.axis(im_range)
        plt.xlabel("x [$\mu$m]")
        plt.ylabel("y [$\mu$m]")
        plt.show(LabFrame)
        return counts_filter,extent
        
    def make_rot_image(self, gfilter=0.0):
        xLab = self.x
        yLab = self.y
        phaseOfWall = 2*pi*self.fw* self.t

        xRot = xLab * np.cos(phaseOfWall) + yLab * np.sin(phaseOfWall)
        yRot = yLab * np.cos(phaseOfWall) - xLab * np.sin(phaseOfWall)
        
        #Make Rotating Frame Image
        counts, xedges, yedges, RotFrame = plt.hist2d(xRot,yRot,bins=self.bins,
                                                   cmap = mpl.cm.Blues,normed=False)
        extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
        counts_filter = ndi.gaussian_filter(counts-self.bckgnd,gfilter)
        RotFrame = plt.imshow(counts_filter,extent=extent, cmap = mpl.cm.Blues,
                              vmin = 0, vmax = np.max(counts_filter))
        plt.xlabel("x [$\mu$m]")
        plt.ylabel("y [$\mu$m]")
        plt.show(RotFrame)
        self.rot_image = counts_filter
        self.extent = extent
        
    def get_ion_positions(self):
        self.coordinates = peak_local_max(self.rot_image, min_distance=3.0,threshold_rel=0.4)
        
    def show_rot_image(self,im_range,low_threshold= 0):
        if np.size(self.rot_image) == 0:
            print("No rotating frame image, must make_rot_image() first")
        else:
            image = np.copy(self.rot_image)
            image[image < low_threshold] = 0
            RotFrame = plt.imshow(image,extent=self.extent, cmap = mpl.cm.Blues,
                              vmin = 0, vmax = np.max(image))
            if im_range!=None: plt.axis(im_range)
            plt.xlabel("x [$\mu$m]")
            plt.ylabel("y [$\mu$m]")
            plt.show(RotFrame)


######### useful functions but shouldn't be stored in the class  ########
def im_extent(mag):
    return np.array([-mag,mag,-mag,mag])
            
        