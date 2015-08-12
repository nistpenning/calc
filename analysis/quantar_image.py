# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 11:09:48 2015

@author: jgb
"""
from scipy.constants import pi
import numpy as np
import scipy.ndimage as ndi
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

import skimage
from skimage.feature import peak_local_max
import skimage.exposure


class NiQuantarFileError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

class QuantarImage:
    """A class for images from quantar .dat files
    **parameters**, **types**, **return** and **return types**::

    :param x0:
    :param y0: 
    :param num_to_read: 
    :param file_time: integration time for each file [s]
    :param fwall: rotating wall freq [kHz]
    :param first_file: string of first file to read
    :param fdir: absolute path to data folder (e.g. 'c:\\my_data')
    :return: return description
    :rtype: the return type description    
    
    """
    bins = 250

    def __init__(self, x0=0, y0=0,fwall=100e3 ):
        self.x0 = x0
        self.y0 = y0
        self.fw = fwall

        '''attempt to scale the data to correct xy distances
        conversion from quantar data to um is (53mm/60)*(97um/38.5mm)
        based on the size of the cloud on 032420125 and the known size of
        the qtxyt3 centering cicle
        '''
        self.scale_xy = (53/60.0)*(97.0/38.5)

    def read_niquantar_file(self, fnum, fdir, verbose=False):
        """Read binary file written by niquantar.exe
        :param fnum: file number
        :param fdir: absolute path to file directory
        :return: (x,y,t)
        """

        fname = str(fnum).zfill(8) + ".dat" # e.g. 00000001.dat
        fpath = os.path.join(fdir,fname)
        if verbose is True:
            print(fpath)
        try:
            fh = open(fpath, 'rb')
        except IOError:
            print("Error: File does not appear to exist.")
            return -2
        try:
            # read file marker
            norf = fh.read(16)
            if norf != b'NumOfRepsFollows': raise NiQuantarFileError(-1)
            # read number of reps per scan
            nrps = np.fromfile(fh, dtype=np.uint32, count=1)
            rep_boundary_message = fh.read(14)
            if rep_boundary_message !=  b'NextRepFollows':
                raise NiQuantarFileError(-3)
            n_this_rep = np.fromfile(fh, dtype=np.uint32, count=1)
            if n_this_rep < 0: raise NiQuantarFileError(-4)
        except IOError:
            print("Error reading norf"); return -5
        except NiQuantarFileError as e:
            print("NiQuantarFileError: {}".format(e)); return -6
        xyt = np.zeros((n_this_rep,3))
        try:
            if n_this_rep > 0:
                xyt[:,0] = np.fromfile(fh, dtype=np.float64, count=n_this_rep, sep="")
                xyt[:,1] = np.fromfile(fh, dtype=np.float64, count=n_this_rep, sep="")
                xyt[:,2] = np.fromfile(fh, dtype=np.float64, count=n_this_rep, sep="")
        except IOError:
            print("Error reading data"); return -7
        fh.close()
        return xyt

    def read_file_range(self, fdir, nfirst, n_to_read):
        """Read range of files from disk
        :param fdir: absolute path to directory [c:\\tmp]
        :param nfirst: first file number
        :param n_to_read: number to read
        :return: xyt
        """
        xyt = self.read_niquantar_file(nfirst, fdir)
        for fnum in range(nfirst+1,nfirst+n_to_read):
            xyt_new = self.read_niquantar_file(fnum, fdir)
            xyt = np.concatenate((xyt, xyt_new), axis=0)
        return xyt

    def rot_frame(self, xyt):
        """Return xyt in rotating frame"""
        x_lab = (xyt[:, 0] + self.x0)*self.scale_xy
        y_lab = (xyt[:, 1] + self.y0)*self.scale_xy
        phase_of_wall = 2*pi*self.fw*xyt[:, 2]

        x_rot = x_lab*np.cos(phase_of_wall) + y_lab * np.sin(phase_of_wall)
        y_rot = y_lab*np.cos(phase_of_wall) - x_lab * np.sin(phase_of_wall)

        return np.column_stack((x_rot, y_rot, xyt[:,2]))

    def crop_image(self, image, c):
        lx, ly = image.shape
        return image[lx/2*(1-c) : lx/2*(1+c), ly/2*(1-c): ly/2*(1+c)]
    
    def plot_img_and_hist(self, img, axes, bins=256):
        """Plot an image along with its histogram and cumulative histogram.

        """
        img = skimage.img_as_float(img)
        ax_img, ax_hist = axes
        ax_cdf = ax_hist.twinx()

        # Display image
        ax_img.imshow(img, cmap=plt.cm.gray)
        ax_img.set_axis_off()

        # Display histogram
        ax_hist.hist(img.ravel(), bins=bins, histtype='step', color='black')
        ax_hist.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))
        ax_hist.set_xlabel('Pixel intensity')
        ax_hist.set_xlim(0, 1)
        ax_hist.set_yticks([])

        # Display cumulative distribution
        img_cdf, bins = skimage.exposure.cumulative_distribution(img, bins)
        ax_cdf.plot(bins, img_cdf, 'r')
        ax_cdf.set_yticks([])

        return ax_img, ax_hist, ax_cdf
    
    def fix_contrast(self, img, debug=False):
        # normalize
        img = img / np.max(img)

        # Contrast stretching
        p2, p98 = np.percentile(img, (2, 98))
        img_rescale = skimage.exposure.rescale_intensity(img, in_range=(p2, p98))

        if debug:
            # Display results
            fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 5))
            ax_img, ax_hist, ax_cdf = self.plot_img_and_hist(img, axes[:, 0])
            ax_img.set_title('Low contrast image')
            y_min, y_max = ax_hist.get_ylim()
            ax_hist.set_ylabel('Number of pixels')
            ax_hist.set_yticks(np.linspace(0, y_max, 5))

            ax_img, ax_hist, ax_cdf = self.plot_img_and_hist(img_rescale, axes[:, 1])
            ax_img.set_title('Contrast stretching')
            # prevent overlap of y-axis labels
            fig.subplots_adjust(wspace=0.4)
            plt.show()
            
        return img_rescale

    def make_image(self, xyt, im_range=[-256,256,-256,256], gfilter=0.0):
        """plot image

        :param im_range: [-256,256,-256,256] is full range for Quantar
        :gfilter: ndi.gaussian_filter() argument
        :return: 2d histogram
        """
        
        #Make Rotating Frame Image
        plt.subplot(111, aspect='equal')
        ax = plt.gca()
        ax.grid(True)
        counts, xedges, yedges, RotFrame = plt.hist2d(xyt[:,0], xyt[:,1],
                                         bins=self.bins,
                                         cmap = mpl.cm.Blues, normed=False)
        extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
        counts_filter = ndi.gaussian_filter(counts,gfilter)
        RotFrame = plt.imshow(counts_filter,extent=extent, cmap = mpl.cm.Blues,
                              vmin = 0, vmax = np.max(counts_filter))
        plt.axis(im_range)
        plt.xlabel("x [$\mu$m]")
        plt.ylabel("y [$\mu$m]")
        plt.show(RotFrame)

        
    def get_ion_positions(self):
        self.coordinates = peak_local_max(self.rot_image, min_distance=3.0,threshold_rel=0.4)

    def show_rot_image(self, im_range, low_threshold= 0):
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

def main():
    fdir = "D:\\tmp\\20150730\\pic_2_load300_188kHz_rotation"
    qi = QuantarImage(x0=55,y0=-15,fwall=188e3)
    xyt = qi.read_file_range(fdir, 2780, 35)
    xytr = qi.rot_frame(xyt)
    xyt_bg = qi.read_file_range(fdir, 3100, 35)
    xytr_bg = qi.rot_frame(xyt_bg)
    qi.make_image(xytr)

if __name__ == "__main__":
    main()