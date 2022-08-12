#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 17:46:56 2022

@author: Pawel Janas (Python Python 3.9.12)

Description: This program performs continuum subtraction on NIRCam images from
JWST for detection of protostellar outflows.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import glob # for filepath handling
import sys

# Astropy:
from astropy.io import fits
from astropy.wcs import WCS
import reproject as rpj # needed for reproject module

# Sci-kit learn for KMeans clustering (outlier removal)
from sklearn.cluster import KMeans

sys.path.append('../code')
from reproject_combine import reproject

def open_fits(fname):
    """ Convenience function for reading in fits files, getting data and
    header information."""
    with fits.open(fname) as hdu:
        data = hdu['SCI'].data
        header = hdu['SCI'].header
    
    return data, header

def remove_outliers(c_dat, n_dat, n_clusters=8, max_labels=2):
    """
    Function for removing outliers using KMeans cluster algorithm.

    Parameters
    ----------
    c_dat : 1D array (masked)
        Masked array of flattened, masked continuum data.
    n_dat : 1D array (masked)
        Masked array of flattened, masked narrowband data.
    n_clusters : int, optional
        The number of clusters to form as well as the number of centroids 
        to generate for the KMeans algorithm. The default is 8.
    max_labels : int, optional
        The max label number (cluster number) to gather final data points 
        from. Higher values result it more outliers. Must be <= n_clusters.
        The default is 2.

    Returns
    -------
    cluster_data[:,0] : array
        1D array of continuum data with outliers removed.
    cluster_data[:,1] : array
        1D array of narrowband data with outliers removed.

    """
    if max_labels > n_clusters:
        print("n_clusters < max_labels. Make sure max_labels < \
n_clusters. Reverting to defaults...")
        n_clusters=8
        max_labels=2
        
    # Perform KMeans cluster algo
    cluster_data = np.ma.vstack((c_dat, n_dat)).T
    k_means = KMeans(n_clusters=n_clusters).fit(cluster_data)
    labels = k_means.labels_
    cluster_data = cluster_data[np.where(labels <= max_labels)]
    
    return cluster_data[:,0], cluster_data[:,1]

class ContinuumSubtract:
    def __init__(self, narrowband_file, continuum_file):
        self.narrowband_file = narrowband_file
        self.continuum_file = continuum_file
        
    def get_data_headers(self, get_data=True, get_headers=True):
        c_dat, c_header = open_fits(self.continuum_file)
        n_dat, n_header = open_fits(self.narrowband_file)
        if get_data and get_headers:
            return c_dat, c_header, n_dat, n_header
        elif get_data:
            return c_dat, n_dat
        elif get_headers:
            return c_header, n_header
        else:
            print("You didn't select get_data or get_headers so the function \
returned None")
        
    def reproject_continuum(self):
        """ Reproject continuum data onto narrowband data and also resizes them 
        to match. Resizes to the smaller image for memory conservation. Returns
        list of reprojected data."""
        list_files = [self.continuum_file, self.narrowband_file]
        # get target wcs
        with fits.open(self.narrowband_file) as hdu:
            target_header = hdu['SCI'].header
            n_data = hdu['SCI'].data
        with fits.open(self.continuum_file) as hdu:
            c_data = hdu['SCI'].data
        
        # select smaller data to be target out shape
        if n_data.size < c_data.size:
            target_data = n_data
        else:
            target_data = c_data
            
        target_wcs = WCS(target_header)
        reprojected = reproject(list_files, target_wcs, target_data)
        return reprojected
    
    def func(x, m, c):
        """ Equation of line function to be called when m (slope) and c 
        (intersect) are obtained."""
        return m * x + c
        
    def get_line_params(self, region1, region2, c_dat, n_dat, \
                        rm_outliers=True, n_clusters=8, max_labels=2):
        """ Calculates the scale factor to be applied for continuum subtraction
        based on line fit model to n_dat vs. c_dat. Regions relate to parts of
        image (one on cloud and one off). N.B: Regions must be string slices.
        Example of region:
            region1 = "slice(3000,4000), slice(1000,1600)" """
        # stack arrays
        c_stack = np.append(c_dat[eval(region1)], c_dat[eval(region2)])
        n_stack = np.append(n_dat[eval(region1)], n_dat[eval(region2)])
        # mask continuum array for 0 and NaN entries
        cmask = (c_stack<=0) | (np.isnan(c_stack))
        nmask = n_stack==0
        c_stack = np.ma.array(c_stack, mask=cmask)
        n_stack = np.ma.array(n_stack, mask=nmask)
        
        c_flat = c_stack.flatten()
        n_flat = n_stack.flatten()
        
        if rm_outliers: # remove outliers and remake c_flat, n_flat arrays
            c_flat, n_flat = remove_outliers(c_flat, n_flat, n_clusters, max_labels)
        else:
            pass
        # fit line
        m, c = np.ma.polyfit(c_flat, n_flat, deg=1)
        return m, c, c_flat, n_flat
        
    
    def continuum_sub(self, reprojected_images, scale_factor):
        """ Subtracts continuum from narrowband data. Requires image data to be
        equal sized. Use reproject_continuum() to get equal sized images."""
        c_dat = reprojected_images[0] * scale_factor
        subc = reprojected_images[1] - c_dat
        return subc
    
# ==============================TESTING========================================
        
# directory = '/mnt/d/st_images/Carina_level3/'
# #files = glob.glob(directory + '*i2d.fits')

# # Continuum file from NIRCam (f444w)
# c_name = glob.glob(directory+'*f444w_i2d.fits')[0]

# # Narrowband file name (NIRCam f444w-f470n)
# n_name = glob.glob(directory+'*f444w-f470n*')[0]

# test = ContinuumSubtract(n_name, c_name)
# reproj = test.reproject_continuum()

# # f = test.get_data_headers()
# # c_dat, c_header = f[:2]
# # n_dat, n_header = f[2:]
# # del f # delete f object for memory consvervation

# def func(x, m, c):
#     """ Equation of line function to be called when m (slope) and c (intersect) are obtained."""
#     return m * x + c

# # setup regions for continuum subtract calculations
# region1 = "slice(3000,4000), slice(1000,1600)"
# region2 = "slice(4000,5000), slice(3000,3600)"

# m, c_, c_flat, n_flat = test.get_line_params(region1, region2, reproj[0], \
#                                             reproj[1], rm_outliers=True, max_labels=2)
# print('Line parameters:')
# print(f'Slope: {m}', f'Intersect: {c_}')

# =============================================================================
