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

sys.path.append('../code')
from reproject_combine import reproject

def open_fits(fname):
    """ Convenience function for reading in fits files, getting data and
    header information."""
    with fits.open(fname) as hdu:
        data = hdu['SCI'].data
        header = hdu['SCI'].header
    
    return data, header

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
    
    def get_line_params(self, region1, region2, c_dat, n_dat):
        """ Calculates the scale factor to be applied for continuum subtraction
        based on line fit model to n_dat vs. c_dat. Regions relate to parts of
        image (one on cloud and one off). N.B: Regions must be string slices.
        Example of region:
            region1 = "slice(3000,4000), slice(1000,1600)" """
        # stack arrays
        c_stack = np.append(c_dat[eval(region1)], c_dat[eval(region2)])
        n_stack = np.append(n_dat[eval(region1)], n_dat[eval(region2)])
        # mask continuum array for 0 and NaN entries
        mask = (c_stack<=0) | (np.isnan(c_stack))
        c_stack = np.ma.array(c_stack, mask=mask)
        
        c_flat = c_stack.flatten()
        n_flat = n_stack.flatten()
        
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

# f = test.get_data_headers()
# c_dat, c_header = f[:2]
# n_dat, n_header = f[2:]
# del f # delete f object for memory consvervation

# mask = (np.isnan(reproj[0])) & (reproj[0] == 0)
# # scale_factor = 2
# # subc = test.continuum_sub(reproj, 2)

# =============================================================================
