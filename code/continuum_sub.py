#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 17:46:56 2022

@author: Pawel Janas (Python Python 3.9.12)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import glob # for filepath handling
import sys

# Astropy:
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import LinearStretch, LogStretch, ZScaleInterval, \
    MinMaxInterval, make_lupton_rgb
import reproject as rpj

sys.path.append('../code')
from reproject_combine import reproject, make_rgb

def open_fits(fname):
    """ Convenience function for reading in fits files, getting data and
    header information."""
    with fits.open(fname) as hdu:
        data = hdu['SCI'].data
        header = hdu['SCI'].header
    
    return data, header

directory = '/mnt/d/st_images/Carina_level3/'
files = glob.glob(directory + '*i2d.fits')

# open continuum file from NIRCam (f444w-f470n)
c_name = glob.glob(directory+'*f444w_i2d.fits')[0]
cont_data, cont_header = open_fits(c_name)
# c_dat_flat = cont_data.flatten() # flatten for plotting

# open narrowband file
n_name = glob.glob(directory+'*f444w-f470n*')[0]
narrow_dat, narrow_head = open_fits(n_name)
# n_dat_flat = narrow_dat.flatten()

list_files = [c_name, n_name]
target_wcs = WCS(cont_header)
reproj = reproject(list_files, target_wcs, cont_data)

c_dat_flat = reproj[0].flatten()
n_dat_flat = reproj[1].flatten()
