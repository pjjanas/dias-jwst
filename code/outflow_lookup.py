#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 14:28:25 2022

@author: Pawel Janas (Python Python 3.9.12)

Description: Program for making RGB images for specific outflows.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd
import glob # for filepath handling
import sys

# Astropy:
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import LinearStretch, LogStretch, ZScaleInterval, \
    MinMaxInterval
import reproject as rpj # needed for reproject module

sys.path.append('../code')
from reproject_combine import reproject, images_transform, make_rgb

def get_coordinates(path_to_excel_file, focus):
    """ Reads in excel file into dataframe and gets RA/Dec coordinates of 
    choice source."""
    
    # Reading in excel file as DataFrame
    df = pd.read_excel(path_to_excel_file, header=18)
    dataset = df.to_dict()

    # Get coordinate data for that focus source
    names = list(dataset['Name'].values())
    index = names.index(focus)
    RA = dataset['RA (hms)'][index] # get str RA coordinates
    Dec = dataset['Dec (dms)'][index]
    coords = SkyCoord(RA+' '+Dec, unit=(u.hourangle, u.deg))
    
    return coords

def pixel_region(pixel_x, pixel_y, x_plusminus=300, y_plusminus=300, 
                 custom_bounds=None):
    """
    Function for defining a pixel region around some focus (pixel_x, pixel_y). 
    Custom bounds are of the format: 
    custom_bounds=[xminus, xplus, yminus, yplus].

    Parameters
    ----------
    pixel_x : int
        The x pixel number of the focus point.
    pixel_y : int
        The y pixel number of the focus point.
    x_plusminus : int, optional
        The area in pixels to keep to the left and right of the focus point. 
        The default is 300.
    y_plusminus : int, optional
        The area in pixels to keep to the top and bottom of the focus point. 
        The default is 300.
    custom_bounds : list, array, optional
        If custom cropping is required. List takes the form of:
        custom_bounds = [left, right, bottom, top]. The default is None.

    Returns
    -------
    region : TYPE
        DESCRIPTION.

    """
    
    if isinstance(custom_bounds, (list, np.ndarray)):
        region = f"slice({pixel_y} - {custom_bounds[2]}, {pixel_y} + {custom_bounds[3]}), \
            slice({pixel_x} - {custom_bounds[0]}, {pixel_x} + {custom_bounds[1]})"
    
    else:  
        region = f"slice({pixel_y} - {y_plusminus}, {pixel_y} + {y_plusminus}), \
        slice({pixel_x} - {x_plusminus}, {pixel_x} + {x_plusminus})"
    
    return region

def crop_images(list_images, region):
    """ Crops each image to cover region defined."""
    cropped_imgs = []
    for img in list_images:
        img = img[eval(region)]
        cropped_imgs.append(img)
    
    return cropped_imgs

# =================================TESTING====================================
# # Set outflow to focus on from table and get coordinates
# focus = 'PJ001'
# coords = get_coordinates(focus)
# 
# # Get images for RGB channels
# img_dir = '/mnt/d/st_images/Carina_level3/'
# files = glob.glob(img_dir + '*nircam*')
# target = files[5]
# 
# with fits.open(target) as hdu:
#     tar_header = hdu['SCI'].header
#     tar_data = hdu['SCI'].data
# 
# # Get target WCS and pixel values of focus outflow for image generation
# target_wcs = WCS(tar_header)
# pixel_x, pixel_y = target_wcs.world_to_pixel(coords)
# pixel_x = int(pixel_x)
# pixel_y = int(pixel_y)
# 
# y_pm = 300 # y plus minus (MUST BE int)
# x_pm = 300 # x plus minus
# 
# # Get pixel region of image
# region = pixel_region(pixel_x, pixel_y, x_pm, y_pm)
# 
# # Now we will reproject, crop and transform images (in that order)
# reproj = reproject(files, target_wcs, tar_data)
# cropped_imgs = crop_images(reproj, region)
# transformed = images_transform(cropped_imgs, stretch=LogStretch(), \
#                                interval=MinMaxInterval())
# 
# =============================================================================



