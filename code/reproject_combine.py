#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 13:40:24 2022

@author: Pawel Janas (Python 3.9.12)

Description: Combine images using their WCS information and scale them to fit
reference image set. Can make rgb image with custom stacking for individual r, g, 
b channels. 
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import glob # for filepath handling

# Astropy:
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import LinearStretch, LogStretch, ZScaleInterval, \
    MinMaxInterval, make_lupton_rgb
import reproject as rpj

# Functions
def plot_image(image, norm=False, axis=True): # convenience function
    plt.figure(figsize=(12,8))
    plt.imshow(image, origin='lower')
    
    if norm:
        plt.imshow(image, origin='lower', norm=norm)
        
    if axis==False:
        plt.axis('off')

def reproject(list_fits_files, target_wcs, target_data):
    """
    Reprojects and resizes list of fits files onto some target image with WCS 
    information.

    Parameters
    ----------
    list_fits_files : list
        List of fits files to reproject onto target WCS.
    target_data : array
        2D array of target image data.
        
    Returns
    -------
    reprojected : list
        List of reprojected and resized image arrays.

    """
    reprojected = [] # init list for reprojected images
    for file in list_fits_files:

        with fits.open(file) as hdu:
            new_image = rpj.reproject_interp(hdu[1].copy(), target_wcs, shape_out=target_data.shape)[0]
            reprojected.append(new_image)

    return reprojected


def images_transform(list_images, stretch=LinearStretch(), \
                   interval=ZScaleInterval()):
    """
    Function for transforming a list of images to some desired stretch/interval.
    If a custom stretch/interval for each image is needed. Use config_stretches.py.

    Parameters
    ----------
    list_images : list
        List of images to be transformed.
    stretch : class, optional
        The stretch to be applied to the images. Defaults to LinearStretch().
    interval : class, optional
        The interval to be applied on the images. Defaults to ZScaleInterval().

    Returns
    -------
    transformed_list: list
        List of transformed images.

    """
    transform = stretch + interval
    transformed_list = []
    for img in list_images:
        new_img = transform(img)
        transformed_list.append(new_img)
    
    return transformed_list

def make_rgb(list_images, stretch=5, Q=8, minimum=0, custom_stack=None):
    """
    Function to make rgb images using astropy.visualisation.make_lupton_rgb.

    Parameters
    ----------
    list_images : list
        List of reprojected image arrays obtained using reproject().
    stretch : float, optional
        The linear stretch of the image to use with make_lupton_rgb(). 
        The default is 5.
    Q : float, optional
        The asinh softening parameter for make_lupton_rgb(). The default is 8.
        Recommend using Q = 0.
    minimum : float, optional
        Intensity to be mapped to black (can use array or list for R, G, B). 
        The default is 0.
    custom_stack : list, str, optional
        Array or list of size len(list_images) to denote stacking order of  each 
        R, G, B bands. 
        Example: Want to stack the first 2 images for blue, 1 image for green 
        and 3 images for red, use custom_stack = [0, 1, 2, 2, 3, 5]. If uneven 
        number of input images and custom_stack=None, then first split is
        largest.
        The default is None.

    Returns
    -------
    rgb_img : ndarray
        RGB (integer, 8-bits per channel) color image as an NxNx3 numpy array.

    """
    # First bundle images into r, g, b channels (must go from shortest to longest wavelength)
    if isinstance(custom_stack, (list, np.ndarray)):
        custom_stack = np.array(custom_stack) # make sure type=np array
        list_images = np.array(list_images)
        list_images = np.split(list_images, len(list_images)) # split into equal chunks
        
        if custom_stack[0]==custom_stack[1]:
            b_stack = np.sum(list_images[custom_stack[0]], axis=0)
        else: 
            b_stack = np.sum(list_images[custom_stack[0]:custom_stack[1]], axis=0)[0,:,:]
        
        if custom_stack[2]==custom_stack[3]:
            g_stack = np.sum(list_images[custom_stack[2]], axis=0)
        else:
            g_stack = np.sum(list_images[custom_stack[2]:custom_stack[3]], axis=0)[0,:,:]
            
        if custom_stack[4]==custom_stack[5]:
            r_stack = np.sum(list_images[custom_stack[4]], axis=0)
        else:
            r_stack = np.sum(list_images[custom_stack[4]:custom_stack[5]], axis=0)[0,:,:]
                
    else: 
        list_images = np.array(list_images)
        list_images = np.array_split(list_images, 3) # split into 3 channels, if uneven; first split is largest
        b_stack = np.sum(list_images[0], axis=0)
        g_stack = np.sum(list_images[1], axis=0)
        r_stack = np.sum(list_images[2], axis=0)
        
    rgb_img = make_lupton_rgb(r_stack, g_stack, b_stack, stretch=stretch, Q=Q, \
                              minimum=minimum)
    return rgb_img
    
