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
def plot_image(image, norm=False, axis=True, dpi=300): # convenience function
    plt.figure(figsize=(10,10), dpi=dpi)
    plt.imshow(image, origin='lower')
    
    if norm:
        plt.imshow(image, origin='lower', norm=norm)
        
    if axis==False:
        plt.axis('off')

def reproject(list_fits_files, target_data, write_new=False, dir_=None):
    """
    Reprojects and resizes list of fits files onto some target image with WCS 
    information.

    Parameters
    ----------
    list_fits_files : list
        List of fits files to reproject onto target WCS.
    target_data : array
        2D array of target image data.
    write_new : bool, optional
        If True, write new image data to new .fits files. Currently names files 
        as reprojectedN, where N is file number unless dir_ is passed. 
        The default is False.
    dir_ : str, optional
        Name of directory to put new .fits files if write_new=True. 
        The default is None.

    Returns
    -------
    reprojected : list
        List of reprojected and resized image arrays.

    """
    reprojected = [] # init list for reprojected images
    for file in list_fits_files:
        img_dat = fits.open(file)[1] # get image data
        new_image = rpj.reproject_interp(img_dat, target_wcs, \
                                         shape_out=target_data.shape)[0]
        reprojected.append(new_image)
    return reprojected
    
    if write_new: # if True write data to new .fits files
        if isinstance(dir_, str):
            for i in range(len(reprojected)):
                fits.writeto(dir_+f'reprojected{i}.fits', reprojected[i], overwrite=True)
        else:
            for i in range(len(reprojected)): # stay in same directory
                fits.writeto(f'reprojected{i}.fits', reprojected[i], overwrite=True)

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
    

# open target image
directory = '/mnt/d/st_images/Carina_level3/' # set your own dir here
target_file = directory + 'jw02731-o001_t017_nircam_clear-f335m_i2d.fits'
target = fits.open(target_file)

# get data and science headers
tar_sci_data = fits.getdata(target_file, 'SCI')
tar_sci_header = fits.getheader(target_file, 1)

list_fits_files = glob.glob(directory+'*nircam*') # get all NIRCam images in list
target_wcs = WCS(tar_sci_header) # get target WCS info

# resize all images to have same dimensions as target image [NOT NEEDED]
# resized_images = []
# for file in list_fits_files:
#     img = fits.getdata(file, 'SCI')
#     new_image = resize(img, (tar_sci_data.shape[0], tar_sci_data.shape[1]), \
#                        anti_aliasing=True)
#     resized_images.append(new_image)
    
reprojected = reproject(list_fits_files, tar_sci_data) # get reprojected image arrays

# ================================TESTING======================================
# 
# # Trying with only two images
# # f1 = fits.open(list_fits_files[0])
# # test_image = rpj.reproject_interp(f1[1], target[1].header)[0]
# 
# # # compare img and new_image
# # plt.figure(figsize=(10,10), dpi=300)
# # ax1 = plt.subplot(1,2,1, projection=WCS(fits.getheader(list_fits_files[2],1)))
# # ax1.imshow(img, origin='lower', norm=LogNorm())
# # ax1.coords['ra'].set_axislabel('Right Ascension')
# # ax1.coords['dec'].set_axislabel('Declination')
# # ax1.set_title('Original Image')
# 
# # ax2 = plt.subplot(1,2,2, projection=target_wcs)
# # ax2.imshow(new_image, origin='lower', norm=LogNorm())
# # ax2.coords['ra'].set_axislabel('Right Ascension')
# # ax2.coords['dec'].set_axislabel('Declination')
# # ax2.set_title('Reprojected Image')
# 
# =============================================================================
