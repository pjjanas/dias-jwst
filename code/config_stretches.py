#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 19:48:18 2022

@author: Pawel Janas (Python Python 3.9.12)

Description: Config program to compute stretches and intervals for NIRcam images.

Normalises all images. Yet to implement non-normalised stretches/intervals.
To be passed to make_rgb in reproject_combine.py.
"""

# excess import if changes needed
from astropy.visualization import LinearStretch, LogStretch, ZScaleInterval, \
    MinMaxInterval,  ManualInterval
    
def transformed_NIRCam(reprojected):
    transformed_imgs = []
    
    # 1st Image (f090w) - blue
    interval = ZScaleInterval(contrast=0.1)
    stretch = LinearStretch()
    transform = stretch + interval
    t = transform(reprojected[0])
    transformed_imgs.append(t)
    
    # 2nd Image (187n) - blue - GOOD
    interval = ZScaleInterval(contrast=0.2)
    stretch = LinearStretch()
    transform = stretch + interval
    t = transform(reprojected[1])
    transformed_imgs.append(t)
    
    # 3rd Image (f200w) - blue
    interval = ZScaleInterval()
    stretch = LinearStretch()
    transform = stretch + interval
    t = transform(reprojected[2])
    transformed_imgs.append(t)
    
    # 4th Image (f335m) - green
    interval = ZScaleInterval(contrast=0.1)
    stretch = LinearStretch()
    transform = stretch + interval
    t = transform(reprojected[3])
    transformed_imgs.append(t)
    
    # 5th Image (f400w) - red
    interval = ZScaleInterval(contrast=0.15)
    stretch = LinearStretch()
    transform = stretch + interval
    t = transform(reprojected[4])
    transformed_imgs.append(t)
    
    # 6th Image (f470n) - red
    interval = ZScaleInterval(contrast=0.3)
    stretch = LinearStretch()
    trans = stretch + interval
    t = trans(reprojected[5])
    transformed_imgs.append(t)
    
    return transformed_imgs

    
    
    


