#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 17:17:56 2022

@author: Pawel Janas (Python Python 3.9.12)
"""
import numpy as np
import matplotlib.pyplot as plt

import glob

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

dir = '/mnt/d/st_images/Carina_level3/'

files = glob.glob(dir+'*nircam*')

target = files[5]

with fits.open(target) as hdu:
    tar_header = hdu[1].header
    tar_data = hdu[1].data
    
wcs = WCS(tar_header)

# 10:36:50.9373 -58:36:20.252 furthest point of outflow
# 10:36:52.3130 -58:36:19.275 closer edge
# 10 36 54.2614 -58 36 26.768 # source

dist = 2800 * u.pc # distance from earth to NGC3324

source = SkyCoord('10 36 54.2614 -58 36 26.768', distance = dist, \
                  unit=(u.hourangle, u.deg), frame='icrs')
edge1 = SkyCoord('10:36:50.9373 -58:36:20.252', distance = dist, \
                unit=(u.hourangle, u.deg))
edge2 = SkyCoord('10:36:52.3130 -58:36:19.275', distance=dist, \
                 unit=(u.hourangle, u.deg))

sep1 = source.separation_3d(edge1) # separation
sep2 = source.separation_3d(edge2)

sep1_au = sep1.to(u.au)
sep2_au = sep2.to(u.au)