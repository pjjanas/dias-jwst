# dias-jwst
Repo for sharing code based files on analysing and carrying out photometry on images from the James Webb Space Telescope (JWST). Currently working on the star forming region of NGC 3324 in the Carina Nebula. Looking for protostellar jets in the surrounding regions, previously masked by molecular clouds. JWST's NIRCam instrument allows us to look through the clouds, uncovering a myriad of young stellar objects (YSOs).

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

<img src="https://user-images.githubusercontent.com/81090178/185171680-a0304640-aa18-4a93-806c-514d284fd013.png" width=100%, height=100%>

# Dependencies
Here are the packages used in this repo that need to be installed prior to use (specific version may not be necessary).
```
numpy=1.23.1
matplotlib=3.5.2
glob
astropy=5.1
scikit-learn=1.1.2
pandas
ipympl # for jupyter notebook plot interaction
```
***Note:*** This list is a work in progress so keep an eye on packages used to see if you need to install new ones.

# Use
Current usage is running `reproject_combine.py` with desired directory, target data etc. to get reprojected images and load functions. Then can run specific desired tasks such as making a 3-colour image using `make_rgb()` from the console. Transformed images can also be made making use of `astropy.visualization` using the `images_transform()` function or use `config_stretches.py` for custom transformations. See the [Example RGB Image Notebook](https://github.com/pjjanas/dias-jwst/blob/main/notebooks/Make_RGB_image_example.ipynb) for implementation.

There is now an implementation for continuum subtraction of images using `continuum_sub.py` in which there is a class with the necessary functions needed for continuum subtraction. See the [Continuum Subtraction Notebook](https://github.com/pjjanas/dias-jwst/blob/main/notebooks/Continuum_subtract.ipynb) for usage and different plots for visualisation of the key steps required.

Examples of images produced using `make_rgb()`:

<img src="https://user-images.githubusercontent.com/81090178/182904964-8670ae3f-53c4-4b83-ba3c-e813f129b6ac.png" width=75% height=75%>
<img src="https://user-images.githubusercontent.com/81090178/182905953-e6aa4cfb-12e1-44f6-99b8-42af8bf85ca2.jpg" width=75% height=75%>
<img src="https://user-images.githubusercontent.com/81090178/182910159-41b22a51-600b-4446-807f-7a0449172516.png" width=75% height=75%>



# Programs 
## `reproject_combine.py`
This program reads in a target image and a list of images to reproject via WCS information onto the target WCS. Need to set your own working directory using `directory` variable. `glob` is utilised here for Unix based file handling; e.g. to get all NIRCam images use `glob.glob(directory+'*nircam*')`. 

Program also includes functions to make rgb images and perform specific transforms on the reprojected images using `astropy.visualization`. If specific transformations on individual images needed, use `config_stretches.py`.

## `config_stretches.py`
Program for custom transformation of NIRCam images from JWST of the Carina Nebula region NGC 3324 (6 images). Requires a list of reprojected image arrays. To change stretch or interval for each image, edit the `transformed_NIRCam()` function and apply any changes needed and run the program to reload the function.

## `continuum_sub.py`
Program made for subtracting continuum from a desired image. In our case, we're using the NIRCam F444W-F470N as our narrowband image and F444W image as the adjacent continuum. The class `ContinuumSubtract()` features all functions necessary to carry out subtraction. `scikit-learn` was utilised for outlier removal for determining the `scale_factor` to apply to the continuum image before subtraction.
