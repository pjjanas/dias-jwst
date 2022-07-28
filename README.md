# dias-jwst
Repo for sharing code based files on analysing and carrying out photometry on images from the James Webb Space Telescope (JWST). Currently working on the star forming region of NGC 3324 in the Carina Nebula.

# Dependencies
Here are the packages used in this repo that need to be installed prior to use (specific version may not be necessary).
```
numpy=1.23.1
matplotlib=3.5.2
glob
astropy=5.1
```
***Note:*** This list is a work in progress so keep an eye on packages used to see if you need to install new ones.

# Use
Current usage is running `reproject_combine.py` with desired directory, target data etc. to get reprojected images and load functions. Then can run specific desired tasks such as making a 3-colour image using `make_rgb()` from the console. Transformed images can also be made making use of `astropy.visualization` using the `images_transform()` function or use `config_stretches.py` for custom transformations.

# Programs 
## `reproject_combine.py`
This program reads in a target image and a list of images to reproject via WCS information onto the target WCS. Need to set your own working directory using `directory` variable. `glob` is utilised here for Unix based file handling; e.g. to get all NIRCam images use `glob.glob(directory+'*nircam*')`. 

Program also includes functions to make rgb images and perform specific transforms on the reprojected images using `astropy.visualization`. If specific transformations on individual images needed, use `config_stretches.py`.

## `config_stretches.py`
Program for custom transformation of NIRCam images from JWST of the Carina Nebula region NGC 3324 (6 images). Requires a list of reprojected image arrays. To change stretch or interval for each image, edit the `transformed_NIRCam()` function and apply any changes needed and run the program to reload the function. 
