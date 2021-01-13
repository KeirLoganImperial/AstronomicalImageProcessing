# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 17:53:40 2021

@author: logan
"""

from astropy.io import fits
import os
import numpy as np

read_path = 'mosaic.fits'
write_path = 'image_cropped_refined.fits'

#opens image data
hdul = fits.open(read_path)

#cropping out edge regions with less subexposures
hdul[0].data = hdul[0].data[470:4400,200:2340]
hdulist = hdul[0].data 


def circular_boolean(centre, radius):
    """
    Function that defines boolean array with 'True' in specified circular region.
    
    centre (tuple) - defines x and y value of the centre of the circular array.
    radius (float or int) - defines size of circular truth array around centre.
    
    returns: circle - boolean array with 'True' values within a circle of specified radius and centre.
    """
    #define size of image
    y = len(hdulist)
    x = len(hdulist[0])
    Y, X = np.ogrid[:y, :x]

    #calculate distance from each point to centre
    disp_from_centre = np.sqrt((X - centre[0])**2 + (Y-centre[1])**2)
    
    #distances less than radius given value 'True'
    circle = disp_from_centre <= radius
    return circle


def rectangular_boolean(centre, dimensions):
    """
    Function that defines boolean array with 'True' in specified rectanglar region.
    
    centre (tuple) - defines x and y value of the centre of the rectangular array.
    dimensions (tuple) - defines x and y dimensions for the rectangular truth array.
    
    returns: rect - boolean array with 'True' values within a circle of specified radius and centre.
    """
    #defines boolean array with 'true' in masked rectangle
    
    #uniform 'false' array
    rect = np.zeros((len(hdulist), len(hdulist[0])), dtype=bool)
    
    #rectangle's coordinates to be masked
    x0=int(centre[0]-dimensions[0]/2)
    x1=int(centre[0]+dimensions[0]/2)
    y0=int(centre[1]-dimensions[1]/2)
    y1=int(centre[1]+dimensions[1]/2)
    
    #set true in rectangle
    rect[y0:y1, x0:x1] = True
    return rect


# Implementing the masking for identified problematic objects

circ1 = circular_boolean(centre=(1230,2746), radius=230)
rect1 = rectangular_boolean(centre=(1238,1981), dimensions=(23,3962))

circ2 = circular_boolean(centre=(704,1820), radius=40)
rect2 = rectangular_boolean(centre=(704,1820), dimensions=(10,140))

circ3 = circular_boolean(centre=(772,2304), radius=45)
rect3 = rectangular_boolean(centre=(772,2304), dimensions=(10,140))

circ4 = circular_boolean(centre=(577,2849), radius=45)
rect4 = rectangular_boolean(centre=(577,2849), dimensions=(10,231))

circ5 = circular_boolean(centre=(1933,3288), radius=30)
rect5 = rectangular_boolean(centre=(1933,3288), dimensions=(10,140))

circ6 = circular_boolean(centre=(1889,957), radius=30)
circ7 = circular_boolean(centre=(1932,1839), radius=30)
circ8 = circular_boolean(centre=(1257,3561), radius=30)
circ9 = circular_boolean(centre=(357,3626), radius=30)
circ10 = circular_boolean(centre=(2060,2827), radius=40)

mask = np.zeros((len(hdulist), len(hdulist[0])), dtype=bool)

circs = circ1 | circ2 | circ3 | circ4 | circ5 | circ6 | circ7 | circ8 | circ9 | circ10
rects = rect1 | rect2 | rect3 | rect4 | rect5

#setting values of the image in masked regions to zero counts
mask = circs | rects
hdulist[mask] = 0

if os.path.exists(write_path):
    os.remove(write_path)
hdul.writeto(write_path)
hdul.close()