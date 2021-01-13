# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 11:00:35 2021

@author: logan
"""

import numpy as np
from astropy.io import fits

masks = np.load('apertures.npy')
objects = np.loadtxt('object_counts.txt')
locations = np.loadtxt('locations.txt')

hdul = fits.open('output.fits')
hdulist = hdul[0].data
hdul.close()

size=[len(hdulist[0]),len(hdulist)]

object_array=np.zeros((size[0],size[1]))
for ind, obj in enumerate(objects):
    loc=locations[ind]
    object_array[int(loc[0]),int(loc[1])] = obj

hdulist_empty = np.zeros((size[1],size[0]), dtype=bool)

def circular_boolean(centre, radius, Hdulist):
    #defines boolean array with 'true' in masked circle
    
    #define size of image
    y = len(Hdulist)
    x = len(Hdulist[0])
    Y, X = np.ogrid[:y, :x]

    #calculate distance from each point to centre
    disp_from_centre = np.sqrt((X - centre[0])**2 + (Y-centre[1])**2)
    
    #distances less than radius made 'true'
    circle = disp_from_centre <= radius
    return circle

extended_locations = []
objects_with_extendeds = []
for i in range(7,18):
    for j in range(0, len(masks[0])):
        for k in range(0, len(masks)):
            mask = circular_boolean((j,k), radius=i, Hdulist=hdulist_empty)
            extended_overlap = np.logical_and(mask, masks)
            mask_overlap = np.logical_xor(mask, extended_overlap).flatten()
            
            if True not in mask_overlap:
                extended_object = circular_boolean((j,k), radius=i-6, Hdulist=hdulist_empty)
                obj=[]
                for m in range(0, len(masks[0])):
                    for n in range(0, len(masks)):
                        if extended_object[n,m]==True:
                            obj.append(object_array[m,n])
                objects_with_extendeds.append(sum(obj))
                extended_locations.append((j,k))
                print(j,k)
            else:
                objects_with_extendeds.append(object_array[j,k])
                extended_locations.append(locations[k,j])
            masks[mask]=False
                            
                            