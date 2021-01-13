# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 22:16:20 2021

@author: logan
"""

### NOTE!! - This script takes around ONE HOUR to perform.

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from datetime import datetime
import os

read_path = 'image_cropped_refined.fits'
write_path = 'output.fits'

#read data and flatten to 1d numpy array
hdul = fits.open(read_path)
hdulist = hdul[0].data
hdulist_flat = hdulist.flatten()
hdul.close()

#plot histogram of pixel counts in image
plt.figure(1)
histmin=3350
histmax=3600
plt.hist(hdulist_flat, range=(histmin,histmax), bins=(histmax-histmin))
plt.savefig('distribution.png', dpi=1000)
plt.title('Pixel value distribution - cropped and masked')
plt.ylabel('Occurence')
plt.xlabel('Pixel value')
plt.show()


def circular_boolean(centre, radius, Hdulist=hdulist):
    """
    Function that defines boolean array with 'True' in specified circular region.
    
    parameters:
        centre (tuple) - defines x and y value of the centre of the circular array.
        radius (float or int) - defines size of circular truth array around centre.
    
    returns:
        circle - boolean array with 'True' values within a circle of specified radius and centre.
    """
    
    #define size of image
    y = len(Hdulist)
    x = len(Hdulist[0])
    Y, X = np.ogrid[:y, :x]

    #calculate distance from each point to centre
    disp_from_centre = np.sqrt((X - centre[0])**2 + (Y-centre[1])**2)
    
    #distances less than radius made 'true'
    circle = disp_from_centre <= radius
    return circle


def aperture_photometry(Hdulist, div):
    """
    Function to perform aperture photometry on objects in image data 'hdulist'.
    
    parameters:
        Hdulist (2D array) - two-dimensional array of pixel counts for analysis.
        div (tuple) - number of x and y subdivisions to split up image for analysis.
    
    returns:
        objects - list of the pixel sums of the measured objects.
        locations - list of tuples of the detected objects' x and y positions in the image.
        obj_err - list of errors in the total object counts.
        backgrounds - list of calculated background count for each detected object.
        backg_errs - list of errors in these calculated backgrounds.
        masks - boolean array of the size of the image with True values where objects were detected.
    """
    
    print('0 % complete - Elapsed time: 0:00:00.000000 - Estimated time left: N/A')
    
    #initialise arrays and values
    StartTime = datetime.now()
    progress = 0
    objects=[]
    obj_err=[]
    locations=[]
    backgrounds=[]
    backg_errs = []
    masks = np.zeros((len(Hdulist),len(Hdulist[0])), dtype=bool)
    
    #calculate x and y dimensions of subdivisions
    xdiv=int(len(Hdulist[0])/div[0])
    ydiv=int(len(Hdulist)/div[1])  
    
    #set scan parameters for background estimation
    rad=100
    max_thresh = 3650
    min_thresh = 3450
    
    #iterate over subdivisions of x and y indice ranges xi to (xi + xdiv) and yi to (yi + ydiv) respectively
    xi=0
    for i in range(0,div[0]):
        yi=0
        for j in range(0,div[1]):    
            
            #crop image data to within subdivision
            Hdulist_cropped=Hdulist[yi:yi+ydiv, xi:xi+xdiv]
            
            #flat array enables speedier finding of max value and indices
            flat = Hdulist_cropped.flatten()
            
            #find x and y values of maximum pixel count in subdivision
            xmax = xi + np.argmax(flat)%len(Hdulist_cropped[0])
            ymax = yi + int(np.argmax(flat)/len(Hdulist_cropped[0]) - xmax/len(Hdulist_cropped[0]))
            
            #(rather crude) method of calculating background and its error within circular region around maximum value
            reference = Hdulist[circular_boolean(centre=(xmax,ymax), radius=rad, Hdulist=Hdulist)].flatten()
            statslist=[]
            for i in reference:
                if i<max_thresh and i!=0:
                    statslist.append(i)
            background = np.mean(statslist)  
            st_dev = np.std(statslist)
            backg_frac_error = st_dev/len(statslist)/background
            
            #define threshold for object to be detected based on local background count but within defined threshold range
            threshold = min(max(background + 2*st_dev, min_thresh),max_thresh)
            
            
            #iterate over subdivision data until no objects left to detect above their individual thresholds
            while max(flat)>threshold:
                
                
                #defines aperture of radius 6 pixels
                aperture = circular_boolean(centre=(xmax,ymax), radius=6, Hdulist=Hdulist)
                
                #records geometry of aperture in 2D truth array
                masks[aperture]=True
                
                #calculation of total pixel count for detected object, minus the background count, plus errors
                count=0
                for i in Hdulist[aperture]:
                    if i-background>0:
                        count+= (i-background)
                combined_frac_err = 0.0985 + backg_frac_error
                
                #save data to lists
                objects.append(count)
                obj_err.append(combined_frac_err*count)
                backgrounds.append(background)
                backg_errs.append(backg_frac_error*background)
                
                #records location of detected object
                locations.append((int(xmax),int(ymax)))
                
                #wipes image data for detected object
                Hdulist[aperture]=0
                
                #recalculates xmax and ymax for next object
                flat = Hdulist_cropped.flatten()
                xmax = xi + np.argmax(flat)%len(Hdulist_cropped[0])
                ymax = yi + int(np.argmax(flat)/len(Hdulist_cropped[0]) - xmax/len(Hdulist_cropped[0]))
                
                #recalculates background, threshold and errors for next object
                reference = Hdulist[circular_boolean(centre=(xmax,ymax), radius=rad, Hdulist=Hdulist)].flatten()
                statslist=[]
                for i in reference:
                    if i<max_thresh and i!=0:
                        statslist.append(i) 
                background = np.mean(statslist)
                st_dev = np.std(statslist)
                backg_frac_error = st_dev/len(statslist)/background
                threshold = min(max(background + 2*st_dev, min_thresh),max_thresh)

            
            #calculates progress and estimates time left to complete scan
            progress += 1/(div[0]*div[1])
            elapsed_time = datetime.now() - StartTime
            time_to_go = elapsed_time/progress - elapsed_time
            print(int(progress*100), '% complete - Elapsed time:', elapsed_time, '- Estimated time left:', time_to_go)
            
            #periodically saves masked data to fits file for viewing (can use to see progress etc.)
            hdul.writeto(write_path, overwrite=True)
            hdul.close()
            
        #updates x and y indices for next subdivision
            yi += ydiv
        xi += xdiv
               
    return objects, locations, obj_err, backgrounds, backg_errs, masks

#performs aperture photometry on input image hdulist using 5 subdivisions in x axis and 10 subdivisions in y axis
objects, locations, obj_err, backgrounds, backg_errs, masks = aperture_photometry(hdulist, div=(5,10))

#saves output data for analysis and cataloguing
np.savetxt('object_counts.txt', objects)
np.savetxt('object_locations.txt', locations)
np.savetxt('object_error.txt', obj_err)
np.savetxt('background.txt', backgrounds)
np.savetxt('background_error.txt', backg_errs)
np.save('apertures.npy', masks)


