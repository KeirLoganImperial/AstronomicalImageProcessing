# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 11:23:38 2020

@author: logan
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

read_path = 'output.fits'
hdul = fits.open(read_path)
objects = np.loadtxt('object_counts.txt')
locations = np.loadtxt('object_locations.txt')
object_err = np.loadtxt('object_error.txt')
background = np.loadtxt('background.txt')
backg_error = np.loadtxt('background_error.txt')
masks = np.load('apertures.npy')

#extracts instrumental zero point from dataset and calculates calibrated magnitudes from the object counts
zero_point = hdul[0].header['MAGZPT']
zero_point_err = hdul[0].header['MAGZRR']
calibrated = zero_point + -2.5*np.log10(objects)
calibrated_err = 0.434*0.0985*2.5 + zero_point_err   #error of 9.85% calulcated externally; factor of 0.434 comes from error propagated through log10 function
err_array = calibrated_err + np.zeros(len(calibrated))

#defines evenly spaced array of magnitude values and iterates over these and the object magnitudes to calculate number count data at these intervals
#calculates lower and upper error for each magnitude value by considering number count for the extremities of the magnitude errors
mags = np.linspace(min(calibrated),max(calibrated),num=50)
number=[]
low_errs=[]
high_errs=[]
for i in mags:
    count=0
    low_count=0
    high_count=0
    xerror=0.0985*i     
    for j in calibrated:
        if j<i:
            count+=1
        if j<i-calibrated_err:
            low_count+=1
        if j<i+calibrated_err:
            high_count+=1
    number.append(count)
    low_errs.append(count-low_count)
    high_errs.append(high_count-count)

#convert lists to numpy arrays
number = np.array(number)

#finds best fit for number count plot between magnitudes of 10 and 18, calculates slope and slope error
mags_linear=[]
number_linear=[]
for i,m in enumerate(mags):
    if m>10 and m<18:
        mags_linear.append(m)
        number_linear.append(number[i]) 
m_line = np.linspace(10,18.5)
fit, cov = np.polyfit(mags_linear,np.log10(number_linear), 1, cov=True)
slope = fit[0]
slope_err = np.sqrt(cov[0][0])

#plots data with logarithmic y axis, with error bars, and saves image of plot
plt.plot(mags,number,'x')
plt.plot(m_line,10**(fit[0]*m_line + fit[1]),'-')
plt.xlim((min(calibrated)-0.5,20))
plt.ylim(1,20000)
plt.yscale('log')
plt.errorbar(mags,number, yerr=np.array([low_errs,high_errs]), fmt='none',capsize=3)#xyerr=0.0985*mags)
plt.xlabel('Magnitude', fontsize=13)
plt.ylabel('Number', fontsize=13)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.grid()
plt.savefig('number_count_plot.png', dpi=1000)
plt.show()

#saves arrays of apparent magnitudes & their errors for cataloguing
np.savetxt('apparent.txt', calibrated)
np.savetxt('apparent_error.txt', err_array)

hdul.close()
