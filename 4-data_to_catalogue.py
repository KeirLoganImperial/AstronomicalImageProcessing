# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 18:01:23 2021

@author: logan
"""

import numpy as np
from astropy.io import ascii
from astropy.table import QTable
import os

filenames = ['object_locations', 'object_counts', 'object_error','background','background_error','apparent','apparent_error']

for file in filenames:
    vars()[file] = np.loadtxt(file+'.txt')
    
xloc=object_locations[:,0]
yloc=object_locations[:,1]

data = QTable([xloc, yloc, object_counts, object_error, background, background_error, apparent, apparent_error], names=['x', 'y', 'count', 'count_error', 'background', 'background_error', 'magnitude', 'magnitude_error'])
ascii.write(data, 'object_catalogue.csv', overwrite=True)

for file in filenames:
    os.remove(file+'.txt')



