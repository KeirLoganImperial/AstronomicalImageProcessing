The scripts in this folder were written and executed within Spyder IDE version 4.1.4 running Python 3.8.3 64-bit within Windows 10.

The scripts should be run in the order of the number in their file name:

'1-image_preparation.py' - prepares image and data for object detection by removing unwanted objects and cropping
'2-source_detection.py' - performs aperture photometry on the image
'3-analysis.py' - analyses the output data and plots the number count
'4-data_to_catalogue.py' - saves object data into an ASCII CSV file

The 5th script - 'DOES_NOT_WORK-extended_objects_conversion.py' - is incomplete and does not form part of the object analysis.
This script was intended to detect extended objects larger than the aperture size, by checking for circular regions (larger than the aperture size) in the 'output.fits' image which contain only detected objects - i.e. a potential extended object.
This was included to be noted for further improvement.

The scripts utilise the following libraries:

numpy
matplotlib
astropy
os
datetime



Input (required to be in the same folder as the 4 scripts):

'mosaics.fits'  - a FITS file containing the raw astronomical image to be analysed



Output (after running all four scripts):

'object_catalogue.csv' - CSV file containing the data collected for objects detected in the 'mosaics.fits' image
Object data stored in this file:
x,y - x and y positions
count - summated pixel count
count_error - error in pixel count
background - calculated local background count
background_error - error in background count
magnitude - calculated apparent magnitude
magnitude_error - error in apparent magnitude value

'image_cropped_refined.fits' - image after removing statistically problematic objects and cropping outer regions with less subexposures

'output.fits' - the "image_cropped_refined.fits" image with detected objects removed

'number_count_plot.png' - high resolution image of the produced number count plot

'distribution.png' - high resolution image of a histogram of the pixel counts in 'image_cropped_refined.fits'

'apertures.pny' - a 2D numpy boolean array of the same dimensions as "output.fits" containing 'True' in the x-y positions within the aperture of a detected object.
This was intended to be used in the "extended_objects_conversion.py" script.