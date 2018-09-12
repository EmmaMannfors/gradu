################################################################################
#	Takes the average value of the pixels in an empty(-ish) space 
#	and subtracts this value from the entire image
################################################################################
#					THEN: 
#			Makes sure that all values are above zero
################################################################################
#		Uses the images in units of MJy/sr
################################################################################
#			Jun 27, 2018
################################################################################
#execfile('Aux.py')

#import aplpy
from astropy.io import fits
import numpy as np
import os

#from functions import distance_on_sphere_deg

########################################
#	Paths & files
########################################
dir_in = '/home/emma/gradu/data/All_data/PACS_160_Herschel/fixed_units_cropped_160/average_subtracted_160/'
dir_out = '/home/emma/gradu/data/All_data/PACS_160_Herschel/fixed_units_cropped_160/average_subtracted_160/no_negatives/'

########################################
#	Functions
########################################
def outName(filename):
	if filename.endswith('norm_1.fits'):
		outName = filename.strip('_norm_1.fits')+('_norm_pos_1.fits')
	elif filename.endswith('norm_2.fits'): 
		outName = filename.strip('_norm_2.fits')+('_norm_pos_2.fits')
	elif filename.endswith('norm_3.fits'): 
		outName = filename.strip('_norm_3.fits')+('_norm_pos_3.fits')
	elif filename.endswith('norm_4.fits'):
		outName = filename.strip('_norm_4.fits')+('_norm_pos_4.fits')
	else:
		print 'AAAA'
	return outName



########################################
#	Program
########################################
for filename in os.listdir(dir_in):
	if filename.endswith('.fits'):
		f = fits.open(dir_in+filename)[0]	
		fData = f.data
		
		# Get smallest value (sm) and takes the absolute value of it (s)
		sm = np.min(fData)
		s = abs(sm)

		#Add abs(smallest value) = s to each pixel
		for val in np.nditer(fData, op_flags=['readwrite']):
			val += s

		f.data = fData
		
		print np.min(f.data)

		out = dir_out+outName(filename)
		
		#print out

	#	f.writeto(out)



