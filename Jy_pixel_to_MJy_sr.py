########################################################################
# Change 160 um map to MJy/sr from Jy/pixel
########################################################################

from astropy.io import fits
import os
import math
import numpy as np
################################
# Paths
################################
in_70 = '/home/emma/gradu/data/All_data/PACS_70_Herschel/cropped_70/'
in_100 = '/home/emma/gradu/data/All_data/PACS_100_Herschel/cropped_100/'
in_160 = '/home/emma/gradu/data/All_data/PACS_160_Herschel/cropped_160/'

out_70 = '/home/emma/gradu/data/All_data/PACS_70_Herschel/MJy_cropped_70/'
out_100 = '/home/emma/gradu/data/All_data/PACS_100_Herschel/MJy_cropped_100/'
out_160 = '/home/emma/gradu/data/All_data/PACS_160_Herschel/fixed_units_cropped_160/'

testFile = in_160+'G195.74-02.29_160_L2.5_J_crop.fits'

################################
# Convert Jy/pixel to MJy/sr
################################
def convert(dir_in,dir_out):
	for filename in os.listdir(dir_in):
		if filename.endswith('.fits'):
	#	if filename == testFile:
			f = fits.open(dir_in+filename)[0]
			fData = f.data
	
			################################
			# Doing the calculations for 
			# each value in array
			################################
			x = (1.0/(3.2)**2)		# pixels / arcsec^2
			y = (1.0/(4.8481368e-6))**2	# arcsec^2 / rad^2
			z = 10.0**(-6.0)			# MJy / Jy
	
	
			for val in np.nditer(fData, op_flags=['readwrite']):
				val *= (x*y*z)
			
			################################
			# Write to output file
			# change units to MJy/sr in hdr
			################################
			f.data = fData
			if filename.endswith('_crop.fits'):
				outName = dir_out+filename.strip('_crop.fits')+'_units_1.fits'
			elif filename.endswith('_crop_2.fits'):
				outName = dir_out+filename.strip('_crop_2.fits')+'_units_2.fits'
			elif filename.endswith('_crop_3.fits'):
				outName = dir_out+filename.strip('_crop_3.fits')+'_units_3.fits'
			elif filename.endswith('_crop_4.fits'):
				outName = dir_out+filename.strip('_crop_4.fits')+'_units_4.fits'
			else:
				outName = dir_out+'whoops'
				print filename
	
			f.header['BUNIT'] = 'MJy/sr'
			f.writeto(outName)
	
	
			continue
		else:
			continue





def check(out_dir):
	n = raw_input('Check BUNIT?\n')
	if n == 'y':
		for filename in os.listdir(out_dir):
			if filename.endswith('.fits'):
				f = fits.open(out_dir+filename)[0]
				hdr = f.header['BUNIT']
				if hdr != 'MJy/sr':
					print hdr




band = raw_input('What wavelength to fix? (70, 100, 160)\n')
if band == '70':
	convert(in_70,out_70)
	check(out_70)
elif band == '100':
	convert(in_100,out_100)
	check(out_100)
elif band == '160':
	convert(in_160,out_160)
	check(out_160)
else:
	print 'meowmeowmeow'






