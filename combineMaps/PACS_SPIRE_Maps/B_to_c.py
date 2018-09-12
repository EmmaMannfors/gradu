#################################################################################
#	Convolve 350, 250, and 160 um's calculate_n code's NH2 maps to 
#				500 um FWHM
#		 		I.E. B to C
#################################################################################
#				DO NOT REPROJECT
#################################################################################
import montage_wrapper as montage
from astropy.io import fits
import os
import numpy as np


from astropy.utils.data import get_pkg_data_filename
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy.convolution import Gaussian1DKernel
########################
# Paths
########################
B_in = '/home/emma/gradu/data/All_data/n_tau_all/350_160_n_tau/'
C_out = '/home/emma/gradu/data/All_data/n_tau_all/350_160_n_tau/convolved_to_500/'

B_C_out = '/home/emma/gradu/data/All_data/n_tau_all/350_160_n_tau/B_minus_C/'
########################
# Functions
########################
def sigma(img):
	sigma500 = 35.4
	sigmaU = 24.2		#350 um resolution
	minus = sigma500**2 - sigmaU**2
	sq = math.sqrt(minus)
	
	CDELT1_deg = img.header['CDELT1']
	CDELT1 = abs(CDELT1_deg*3600.0)
	print 'C= ',CDELT1
	new_FWHM = sq / CDELT1
#	print 'new FWHM= ', new_FWHM
	sigma = new_FWHM/(math.sqrt(CDELT1*8*math.log(2)))
#	print sigma
	return sigma



########################
# Code
########################
for filename in os.listdir(B_in):
	if filename.endswith('NH2.fits'):
		current = get_pkg_data_filename(B_in+str(filename))
		img_2 = fits.open(current)[0]
		img = img_2.data
		sKernel2 = sigma(img_2)
		#sKernel2 = (0.981339795226)
		# Create kernel and convolve
		kernel = Gaussian2DKernel(sKernel2)
		conv = convolve(img,kernel)
		# Write to file
		img_2.data = conv
		p = filename.split('_')
		if p[1] == 'test':
			outPutFile = C_out+str(filename).strip('test_NH2.fits')+'_conv_1.fits'
		elif p[1] == '2':
			outPutFile = C_out+str(filename).strip('_2_test_NH2.fits')+'_conv_2.fits'
		elif p[1] == '3':
			outPutFile = C_out+str(filename).strip('_3_test_NH2.fits')+'_conv_3.fits'
		img_2.writeto(outPutFile)
			continue
	else:
		continue

















