########################################################################
#				Calculating 
#				A + (B-C) + (D-E)
########################################################################
#				Jul 4, 2018
########################################################################
import numpy as np
from astropy.io import fits
import math
import os
import montage_wrapper as montage
###
from astropy.utils.data import get_pkg_data_filename
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy.convolution import Gaussian1DKernel
######################################
#	Paths
######################################

# Input files
Adir = '/home/emma/gradu/data/All_data/n_tau_all/500_160_n_tau/reproject_to_250/'
BC_in = '/home/emma/gradu/data/All_data/n_tau_all/350_160_n_tau/B_minus_C/reproject_to_250/'
DE_in = '/home/emma/gradu/data/All_data/n_tau_all/250_160_n_tau/D_minus_E/'

# Output
outDir = '/home/emma/gradu/data/All_data/combinations/PACS_SPIRE/'

# Remove bad pixels
zerosOut = '/home/emma/gradu/data/All_data/combinations/PACS_SPIRE/removedBadPixels/'


################################################
#	DO NOT CONVOLVE A OR (B-C)!!!
################################################
def sigma(img,band):
	if band == 500:
		sigmaU = 35.4		# A FWHM
	elif band == 350:
		sigmaU = 24.2		# BC FWHM
	sigma250 = 17.9		# Result FWHM
	
	minus = sigmaU**2 - sigma250**2
	sq = math.sqrt(minus)
	
	CDELT1_deg = img.header['CDELT1']
	CDELT1 = abs(CDELT1_deg*3600.0)
	print 'C= ',CDELT1
	new_FWHM = sq / CDELT1
	sigma = new_FWHM/(math.sqrt(CDELT1*8*math.log(2)))
	return sigma


def convolve(filename):
	print 'NO DO NOT'
#	current = get_pkg_data_filename(filename)
#	img_2 = fits.open(current)[0]
#	img = img_2.data
#	sKernel2 = sigma(img_2)
#	# Create kernel and convolve
#	kernel = Gaussian2DKernel(sKernel2)
#	conv = convolve(img,kernel)
#	# Write to file
#	img_2.data = conv
#	return img_2


# A: G202.16+02.64_conv_rep_1.fits
# BC: G202.16+02.64_conv_rep_BC_1.fits
# DE: G202.16+02.64_B_C_1.fits

######################################
#	A + (B - C) + (D - E)
######################################
for filename in os.listdir(Adir):

	if filename.endswith('.fits'):
	
		f = fits.open(Adir+filename)[0]
	
		# Naming files
		p = filename.split('_')
		aFile = Adir+filename
		bcFile = BC_in+p[0]+'_rep_BC_'+p[2]
		deFile = DE_in+p[0]+'_B_C_'+p[2]
		outPut = outDir+p[0]+'_AE_'+p[2]
	
		# Open data
		a = fits.open(aFile)[0]
		bc = fits.open(bcFile)[0]
		de = fits.open(deFile)[0]
		
		# Adding
		plus = a.data + bc.data + de.data
		
		# Saving to output
		f.data = plus
		f.writeto(outPut,clobber=True)
	
	
	
	
		continue
		
# Old output name = G202.16+02.64_ABC_1.fits




for filename in os.listdir(outDir):
	if filename.endswith('.fits'):
		f = fits.open(outDir+filename)[0]
		fData = f.data
		f.data[np.isnan(f.data)] = 0
		f.data[np.isinf(f.data)] = 0
		f.data[f.data<0] = 0
		f.data[f.data>10**24] = 0
		outFile = zerosOut+p[0]+'_AE_zeros_'+p[2]
		f.writeto(outFile,clobber=True)



#S[S < 0] = 0




















