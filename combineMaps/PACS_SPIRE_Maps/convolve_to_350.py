########################################################################
#			Convolving the N maps from the 
#			calculate_N code for 500, 350, 
#			250  and 160 um and convolving them to
#			350 um resolution
########################################################################
#				Jun 15, 2018
########################################################################
#		Needs to be run in same directory as input files
#		REAL CODE IN DIR_IN
######################################################################## 


NOT BEEN EDITED YET!!!!

import numpy as np
from astropy.io import fits
import montage_wrapper as montage
import math
import os

from astropy.utils.data import get_pkg_data_filename
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy.convolution import Gaussian1DKernel


###################################
#	Calculate sigma
#	Same for all images
###################################

# Should this be the other way around? 
# If so, what to do with sqrt(-x)?
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
	sigma = new_FWHM/(math.sqrt(CDELT1*math.log(2)))
#	sigma = new_FWHM/(math.sqrt(CDELT1*8*math.log(2)))
#	print sigma
	return sigma


###################################
#    Paths to convolve folders
###################################
#Directory of files
dir_in = '/home/emma/gradu/data/All_data/n_tau_all/n_tau_maps'

#Location of initial files
path_in = '/home/emma/gradu/data/All_data/n_tau_all/n_tau_maps/'

#Location of output files
conv_out = '/home/emma/gradu/data/All_data/n_tau_all/n_tau_maps/convolved_to_350/'

#Header files
hdrFiles = '/home/emma/gradu/data/All_data/SPIRE_PMW_350_Herschel/headers_350/'


###################################
#	Function for convolution
###################################
def band_350():
	for filename in os.listdir(dir_in): 
		if filename.endswith('NH2.fits'):
			# Get data
			current = get_pkg_data_filename(path_in+str(filename))
			img_2 = fits.open(current)[0]
			img = img_2.data
			sKernel2 = sigma(img_2)
		#	sKernel350 = 0.981339795226
			# Create kernel and convolve
			kernel = Gaussian2DKernel(sKernel2)
			conv = convolve(img,kernel)
			# Write to file
			img_2.data = conv
			p = filename.split('_')
			if p[1] == 'test':
				outPutFile = conv_out+str(filename).strip('_test_NH2.fits')+'_conv.fits'
			elif p[1] == '2':
				outPutFile = conv_out+str(filename).strip('_2__test_NH2.fits')+'_conv_2.fits'
			elif p[1] == '3':
				outPutFile = conv_out+str(filename).strip('_3__test_NH2.fits')+'_conv_3.fits'
			else: 
				print 'ERROR ABORT BEEP BEEP'
			img_2.writeto(outPutFile,clobber=True)
			continue
		else:
			continue

###################################
#	Reprojection
###################################
L25list = ['G017.69-00.15','G070.07-01.60','G094.11-04.86','G158.55-20.92','G178.31+00.28']
def rep():
	for filename in os.listdir(conv_out):
#	for filename in os.listdir(hdrFiles):
		if filename.endswith('.fits'):
#		if filename.endswith('.txt'):


			inFile = conv_out+filename
			p = filename.split('_')
			if p[1] == 'rep':
				continue
			l = 'L2.0'
			for i in range(0,len(L25list)):
				now = L25list[i]
				if p[0] == now:
					l = 'L2.5'
					break
		#	print p[0],'	',l
			if filename.endswith('conv.fits'):
				hdr = hdrFiles+p[0]+'_PMW_'+l+'.txt'
				outPut = conv_out+filename.strip('_conv.fits')+'_rep_1.fits'
			elif filename.endswith('_conv_2.fits'):
				hdr = hdrFiles+p[0]+'_PMW_'+l+'_2.txt'
				outPut = conv_out+filename.strip('_conv_2.fits')+'_rep_2.fits'
			elif filename.endswith('_conv_3.fits'):
				hdr = hdrFiles+p[0]+'_PMW_'+l+'_3.txt'
				outPut = conv_out+filename.strip('_conv_3.fits')+'_rep_3.fits'
			else: 
				print filename
				print 'AAAA EVERYBODY PANIC THERE IS AN ERROR IN THE CODE'
			montage.reproject(inFile,outPut,header=hdr,exact_size=True)

###################################
#	Run functions
###################################
do = raw_input('What to do? (conv/rep/both)\n')
if do == 'both':
	band_350()
	rep()
elif do == 'conv':
	band_350()
elif do == 'rep':
	rep()
else: 
	print 'cat says meow'




for filename in os.listdir(conv_out):
	q = filename.split('_')
	if q[1] == 'rep':
		f = fits.open(conv_out+filename)[0]
		if f.data.shape != (74,73):
			print filename,'	',f.data.shape
		else: 
			print 'meow'













