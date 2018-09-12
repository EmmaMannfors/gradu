########################################################################
#			Same as convolve_to_350.py, but 
#			with T and I250 and tau maps
########################################################################
#				Aug 7, 2018
########################################################################
#	DO NOT CONVOLVE A MAPS!!!!!!!!!!! ONLY REPROJECT!!
######################################################################## 
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
#	sigma = new_FWHM/(math.sqrt(CDELT1*math.log(2)))
	sigma = new_FWHM/(math.sqrt(CDELT1*8*math.log(2)))
#	print sigma
	return sigma


###################################
#    Paths to convolve folders
###################################
#Directory of files
dir_in = '/home/emma/gradu/data/All_data/n_tau_all/n_tau_maps/'

#Location of initial files
path_in = '/home/emma/gradu/data/All_data/n_tau_all/n_tau_maps/'

#Location of output files
conv_out = '/home/emma/gradu/data/All_data/n_tau_all/n_tau_maps/reprojected_to_350/'

#Header files
hdrFiles = '/home/emma/gradu/data/All_data/SPIRE_PMW_350_Herschel/headers_350/'


###################################
#	Reprojection
###################################
L25list = ['G017.69-00.15','G070.07-01.60','G094.11-04.86','G158.55-20.92','G178.31+00.28']
def rep():
	for filename in os.listdir(dir_in):
#	for filename in os.listdir(hdrFiles):
		if filename.endswith('.fits'):
#		if filename.endswith('.txt'):

			if filename.endswith('NH2.fits'):
				continue		#You already did these
			if filename.endswith('BETA.fits'): #Not necessary
				continue	

			inFile = dir_in+filename
			p = filename.split('_')
			if p[1] == 'rep':
				continue
			l = 'L2.0'
			for i in range(0,len(L25list)):
				now = L25list[i]
				if p[0] == now:
					l = 'L2.5'
					break


			field = p[0]
			if p[1] == 'test':
				num = '1'
				typ=p[2]
				hdr = hdrFiles+field+'_PMW_'+l+'.txt'
				if typ == 'tau.fits':
					end='t'+typ.strip('.fits')+'_rep_'+num+'.fits'	
				else:
					end=typ.strip('.fits')+'_rep_'+num+'.fits'
			else: 
				num = p[1]
				hdr = hdrFiles+p[0]+'_PMW_'+l+'_'+num+'.txt'
				typ=p[3]
				if typ == 'tau.fits':
					end = 't'+typ.strip('.fits')+'_rep_'+num+'.fits'
				else: 
					end = typ.strip('.fits')+'_rep_'+num+'.fits'


			outPut=conv_out+field+'_'+end
			outNAME=field+'_'+end


			
			
			print filename, '\n', outNAME, '\n', hdr
		#	print p[0],'	',l

		#	if filename.endswith('_2_test_NH2.fits'):
		#		hdr = hdrFiles+p[0]+'_PMW_'+l+'_2.txt'
		#		outPut = conv_out+filename.strip('_2_test_NH2.fits')+'_rep_2.fits'
		#	elif filename.endswith('_3_test_NH2.fits'):
		#		hdr = hdrFiles+p[0]+'_PMW_'+l+'_3.txt'
		#		outPut = conv_out+filename.strip('_3_test_NH2.fits.fits')+'_rep_3.fits'
		#	elif filename.endswith('NH2.fits'):
		#		hdr = hdrFiles+p[0]+'_PMW_'+l+'.txt'
		#		outPut = conv_out+filename.strip('_test_NH2.fits')+'_rep_1.fits'
		#	else: 
		#		print filename
		#		print 'AAAA EVERYBODY PANIC THERE IS AN ERROR IN THE CODE'
		#	print filename, '	',hdr
			print 5*'*'
			montage.reproject(inFile,outPut,header=hdr,exact_size=True)
			print 10*'*'

###################################
#	Run functions
###################################
#do = raw_input('What to do? (conv/rep/both)\n')
#if do == 'both':
#	band_350()
#	rep()
#elif do == 'conv':
#	band_350()
#elif do == 'rep':
rep()
#else: 
#	print 'cat says meow'




for filename in os.listdir(conv_out):
	q = filename.split('_')
	#if q[1] == 'rep':
	#	f = fits.open(conv_out+filename)[0]
	#	if f.data.shape != (74,73):
	#		print filename,'	',f.data.shape
	#	else: 
	#		print 'meow'
	if filename.endswith('.fits'):
		f = fits.open(conv_out+filename)[0]
		if f.data.shape != (74,73):
			print filename,'	',f.data.shape
		else: 
			print 'meow'












