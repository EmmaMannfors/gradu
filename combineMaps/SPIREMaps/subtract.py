########################################################################
#				Calculating B - C
#			Where B = 250 + 350 um map
#		And C = B convolved to resolution of 500 um map
########################################################################
#			FIRST reproject A maps to 350 um
########################################################################
#			Then calculating A + (B - C)
#				Jul 5, 2018
########################################################################
import numpy as np
from astropy.io import fits
import math
import os
import montage_wrapper as montage

######################################
#	Paths
######################################


Adir = '/home/emma/gradu/data/All_data/n_tau_all/n_tau_maps/reprojected_to_350/'
Bdir = '/home/emma/gradu/data/All_data/n_tau_all/two_map_N_tau/'
Cdir = '/home/emma/gradu/data/All_data/n_tau_all/two_map_N_tau/convolved_to_500/'


Apath = Adir
Bpath = Bdir
Cpath = Cdir

A_B_C_out = '/home/emma/gradu/data/All_data/combinations/A_plus_B_min_C/'
B_C_out = '/home/emma/gradu/data/All_data/n_tau_all/two_map_N_tau/B_minus_C/'

######################################
#	B - C
######################################
missedList = []

for filename in os.listdir(Cdir):
##	if filename == 'G006.04+36.73_conv_1.fits':
	if filename.endswith('.fits'):
		c = fits.open(Cpath+str(filename))[0]
		cData = c.data
		p = filename.split('_')
		if filename.endswith('conv_1.fits'):
			bFile = p[0] + '_test_NH2.fits'
			outPutFile = B_C_out + p[0] + '_B_C_1.fits'
		elif filename.endswith('conv_2.fits'):
			bFile = p[0] + '_2_test_NH2.fits'
			outPutFile = B_C_out + p[0] + '_B_C_2.fits'
		elif filename.endswith('conv_3.fits'):
			bFile = p[0] + '_3_test_NH2.fits'
			outPutFile = B_C_out + p[0] + '_B_C_3.fits'
		elif filename.endswith('conv_4.fits'):
			bFile = p[0] + '_4_test_NH2.fits'
			outPutFile = B_C_out + p[0] + '_B_C_4.fits'
		elif filename.endswith('conv_5.fits'):
			bFile = p[0] + '_5_test_NH2.fits'
			outPutFile = B_C_out + p[0] + '_B_C_5.fits'
		else:
			print 'wut'
			missedList.append(filename)
		b = fits.open(Bpath+bFile)[0]
		bData = b.data
		minus = bData - cData
		c.data = minus
#		print minus
##		write = raw_input('write to file? (y/n)')
#		write = 'n'
#		if write == 'y':
		c.writeto(outPutFile,clobber=True)
		

print missedList
########################################################################
#			Calculating A + (B - C)
#		where A = the 250, 350, and 500 um map
#			And (B - C) is calculated above
########################################################################
resList = []
for filename in os.listdir(B_C_out):
	if filename.endswith('res_1.fits'):
		resList.append(filename)
	elif filename.endswith('res_2.fits'):
		resList.append(filename)
	elif filename.endswith('res_3.fits'):
		resList.append(filename)
	elif filename.endswith('res_4.fits'):
		resList.append(filename)
	elif filename.endswith('res_5.fits'):
		resList.append(filename)
	else: 
		print 'uhh'

i = 0

for filename in os.listdir(B_C_out):
	if filename.endswith('.fits'):
		bc = fits.open(B_C_out + str(filename))[0]
		bcData = bc.data
		p = filename.split('_')
		if filename.endswith('_1.fits'):
			aFile = p[0] + '_rep_1.fits'
			outPut = p[0] + '_ABC_1.fits'
		elif filename.endswith('_2.fits'):
			aFile = p[0] + '_rep_2.fits'
			outPut = p[0] +	'_ABC_2.fits'		
		elif filename.endswith('_3.fits'):
			aFile = p[0] + '_rep_3.fits'
			outPut = p[0] +	'_ABC_3.fits'		
		elif filename.endswith('_4.fits'):
			aFile = p[0] + '_rep_4.fits'
			outPut = p[0] +	'_ABC_4.fits'		
		elif filename.endswith('_5.fits'):
			aFile = p[0] + '_rep_5.fits'
			outPut = p[0] +	'_ABC_5.fits'		
		else:
			print 'meow',filename
			i += 1
			continue

		a = fits.open(Apath+aFile)[0]
		aData = a.data
		plus = aData + bcData
		bc.data = plus
		bc.writeto(A_B_C_out + outPut,clobber=True)
		continue



#print i # i counts number of missed images
if i == 0:
	print 'jee'












