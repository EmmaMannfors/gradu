#################################################################################
#		Read data from extractclumps results
#################################################################################
#				Jul 31, 2018
#################################################################################


import os
from astropy.io import fits
import numpy as np



#############################
# paths
#############################
dir_in = '/home/emma/gradu/data/All_data/SCUBA_850/clumps/extractResults/'
data_in = '/home/emma/gradu/data/All_data/SCUBA_850/clumps/fits_clumps/'
data_out = '/home/emma/gradu/data/All_data/SCUBA_850/clumps/extractResults/results/'

#############################
# Init
#############################


#############################
# Functions
#############################

for filename in os.listdir(dir_in):
	#if filename.endswith('.FIT'):
	if filename.endswith('ttt'):
	#if filename.endswith('G001_36+20_96_850_cat.FIT'):
		field = filename.strip('_cat.FIT')
		#print field
		outName = field+'_clump_info.txt'
		fitsFile = field+'_clump.fits'

		#print data_in+fits
		###############
		# Find nClumps
		###############
		ft = fits.open(data_in+fitsFile)[0]
		nClumps = int(np.max(ft.data[np.isfinite(ft.data)]))
		#print nClumps
		###############
		# Get data
		###############
		f = fits.open(dir_in+filename)[1]


		###############
		# Write to file
		###############

		text = open(data_out+outName,'w')
		text.write(str(f.data.dtype))
		text.write('\n')
		for i in range(nClumps):
			text.write(str(f.data[i]))
			text.write('\n')
		






text = open(data_out+'header_meanings.txt','w')
filename = 'G001_36+20_96_850_cat.FIT'
f = fits.open(dir_in+filename)[1]









