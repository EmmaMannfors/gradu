######################################################################
# To draw clump contours on SCUBA-2 images
######################################################################

import numpy as np
import aplpy
import os
import cmd
import matplotlib.pyplot as plt
from astropy.io import fits
##########################
# Paths
##########################
clump_in = '/home/emma/gradu/data/All_data/SCUBA_850/clumps/fits_clumps/'
scuba_in = '/home/emma/gradu/data/All_data/SCUBA_850/changedNames/'

outDir = '/home/emma/gradu/data/All_data/SCUBA_850/clumpContours/'

testDir = '/home/emma/gradu/data/All_data/SCUBA_850/clumps/testClumps/'
##########################
# 
##########################

for filename in os.listdir(scuba_in):
	if filename.endswith('.fits'):

		#if filename != 'G001_36+20_96_850.fits':
		#	continue
		
		fig = plt.figure()
		f = aplpy.FITSFigure(scuba_in+filename,figure=fig)
		f.show_colorscale(cmap='inferno',vmin=0)
		clump = filename.strip('.fits')+'_clump.fits'

		#d = fits.open(clump_in+clump)[0]
		#dData = d.data
		#dData[np.isinf(dData)] = 0
		#d.data = dData
		#d.writeto(testDir+clump,clobber=True)
		
		#level = 5
		f.show_contour(clump_in+clump,colors='white')#,levels=level,colors='white',smooth=1)
		
		outName = filename.strip('.fits')+'_clumps.png'
		plt.savefig(outDir+outName,bbox_inches='tight')
		plt.close()
		

