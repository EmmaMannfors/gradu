######################################################################
# To draw clump contours next to SCUBA-2 images
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

#testDir = '/home/emma/gradu/data/All_data/SCUBA_850/clumps/testClumps/'


for filename in os.listdir(scuba_in):
	if filename.endswith('.fits'):

		#if filename != 'G001_36+20_96_850.fits':
		#	continue

		fig = plt.figure()
		f = aplpy.FITSFigure(scuba_in+filename,figure=fig,subplot=[0.1,0.1,0.35,0.35])
		f.show_colorscale(cmap='inferno',vmin=0)

		clump = clump_in+filename.strip('.fits')+'_clump.fits'

		c = aplpy.FITSFigure(clump,figure=fig,subplot=[0.55,0.1,0.35,0.35])
		c.show_colorscale(cmap='inferno',vmin=0)
		c.tick_labels.hide_y()
		c.axis_labels.hide_y()
		
		c.tick_labels.set_xformat('dd.d')
		f.tick_labels.set_xformat('dd.d')
		f.tick_labels.set_yformat('dd.dd')
		
		outName = filename.strip('.fits')+'_clumps.png'
		plt.savefig(outDir+outName,bbox_inches='tight')
		plt.close()


#fig.tick_labels.set_xformat('hh:mm:ss.ss')

#Set the format for the y-axis labels (e.g dd:mm, dd:mm:ss.s, etc.):

#fig.tick_labels.set_yformat('dd:mm:ss.s')


