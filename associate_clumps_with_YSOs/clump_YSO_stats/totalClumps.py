###############################################################################
#		 Get max number of clumps
###############################################################################
#				jul 30, 2018
###############################################################################

import os
from astropy.io import fits
import numpy as np

######################################
# Paths
######################################
fits_in = '/home/emma/gradu/data/All_data/SCUBA_850/clumps/fits_clumps/'
out = '/home/emma/gradu/data/All_data/YSO/clump_and_YSO_list/clumps_vs_all/'

outFile = open(out+'total_clumps.txt','w')
outFile.write('# field			clumps\n')
######################################
# Code
######################################

for filename in os.listdir(fits_in):
	if filename.endswith('.fits'):
		field = filename.strip('_clump.fits')
		f = fits.open(fits_in+filename)[0]
		fData = f.data
		#m = nonzero(np.isfinite[fData])
		m = np.max(fData[np.isfinite(fData)])
		ma = int(m)
		print field
		print ma
		outFile.write(str(field)+'	'+str(ma)+'\n')


#m  = nonzero((F1[0].data!=0.0)&(isfinite(F1[0].data))&(F2[0].data!=0.0)&(isfinite(F2[0].data))&(F3[0].data!=0.0)&(isfinite(F3[0].data)))
