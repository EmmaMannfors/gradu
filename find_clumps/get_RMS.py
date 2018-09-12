#################################################################################
#		Perform FellWalker on the 850 um files
#			Convert the results into .fits files
#################################################################################
#				Jul 13, 2018
#################################################################################
# 	Each instance of os.system() makes a new shell, and so Starlink-
#			commands must be run in each excecution
#################################################################################

import os
from astropy.io import fits
import numpy as np
import aplpy


# Getting the path thingy to work
import sys
sys.path.append('/home/emma/star-2017A/lib/')
os.environ['PATH'] += os.path.pathsep + '/home/emma/star-2017A/lib/'
os.environ['PATH'] += os.pathsep + '/home/emma/star-2017A/lib/'


from functions import distance_on_sphere_deg,ra_to_deg,dec_to_deg
from FixFits import FixFits
#############################
# paths
#############################
scubaFits_in = '/home/emma/gradu/data/All_data/SCUBA_850/'
in_850 = '/home/emma/gradu/data/All_data/SCUBA_850/sdfFIles/'
out_850 = '/home/emma/gradu/data/All_data/SCUBA_850/clumps/'

#############################
# Write RMS values to file
#############################
rmsVal = open('rms.txt','w')
rmsVal.write('SCUBA		rms\n')


#############################
# Functions
#############################
def outName(name):
	out = name.strip('.sdf')+'_clumps.sdf'
	return out

def outFits(name):
	out = name.strip('.sdf')+'_clumps.fits'
	return out


#def pixel_to_wcs(filename):
#	if filename != 'NA':
#		fig = aplpy.FITSFigure(filename)
#		fData = fits.open(filename)[0].data
#		S = fData.reshape(fData.shape[1],fData.shape[2])
#
#		I,J = np.indices(S.shape)
#		print I.shape # = (243, 241)
#		print J.shape # = (243, 241)
#
#		# THIS IS WHERE IT FUCKS UP! 
#		ra,dec = fig.pixel2world(J,I)
#
#	#	fig.close()
#		#return ra,dec
#		return 1,2

		
i = 0
#############################
# Get FellWalker arguments
#############################
# Get RMS
with open('centers_scuba.txt') as c: 
	next(c)
	for line in c:
		p = line.split()
		e = p[0].split('.')
		field = e[0]+'_'+e[1]+'_'+e[2]
		name = field+'_850.sdf'
		scuba = p[0]+'_850.fits'
	#	print name,'\n',scuba

		ra = p[1]
		dec = p[2]

		rms = p[3]

	#	if field != 'G001_36+20_96': #and field != 'G005_93-01_00':	
	#		continue

	#	if  field != 'G005_93-01_00':	
	#		continue

	
		center_deg = (ra_to_deg(ra),dec_to_deg(dec))
	
		ra0,dec0 = center_deg

		print ra,' ',dec
		print center_deg



		scubaFile_3D = scubaFits_in+scuba
		f = fits.open(scubaFile_3D)#[0]
		sc = FixFits(f)   # sc is same type as f

		

		fData = sc[0].data

		
		fig = aplpy.FITSFigure(sc)
		I,J = np.indices(fData.shape)
		RA,D = fig.pixel2world(J,I)



		dist = distance_on_sphere_deg(ra0,dec0,RA,D)
		m = np.nonzero((dist<0.02)&(np.isfinite(fData)))

		e = np.std(fData[m])

		rmsVal.write(p[0]+'	'+str(e)+'\n')
		print e
		print 5*'*'
		fig.close()



def getRMS(p,ra,dec):
	e = p.split('.')
	field = e[0]+'_'+e[1]+'_'+e[2]
	name = field+'_850.sdf'
	scuba = p+'_850.fits'
	#	print name,'\n',scuba

#	ra = p[1]
#	dec = p[2]



	#	if field != 'G001_36+20_96': #and field != 'G005_93-01_00':	
	#		continue

	#	if  field != 'G005_93-01_00':	
	#		continue

	
	center_deg = (ra_to_deg(ra),dec_to_deg(dec))
	
	ra0,dec0 = center_deg

	print ra,' ',dec
	print center_deg


	scubaFile_3D = scubaFits_in+scuba
	f = fits.open(scubaFile_3D)#[0]
	sc = FixFits(f)   # sc is same type as f

		

	fData = sc[0].data

	
	fig = aplpy.FITSFigure(sc)
	I,J = np.indices(fData.shape)
	RA,D = fig.pixel2world(J,I)



	dist = distance_on_sphere_deg(ra0,dec0,RA,D)
	m = np.nonzero((dist<0.02)&(np.isfinite(fData)))

	rms = np.std(fData[m])	
	rmsVal.write(p[0]+' '+str(rms)+'\n')
	fig.close()
	return rms
	#print 5*'*'






	#	f = fits.open(scubaFile)[0]
	#	print fData.shape
	#	fD = fData.reshape(fData.shape[1],fData.shape[2])  # THE APLPY DATA NEEDS TO BE RESHAPED TOO 
	


#####################################
		# just make a temp.fits file! 
	#	f.data = fD
	#	f.writeto('temp.fits',clobber=True)
	#	f = fits.open('temp.fits')[0]
	#	fData = f.data
		


#####################################
		# THIS IS WHERE IT FUCKS UP
	#	RA,D = pixel_to_wcs(scubaFits_in+scuba)


#####################################
		
		# I RESHAPE F.DATA, NOT FIG.DATA!!

	#	fig = aplpy.FITSFigure(scubaFits_in+scuba,dimensions=[0,1],slices=[1])

#The cube dimensions
 #   are:

  #         0 RA---TAN 241
   #        1 DEC--TAN 243
    #       2 FREQ-W2F 1
# Slices should be array of shape n-1 = 3-1 = 1

	#	I,J = np.indices(fD.shape)#,dtype='float')

	#	A,I,J = np.indices(fData.shape)

	#	print J
	#	print A.shape
	#	print I.shape # = (243, 241)
	#	print J.shape # = (243, 241)

	#	i = I.reshape(I.shape[1],I.shape[2])
	#	j = J.reshape(J.shape[1],J.shape[2])
		
	#	print i.shape # = (243, 241)
	#	print j.shape # = (243, 241)	

		# THIS IS WHERE IT FUCKS UP! 
	#	RA,D = fig.pixel2world(J,I)

	#	RA,D = fig.pixel2world(I,A)
		
#####################################



