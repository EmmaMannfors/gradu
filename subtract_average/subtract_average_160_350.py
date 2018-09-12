################################################################################
#	Takes the average value of the pixels in an empty(-ish) space 
#	and subtracts this value from the entire image
################################################################################
#		Uses the images in units of MJy/sr
################################################################################
# 		FOR 350-CONVOLVED IMAGES
################################################################################
#			Jun 26, 2018
################################################################################
#execfile('Aux.py')

import aplpy
from astropy.io import fits
import numpy as np

from functions import distance_on_sphere_deg

########################################
#	Paths & files
########################################
#Assuming files are in this folder
coordinateFile = 'emptySpaces.txt'
PACS_coordinateFile = 'emptySpaces_only_PACS.txt'


in_160 = '/home/emma/gradu/data/All_data/PACS_160_Herschel/reproject_160_to_350/'
in_250 = '/home/emma/gradu/data/All_data/SPIRE_PSW_250_Herschel/TwoMapReproject_250/'
in_350 = '/home/emma/gradu/data/All_data/SPIRE_PMW_350_Herschel/cropped_350/'
in_850 = '/home/emma/gradu/data/All_data/SCUBA_850/'

matchPath = '/home/emma/gradu/codes/matches/'

means = open('Means.txt','w')
means.write('SCUBA		mean\n')
########################################
#	Functions
########################################
def ra_to_deg(ra):
	g = ra.split(':')
	h = float(g[0])
	m = float(g[1])
	s = float(g[2])
	h += (m/60.0)
	h += (s/3600.0)
	deg = h*(360.0/24.0)
	return round(deg,6)

def dec_to_deg(dec):
	y = dec.split(':')
	deg,amin,asec = float(y[0]),float(y[1]),float(y[2])
	if deg < 0:
		deg -= (amin/60.0)
		deg -= (asec/3600.0)
		return round(deg,6)
	deg += (amin/60.0)
	deg += (asec/3600.0)
	return round(deg,6)
	

def change_160(filename):
	if filename == 'NA':
		return 'NA'
	if filename.endswith('_crop.fits'):
		n = filename.strip('_crop.fits')+'_res_1.fits'
		return n
	elif filename.endswith('_crop_2.fits'):
		n = filename.strip('_crop_2.fits')+'_res_2.fits'
		return n
	elif filename.endswith('_crop_3.fits'):
		n = filename.strip('_crop_3.fits')+'_res_3.fits'
		return n
	elif filename.endswith('_crop_4.fits'):
		n = filename.strip('_crop_4.fits')+'_res_4.fits'
		return n

def change(filename):
	if filename == 'NA':
		return 'NA'
	if filename.endswith('_crop.fits'):
		n = filename.strip('_crop.fits')+'_rep_1.fits'
		return n
	elif filename.endswith('_crop_2.fits'):
		n = filename.strip('_crop_2.fits')+'_rep_2.fits'
		return n
	elif filename.endswith('_crop_3.fits'):
		n = filename.strip('_crop_3.fits')+'_rep_3.fits'
		return n
	elif filename.endswith('_crop_4.fits'):
		n = filename.strip('_crop_4.fits')+'_rep_4.fits'
		return n


def path(band):
	if band == 70:
		path = in_70
	elif band == 100:
		path = in_100
	elif band == 160:
		path = in_160
	elif band == 250:
		path = in_250
	elif band == 350:
		path = in_350
	elif band == 500:
		path = in_500
	elif band == 850:
		path = in_850
	return path

#def pixel_to_wcs(filename):
def pixel_to_wcs(filename):
	if filename != 'NA':
		fig = aplpy.FITSFigure(filename)
		fData = fits.open(filename)[0].data
	#	print fData
	#	print fData.shape
	#	f.pixel2world(fData)
	#	print filename
#table of all pixel coordinates, out comes a table of all pixel coordinates in wcs
	#	J,I = np.indices(fData.shape)
		I,J = np.indices(fData.shape)
		ra,dec = fig.pixel2world(J,I)
		fig.close()
		return ra,dec

def outName(filename):
	if filename.endswith('res_1.fits'):
		outName = filename.strip('_res_1.fits')+('_norm_1.fits')
	elif filename.endswith('rep_1.fits'):
		outName = filename.strip('_rep_1.fits')+('_norm_1.fits')
	elif filename.endswith('res_2.fits'): 
		outName = filename.strip('_res_2.fits')+('_norm_2.fits')
	elif filename.endswith('rep_2.fits'):
		outName = filename.strip('_rep_2.fits')+('_norm_2.fits')
	elif filename.endswith('res_3.fits'): 
		outName = filename.strip('_res_3.fits')+('_norm_3.fits')
	elif filename.endswith('rep_3.fits'):
		outName = filename.strip('_rep_3.fits')+('_norm_3.fits')
	elif filename.endswith('res_4.fits'): 
		outName = filename.strip('_units_4.fits')+('_norm_4.fits')
	elif filename.endswith('crop.fits'):
		outName = filename.strip('_crop.fits')+('_norm_1.fits')
	elif filename.endswith('crop_2.fits'):
		outName = filename.strip('_crop_2.fits')+('_norm_2.fits')
	elif filename.endswith('crop_3.fits'):
		outName = filename.strip('_crop_3.fits')+('_norm_3.fits')
	else:
		outName = 'AAA.fits'
	#	outName = filename.strip('_units_4.fits')+('_norm_4.fits')
	return outName


#Dropbox: Readme
normMatches = open('normalized_matches_350_160.txt','w')
normMatches.write('PACS 160		SPIRE 250		SPIRE 350\n')
########################################
#	Add matches to dictionary
########################################
d = dict()

with open(matchPath+'cropMatches.txt') as matches:
	next(matches)
	for line in matches:
		t = line.split()
	#	print t[0]
		if t[4] == 'G013.90-00.51_PSW_L2.0_crop.fits':
			continue
		elif t[4] == 'G035.19-01.75_PSW_L2.0_crop.fits':
			continue
		else:
			d[str(t[0])] = [t[1],t[2],t[3],t[4],t[5],t[6]]
			continue




########################################
#	Run code
########################################


with open(coordinateFile) as co:
	next(co)
	for line in co:
		# Getting all the matches
		locList = []
		p = line.split()
		scuba = p[0]+'_850.fits'
		if scuba == 'G035.38-01.77_850.fits':
			continue
		if scuba == 'G014.15-00.55_850.fits':
			continue


	####
	# TEST
#
#		if scuba != 'G039.74+01.98_850.fits':
#			continue
	####


		matchList = d[scuba]

		pacs160 = in_160+change_160(matchList[2])
		psw = in_250+change(matchList[3])
		pmw = in_350+matchList[4]
		SCUBA = in_850+scuba
		locList = [pacs160,psw,pmw]

		# Getting coordinates
		ra, dec, r_deg = p[1], p[2], p[3]		# radius in degrees
		if ra == 'NA':
			continue

		
		center = (ra,dec)

		
		center_deg = (ra_to_deg(ra),dec_to_deg(dec))
		
		ra0,dec0 = center_deg


		print center
		print center_deg

		#ra0 = ra_to_deg(ra)
		#dec0 = dec_to_deg(dec)

		for filename in locList:
			if filename.endswith('NA'):
				continue


		####
		# TEST
		####
#			if not filename.endswith('_160_L2.5_J_res_1.fits'):
#				continue
		####


			f = fits.open(filename)[0]
			fData = f.data
	
######################
	#		fig = aplpy.FITSFigure(filename)
	#		J,I = np.indices(f.data.shape)
	#		RA,D = fig.pixel2world(J,I)
	#		fig.close()





			RA,D = pixel_to_wcs(filename)
######################			
			

			

			# Calculate mean
			dist = distance_on_sphere_deg(ra0,dec0,RA,D)


			m = np.nonzero((dist<0.02)&(np.isfinite(f.data)))

		#	m = np.nonzero((dist<0.02)&(np.isfinite(fData)))

		#	m = np.nonzero(np.isfinite(f.data))
		#	m = np.nonzero(dist<0.02)
		#	if scuba == 'G039.74+01.98_850.fits' and filename.endswith('_160_L2.5_J_res_1.fits'):
		#		print f.data[m]
				
	#		print f.data.shape
	#		print f.data[m].shape


			e = np.mean(f.data[m])
		#	print e
			f.data -= e

		#	fData -= np.mean(fData[m])

		#	f.data = fData

			out = outName(filename)

			f.writeto(out,clobber=True)
			normMatches.write(out+'	')
			print 5*'*'
			
		normMatches.write('\n')
		print 5*'*****'
		continue
		
		




# ra0, de0 = center
# dist = distance_on_sphere_deg(ra0,dec0,ra,dec)
#m= nonzero(dist<0.02)
#mean(f[0].data[m])




