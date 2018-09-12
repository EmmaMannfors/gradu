################################################################################
#	Takes the average value of the pixels in an empty(-ish) space 
#	and subtracts this value from the entire image
################################################################################
#		Uses the images in units of MJy/sr
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
coordinateFile = 'testCoordinates.txt'

testdir = '/home/emma/gradu/codes/subtract_average/test/'

matchPath = '/home/emma/gradu/codes/matches/'

means = open('testMeans.txt','w')
means.write('SCUBA-2		mean\n')
########################################
#	Functions
########################################
def ra_to_deg(ra):
	g = ra.split(':')
	h = float(g[0])
	m = float(g[1])
	s = float(g[2])
	h += (s/3600.0)
	h += (m/60.0)
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
	

def change(filename):
	if filename == 'NA':
		return 'NA'
	if filename.endswith('_crop.fits'):
		n = filename.strip('_crop.fits')+'_units_1.fits'
		return n
	elif filename.endswith('_crop_2.fits'):
		n = filename.strip('_crop_2.fits')+'_units_2.fits'
		return n
	elif filename.endswith('_crop_3.fits'):
		n = filename.strip('_crop_3.fits')+'_units_3.fits'
		return n
	elif filename.endswith('_crop_4.fits'):
		n = filename.strip('_crop_4.fits')+'_units_4.fits'
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

def pixel_to_wcs(filename):
	if filename != 'NA':
		f = aplpy.FITSFigure(filename)
		fData = fits.open(filename)[0].data
	#	print fData
	#	print fData.shape
	#	f.pixel2world(fData)
	#	print filename
#table of all pixel coordinates, out comes a table of all pixel coordinates in wcs
		J,I = np.indices(fData.shape)
		ra,dec = f.pixel2world(J,I)
		f.close()
		return ra,dec

def outName_test(filename):
	if filename.endswith('crop.fits'):
		outName = filename.strip('_crop.fits')+('_norm_1.fits')
	elif filename.endswith('units_1.fits'):
		outName = filename.strip('_units_1.fits')+('_norm_1.fits')
	elif filename.endswith('crop_2.fits'): 
		outName = filename.strip('_crop_2.fits')+('_norm_2.fits')
	elif filename.endswith('units_2.fits'):
		outName = filename.strip('_units_2.fits')+('_norm_2.fits')
	elif filename.endswith('crop_3.fits'): 
		outName = filename.strip('_crop_3.fits')+('_norm_3.fits')
	elif filename.endswith('units_3.fits'):
		outName = filename.strip('_units_3.fits')+('_norm_3.fits')
	else:
		outName = filename.strip('_units_4.fits')+('_norm_4.fits')
	return outName

########################################
#	Add matches to dictionary
########################################
d = dict()

with open(matchPath+'cropMatches.txt') as matches:
	next(matches)
	for line in matches:
		t = line.split()
	#	print t[0]
		d[str(t[0])] = [t[1],t[2],t[3],t[4],t[5],t[6]]




########################################
#	Run code
########################################
locList = []

i = 0
with open(coordinateFile) as co:
	next(co)
	for line in co:
		# Getting all the matches
		p = line.split()
		scuba = p[0]+'_850.fits'
		matchList = d[scuba]
		pacs70 = testdir+change(matchList[0])
		pacs100 = testdir+change(matchList[1])
		pacs160 = testdir+change(matchList[2])
		psw = testdir+matchList[3]
		pmw = testdir+matchList[4]
		plw = testdir+matchList[5]
		SCUBA = testdir+scuba
		locList = [pacs70,pacs100,pacs160,psw,pmw,plw]
	#	for item in locList:
	#		print item,'\n'
				# Getting coordinates
		ra, dec, r_deg = p[1], p[2], p[3]		# radius in degrees
		if ra == 'NA':
			continue
		center = (ra,dec)
		# get values for every pixel inside circle, and then take average (np.avg? np.mean?)
		# Subtract this value from every pixel inside the image
		# Write normalized file to output
		
		center_deg = (ra_to_deg(ra),dec_to_deg(dec))
		
		ra0,dec0 = center_deg
		
		for filename in locList:
			if filename.endswith('NA'):
				continue
			f = fits.open(filename)[0]
			fData = f.data
			RA,D = pixel_to_wcs(filename)
			# Calculate mean
			dist = distance_on_sphere_deg(ra0,dec0,RA,D)
			m = np.nonzero(dist<0.02)
			
			
			e = np.mean(f.data[m])
			means.write(p[0]+'	'+str(e)+'\n')
			print e
			# 
			
			print f.data.shape
		#	print fData
			for val in np.nditer(fData, op_flags=['readwrite']):
				val -= e
				
			f.data = fData
			print filename
			outName = outName_test(filename)
			print outName
			
		#	print f.data
			f.writeto(outName)
			#f.close()
			print 5*'*****'
			
			
	

means.close()

# ra0, de0 = center
# dist = distance_on_sphere_deg(ra0,dec0,ra,dec)
#m= nonzero(dist<0.02)
#mean(f[0].data[m])


#"


