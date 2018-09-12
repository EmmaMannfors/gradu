##!/usr/bin/sh
#!/usr/bin/python

#################################################################################
#		Perform FellWalker on the 850 um files
#			Convert the results into .fits files
#################################################################################
#				Jul 13, 2018
#################################################################################
# 	Each instance of os.system() makes a new shell, and so Starlink-
#			commands must be run in each excecution
#################################################################################

#!/bin/dash
#!/bin/bash
#!/bin/sh
#!/usr/bin/python
#!/usr/bin/env python
#!/bin/tcsh

import os
from astropy.io import fits
import numpy as np
import aplpy


# Getting the path thingy to work
import sys
sys.path.append('/home/emma/star-2017A/lib/')#libkpg_adam.so.0.0.0')
#sys.path.append('/home/emma/star-2017A/lib/libkpg_adam.so.0')
#sys.path.append('/home/emma/star-2017A/lib/libkpg_adam.so')
os.environ['PATH'] += os.path.pathsep + '/home/emma/star-2017A/lib/'
os.environ['PATH'] += os.pathsep + '/home/emma/star-2017A/lib/'
#>>> sys.path.append("/opt/local/bin")
#>>> os.system("wget")
#sh: wget: command not found
#32512
#>>> os.environ['PATH'] += os.pathsep + '/opt/local/bin'
#>>> os.system("wget")
#wget: missing URL
#locate libkpg_adam.so.0
#/home/emma/star-2017A/lib/libkpg_adam.so.0
#/home/emma/star-2017A/lib/libkpg_adam.so.0.0.0


#from FixFits import FixFits
#from functions import distance_on_sphere_deg,ra_to_deg,dec_to_deg
#from get_RMS import getRMS
#############################
# paths
#############################
scubaFits_in = '/home/emma/gradu/data/All_data/SCUBA_850/'
in_850 = '/home/emma/gradu/data/All_data/SCUBA_850/sdfFIles/'
out_850 = '/home/emma/gradu/data/All_data/SCUBA_850/clumps/'


#############################
# Functions
#############################
def outName(name):
	out = name.strip('.sdf')+'_clumps.sdf'
	return out

def outFits(name):
	out = name.strip('.sdf')+'_clumps.fits'
	return out


def pixel_to_wcs(filename):
	if filename != 'NA':
		fig = aplpy.FITSFigure(filename)
		fData = fits.open(filename)[0].data
		# fData is (1,228,229) array: reshape into (228,229) array
		S = fData.reshape(fData.shape[1],fData.shape[2])

		# To go back, if needed
		# S = S.reshape((1,S.shape[0],S.shape[1]))
		I,J = np.indices(S.shape)

		ra,dec = fig.pixel2world(J,I)
		fig.close()
		return ra,dec

def getRMS(field):
	with open('rms.txt') as r: 
		next(r)
		for line in r: 
			t = line.split()
			if t[0] == field:
				rms = t[1]
	#print rms
	#rms *= 2
	return rms
#############################
# Init
#############################


cmd1 = "export STARLINK_DIR=/home/emma/star-2017A;"# echo $STARLINK_DIR;"
#cmd2 = "source $STARLINK_DIR/etc/profile;"
cmd2 = "source $STARLINK_DIR/etc/profile;"

convert = '/home/emma/star-2017A/bin/convert/convert.sh'
cupid = '/home/emma/star-2017A/bin/cupid/cupid.sh'

#init = cmd1+cmd2+cmd4

# Individual routines
fits2ndf = '/home/emma/star-2017A/bin/convert/fits2ndf'
ndf2fits = '/home/emma/star-2017A/bin/convert/fits2ndf'
#findclumps = '/home/emma/star-2017A/bin/cupid/findclumps'

smurf      =  "/home/emma/star-2017A/bin/smurf/smurf.csh"
#convert    =  "/home/emma/star-2017A/bin/convert/convert.csh"
#cupid      =  "/home/emma/star-2017A/bin/cupid/cupid.csh"
kappa      =  "/home/emma/star-2017A/bin/kappa/kappa.csh"
makemap       =  "/home/emma/star-2017A/bin/smurf/makemap"
makesnr       =  "/home/emma/star-2017A/bin/kappa/makesnr"
findclumps    =  "/home/emma/star-2017A/bin/cupid/findclumps"
extractclumps =  "/home/emma/star-2017A/bin/cupid/extractclumps"
thresh        =  "/home/emma/star-2017A/bin/kappa/thresh"
#fits2ndf      =  "/home/emma/star-2017A/bin/convert/fits2ndf"
#ndf2fits      =  "/home/emma/star-2017A/bin/convert/ndf2fits"
findback      =  "/home/emma/star-2017A/bin/cupid/findback"
stats         =  "/home/emma/star-2017A/bin/kappa/stats"
ndfcopy       =  "/home/emma/star-2017A/bin/kappa/ndfcopy"
cmd3 = convert
cmd4 = cupid
init = cmd1+cmd2+cmd3+';'+cmd4+';'


#############################
# Write RMS values to file
#############################
#rmsVal = open('rms.txt','w')
#rmsVal.write('SCUBA			rms\n')

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
	#	print name

		ra = p[1]
		dec = p[2]

	#	rms = p[3]
		rms_small = getRMS(p[0])
		rms = str((float(rms_small))*2.0)

	#	if field != 'G001_36+20_96': #and field != 'G005_93-01_00':	
	#		continue
		
#############################
# get_RMS.py
#############################
	#	center_deg = (ra_to_deg(ra),dec_to_deg(dec))
	#
	#	ra0,dec0 = center_deg
#
#		print ra,' ',dec
	#	print center_deg



	#	scubaFile_3D = scubaFits_in+scuba
	#	f = fits.open(scubaFile_3D)#[0]
	#	sc = FixFits(f)   # sc is same type as f

		

	#	fData = sc[0].data

		
	#	fig = aplpy.FITSFigure(sc)
	#	I,J = np.indices(fData.shape)
	#	RA,D = fig.pixel2world(J,I)



	#	dist = distance_on_sphere_deg(ra0,dec0,RA,D)
	#	m = np.nonzero((dist<0.02)&(np.isfinite(fData)))

	#	rms = getRMS(p[0])#,ra,dec)#np.std(fData[m])

		
		
		#print e
		#print 5*'*'
	#	rmsVal.write(p[0]+' '+str(rms)+'\n')

	#	rms = float(rms)

	#	print type(rms)
#############################	
	#	center_deg = (ra_to_deg(ra),dec_to_deg(dec))
		
	#	ra0,dec0 = center_deg

	#	print ra,' ',dec
	#	print center_deg



		

	#	
	#	f = fits.open(scubaFits_in+scuba)[0]
	#	fData = f.data
	#	fD = fData.reshape(fData.shape[1],fData.shape[2])
	#	RA,D = pixel_to_wcs(scubaFits_in+scuba)

	#	dist = distance_on_sphere_deg(ra0,dec0,RA,D)
	#	m = np.nonzero((dist<0.02)&(np.isfinite(f.data)))
#
	#	e = np.std(f.data[m])
	#	print e
	#	print 5*'*'


		
#############################
# Do FellWalker
#############################
		clumps_out = outName(name)

		#print init+findclumps+' '+in_850+name+' '+out_850+clumps_out+' '+'FellWalker'
		
		name_in = in_850+name
		name_out = out_850+clumps_out

		print name_in,'\n',name_out,'\n',rms

		cat = out_850+name.strip('.sdf')+'_cat'

		command = init+findclumps+' '+name_in+' outcat= '+cat+' out='+name_out+' method=FellWalker'+' RMS='+rms+' deconv=false config=^config_file_original.txt'
		#command = 'cd '+in_850+';'+init+findclumps+' '+name_in+' '+name_out+' '+'FellWalker'
#		print command
		os.system(command)

		

		#os.system(cmd1+cmd2+cmd3+';'+cmd4+';'+findclumps+' '+in_850+name+' '+out_850+clumps_out+' '+'FellWalker')

#############################
# sdf2fits
#############################



















# config = config_file.txt # fellWalker commands 
