######################################################################
#		Find YSOs that are associated with clumps
#			Return a list? Map?
######################################################################
#			Jul 25, 2018
######################################################################

import os
from astropy.io import fits
#from astropy import units as u
#from astropy.coordinates import SkyCoord
#from astropy.wcs import WCS
import aplpy
from functions import distance_on_sphere_deg
import numpy as np
from FixFits import FixFits

###########################
# Paths
###########################
clumps_in = '/home/emma/gradu/data/All_data/SCUBA_850/clumps/fits_clumps/'
files_out = '/home/emma/gradu/data/All_data/YSO/clump_and_YSO_list/'

###########################
# YSO list directories
###########################
vizier_I_II_dir = '/home/emma/gradu/data/All_data/catalogs/vizier/class_I_II/'
# G001.36+20.96.txt
vizier_III_dir = '/home/emma/gradu/data/All_data/catalogs/vizier/class_III/'
# G010.21+02.40.txt
simbad_dir = '/home/emma/gradu/data/All_data/catalogs/simbad/'
# G001.36+20.96.txt
Akari_dir = '/home/emma/gradu/data/All_data/catalogs/Akari/'
# G001.36+20.96_FIS.txt
# G001.36+20.96_PSC.txt


###########################
# Functions
###########################
def vI(vizierIFile):
	vizI = open(vizierIFile,'r')
	for line in vizI: 
		v = line.split()
		
def get_txt_name(field):
	return files_out+field+'.txt'

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

#i = 0

###########################
# Missing fields
###########################
#vizierDoc = open('vizier_class_I_II.txt','w')
missList = ['G006.01+36.74','G023.35-00.26','G034.74-01.39','G082.40-01.84','G082.42-01.84','G109.81+02.72','G202.31+02.53']

###########################
# Code
###########################

for filename in os.listdir(clumps_in):
	if filename.endswith('.fits'):
	
	#	if i > 1: 
	#		break
		
	#	if filename != 'G001_36+20_96_850_clump.fits':
	#		continue

		p = filename.split('_')
		field = p[0]+'.'+p[1]+'.'+p[2]
#########
#for field in missList:
#		p = field.split('.')
#		filename = p[0]+'_'+p[1]+'_'+p[2]+'_850_clump.fits'
########

		doc = open(get_txt_name(field),'w')
		
		f = fits.open(clumps_in+filename)
		sc = FixFits(f)
		fData = sc[0].data
		

		#vizierDoc.write('# '+field+'\n')
		doc.write('# '+field+'\n')

		####################
		# catalog names
		####################
		# All but akari
		name = field+'.txt'

		vizierIFile = vizier_I_II_dir+name
		vizierIIIFile = vizier_III_dir+name
		simbadFile = simbad_dir+name
		Akari_FIS_File = Akari_dir+field+'_FIS.txt'
		Akari_PSC_File = Akari_dir+field+'_PSC.txt'
		
		###########################
		# open Vizier I/II catalog
		###########################
		# Try: vizier I/II
		try: 
			vizI = open(vizierIFile,'r')
			doc.write('# Vizier I/II \n# Name			ra		dec		clump\n')
			next(vizI)
			for line in vizI: 
				if len(line.strip()) == 0:
					continue
				v = line.split()
				# v[2] = Name, v[3] = ra, v[4] = dec
				n = v[2]
				ra0 = float(v[3])
				dec0 = float(v[4])
	
				print 'Vizier I/II: '+ra0,dec0
	
				
	
				fig = aplpy.FITSFigure(sc)
				
				# Some function that calculates whether there is a clump within a few pixels of this area	
				# r = 0.005 deg
				I,J = np.indices(fData.shape)
				RA,D = fig.pixel2world(J,I)
				dist = distance_on_sphere_deg(ra0,dec0,RA,D)
				m = np.nonzero((dist<0.005)&(np.isfinite(fData)))
				if fData[m].shape == (0,):
						e = 'nan'
				else: 
					a = np.mean(fData[m])
					e = round(a,0)	
				
				print a
				
				####################
				# Round a to nearest whole number: 
				# Any stars that fall onto
				# two clumps will be assigned to 
				# the clump they share the most area with. 
				####################
				#e = round(a,0)
				
				doc.write(n+'	'+str(ra0)+'	'+str(dec0)+'	'+str(e)+'\n')
	
	
				fig.close()

		except: 
			print 'Vizier I/II: No file exists'
			doc.write('# Vizier I/II \n# No associated YSOs\n')
			
			
		

		###########################
		# open Vizier III catalog
		###########################
		try: 
			vizIII = open(vizierIIIFile,'r')
			doc.write('# Vizier III \n# Name			ra		dec		clump\n')
		
			next(vizIII)
			for line in vizIII:
				if len(line.strip()) == 0:
					continue
				v = line.split()
				# v[2] = Name, v[3] = ra, v[4] = dec
				n = v[2]
				ra0 = float(v[3])
				dec0 = float(v[4])
	
				fig = aplpy.FITSFigure(sc)
				
				I,J = np.indices(fData.shape)
				RA,D = fig.pixel2world(J,I)
				dist = distance_on_sphere_deg(ra0,dec0,RA,D)
				m = np.nonzero((dist<0.005)&(np.isfinite(fData)))
				if fData[m].shape == (0,):
					e = 'nan'
				else: 
					a = np.mean(fData[m])
					e = round(a,0)	
							
				####################
				# Round a to nearest whole number: 
				####################
				#e = round(a,0)
				
				doc.write(n+'	'+str(ra0)+'	'+str(dec0)+'	'+str(e)+'\n')
	
				fig.close()
		except: 
			print 'Vizier III: No file exists'
			doc.write('# Vizier III \n# No associated YSOs\n')
			
		
		###########################
		# open Akari FIS catalog
		###########################
		
		
		try: 
			fis = open(Akari_FIS_File,'r')
			doc.write('# Akari FIS \n# Name		ra		dec		clump\n')

		

			for l in range(0,7):
				next(fis)
	
			for line in fis: 
				if len(line.strip()) == 0:
					continue
				s = line.split()
				
				# s[2] = ra, s[3] = dec
	
				
				ra0 = float(s[2])
				dec0 = float(s[3])
	
				fig = aplpy.FITSFigure(sc)
				
				I,J = np.indices(fData.shape)
				RA,D = fig.pixel2world(J,I)
				dist = distance_on_sphere_deg(ra0,dec0,RA,D)
				m = np.nonzero((dist<0.005)&(np.isfinite(fData)))
				if fData[m].shape == (0,):
					e = 'nan'
				#print fData[m].shape
				else: 
					a = np.mean(fData[m])
					e = round(a,0)	
				####################
				# Round a  
				####################
				print 'Akari FIS'
				
				doc.write('	-	'+str(ra0)+'	'+str(dec0)+'	'+str(e)+'\n')
	
				fig.close()
	
		except: 
			print 'No file exists'
			doc.write('# Akari FIS \n# No associated YSOs\n')
			
		
		###########################
		# open Akari PSC catalog
		###########################
		try: 
			psc = open(Akari_PSC_File,'r')
			doc.write('# Akari PSC \n# Name		ra		dec		clump\n')

			for l in range(0,7):
				next(psc)
	
			for line in psc: 
				if len(line.strip()) == 0:
					continue
				c = line.split()
				
				# c[2] = ra, c[3] = dec
				ra0 = float(c[2])
				dec0 = float(c[3])
	
				fig = aplpy.FITSFigure(sc)
				
				I,J = np.indices(fData.shape)
				RA,D = fig.pixel2world(J,I)
				dist = distance_on_sphere_deg(ra0,dec0,RA,D)
				m = np.nonzero((dist<0.005)&(np.isfinite(fData)))
				if fData[m].shape == (0,):
					e = 'nan'
				else: 
					a = np.mean(fData[m])
					e = round(a,0)	
					
				doc.write('	-	'+str(ra0)+'	'+str(dec0)+'	'+str(e)+'\n')

				print 'Akari PSC'

				fig.close()


		except: 
			print 'No file exists'
			doc.write('# Akari PSC \n# No associated YSOs\n')
			

		###########################
		# open Simbad catalog
		###########################
		try: 
			sim = open(simbadFile,'r')
			doc.write('# Simbad \n# Name		ra		dec		clump	type\n')
			

		
			for l in range(0,9):
				next(sim)
	
			for line in sim: 
				if len(line.strip()) == 0:
					continue
				
				if line[0] == '=':
					continue
				d = line.split()
	
				
				# c[2] = ra, c[3] = dec
				n = d[2]+'_'+d[3]
				typ = d[4]
				if typ != 'TT*' and typ != 'Y*O' and typ != 'Y*?':
					continue
				ra_s = d[5]+':'+d[6]+':'+d[7]
				dec_s = d[8]+':'+d[9]+':'+d[10]
	
	
				print ra_s, dec_s
				ra0 = ra_to_deg(ra_s)
				dec0 = dec_to_deg(dec_s)
#	
				fig = aplpy.FITSFigure(sc)
				
				I,J = np.indices(fData.shape)
				RA,D = fig.pixel2world(J,I)
				dist = distance_on_sphere_deg(ra0,dec0,RA,D)
				m = np.nonzero((dist<0.005)&(np.isfinite(fData)))
				if fData[m].shape == (0,):
					e = 'nan'
				else: 
					a = np.mean(fData[m])
					e = round(a,0)	
					
				doc.write(n+'	'+str(ra0)+'	'+str(dec0)+'	'+str(e)+'	'+typ+'\n')
	
				fig.close()

		except: 
			print 'No file exists'
			doc.write('# Simbad \n# No associated YSOs\n')
			
	#	i += 1
		continue
		
		
#TT*, Y*O,Y*?



			



