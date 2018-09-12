


##############################################################
#		Draw YSO sources onto SCUBA maps
##############################################################
#			Akari sources
##############################################################

import os

import aplpy
import matplotlib.pyplot as plt
#####################################
# Paths
#####################################
dir_in = '/home/emma/gradu/data/All_data/catalogs/Akari/'

scuba_in = '/home/emma/gradu/data/All_data/SCUBA_850/'

outDir = '/home/emma/gradu/data/All_data/YSO_maps/'

#####################################
# Functions
#####################################
def ra_to_deg(h,m,s):
	h += (m/60.0)
	h += (s/3600.0)
	deg = h*(360.0/24.0)
	return deg

def dec_to_deg(deg,amin,asec):
	if deg < 0:
		deg -= (amin/60.0)
		deg -= (asec/3600.0)
		return round(deg,6)
	deg += (amin/60.0)
	deg += (asec/3600.0)
	return deg

def outName(field):
	out = field+'_akari_spire.png'
	return out


#####################################
# Run code
#####################################

for filename in os.listdir(dir_in):
	if filename.endswith('FIS.txt'):
		if filename == 'table.txt':
			continue
		field = filename.split('_')[0]
		#print field

		scuba = field+'_850.fits'
		f = aplpy.FITSFigure(scuba_in+scuba)


		f.show_colorscale(cmap='gist_earth')
		print 10*'*',field,10*'*'

		with open(dir_in+filename) as coList: 
			for i in range(0,7):
				next(coList)
			for line in coList:
				p = line.split()
				# Already in degrees
				ra = float(p[2])
				dec = float(p[3])
				if field == 'G016.37-00.61':
					print ra,'	',dec
		
				try: 
					f.add_label(ra,dec,'*',color='white',size=35)
				except:
					print 'not possible'
		pscFile = field+'_PSC.txt'
		#print pscFile
		
		with open(dir_in+pscFile) as psc:
			for i in range(0,7):
				next(psc)
			for line in psc:
				p = line.split()
				# Already in degrees
				ra = float(p[2])
				dec = float(p[3])
				
				try: 
					f.add_label(ra,dec,'*',color='black',size=25)
				except:
					print 'not possible'




		name_out = outName(field)
		f.save(outDir+name_out)
		f.close()




