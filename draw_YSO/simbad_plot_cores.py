##############################################################
#		Draw cores onto maps
##############################################################
#		(symbol: cor) = dense core
##############################################################
#			SIMBAD 
##############################################################

import os
#import montage_wrapper as montage
#from astropy.io import fits

import aplpy
import matplotlib.pyplot as plt
#####################################
# Paths
#####################################
scuba_in = '/home/emma/gradu/data/All_data/SCUBA_850/'

list_in = '/home/emma/gradu/data/All_data/catalogs/simbad/'

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
	out = field+'_cores_scuba.png'
	return out



#####################################
# Run code
#####################################

for filename in os.listdir(list_in):
	if filename.endswith('.txt'):
		field = filename.strip('.txt')
		print field
		
		# For plotting on scuba maps
		scuba = field+'_850.fits'

	#	if filename != 'G089.66-06.62.txt':	# several candidates, 1 YSO, 1 TT
	#	if filename != 'G202.31+02.53.txt':	# 18 YSOs + candidates
	#		continue

		f = aplpy.FITSFigure(scuba_in+scuba)#,figure=fig)


		

		
		f.show_colorscale(cmap='magma',vmin=0) #'jet_r', vmin=0,vmax=(number)
		# 'color_r' flips the direction of the colorbar around
		#f.save('test.png')


		with open(list_in+filename) as coList:
			for i in range(0,9):
				next(coList)
			for line in coList:
				if line[0] == '=':
					break
				p = line.split()


				# Only collecting cores
				starType = p[4]
				if starType != 'cor':
					continue
				#elif starType == 'Y*?':  # CURRENTLY just draw known YSOs
					
				#	continue
				#print starType
		
				# Coordinates
				h,m,s = float(p[5]),float(p[6]),float(p[7])
				deg,amin,asec = float(p[8]),float(p[9]),float(p[10])
				ra_s = p[5]+':'+p[6]+':'+p[7]
				dec_s = p[8]+':'+p[9]+':'+p[10]
				ra = ra_to_deg(h,m,s)
				dec = dec_to_deg(deg,amin,asec)
				#print ra
				#print dec
				if starType == 'cor': # Plot cores
					try:
						f.show_circles(ra, dec, 0.005,color='white')
						#f.add_label(ra,dec,'',color='white',size=10)	
					except:
						print '\n'+'Unable to draw plot\n'+starType+'	'+str(ra)+'	'+str(dec)+'\n'+10*'*'
				

		name_out = outName(field)
		f.save(outDir+name_out)	
		f.close()	
	#	plt.close()
		
		



# YSO = Y*O 
# YSO candidate = Y*?

# T tauri = TT*
# TT candidate = TT?

# pre-main sequence = pr*
# PMS candidate = pr?


