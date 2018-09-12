##############################################################
#		Draw YSO sources onto maps
##############################################################
#		onto SCUBA maps? NH2 maps????
##############################################################
#			Vizier sources
##############################################################

import os

import aplpy
import matplotlib.pyplot as plt
#####################################
# Paths
#####################################
scuba_in = '/home/emma/gradu/data/All_data/SCUBA_850/'

class_I_in = '/home/emma/gradu/data/All_data/catalogs/vizier/class_I_II/'
class_III_in = '/home/emma/gradu/data/All_data/catalogs/vizier/class_III/'
outDir = '/home/emma/gradu/data/All_data/YSO_maps/'


#n_tau_in = '/home/emma/gradu/data/All_data/combinations/A_plus_B_min_C/'
#n_tau__all_in = '/home/emma/gradu/data/All_data/combinations/PACS_SPIRE/removedBadPixels/'
#matches = '/home/emma/gradu/codes/matches/cropMatches.txt'
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
	out = field+'_vizier_spire.png'
	return out

#def outName_2(field):
#	out = field+'_n_tau_SPIRE_2.png'
#	return out
#####################################
# Run code
#####################################

for filename in os.listdir(class_I_in):
	if filename.endswith('.txt'):
		if filename == 'table.txt':
			continue
		field = filename.strip('.txt')
		#print field

		#if filename != 'G171.53-14.91.txt':	# 2 YSO, 0 candidates
		#	continue
		
		# For plotting on scuba maps
		scuba = field+'_850.fits'
		#print scuba
		f = aplpy.FITSFigure(scuba_in+scuba)#,figure=fig)

		
	#	fig = plt.figure()
		


		
		f.show_colorscale(cmap='gist_stern') #'jet_r', vmin=0,vmax=(number)
		# 'color_r' flips the direction of the colorbar around
		
		with open(class_I_in+filename) as coList: 
			next(coList)
			for line in coList:
				if line[0] == '' or line == 'Result truncated to 50 rows on a total of 62 matching rows':
					continue
				p = line.split()

				# Already in degrees
				try: 
					ra = float(p[3])
					dec = float(p[4])
				except:  # For blank lines etc
					print 'This line is blank',line
					continue
				#print ra,'	',dec

				try: 
					f.add_label(ra,dec,'*',color='black',size=25)
				except: 
					print 'unable to draw\n'+scuba+'	'+str(ra)+'	'+str(dec)
		try:
			with open(class_III_in+filename) as coList: 
				next(coList)
				for line in coList:
					if line[0] == '':
						continue
					p = line.split()
	
					# Already in degrees
					try: 
						ra = float(p[3])
						dec = float(p[4])
					except:  # For blank lines etc
						print 'This line is blank',line
						continue
					#print ra,'	',dec
	
					try: 
						f.add_label(ra,dec,'*',color='white',size=35)
					except: 
						print 'unable to draw\n'+scuba+'	'+str(ra)+'	'+str(dec)
		except: 
			print 'file does not have class III protostars'

		name_out = outName(field)
		f.save(outDir+name_out)
		f.close()





#f.save(outDir+name_out)	

		
		













