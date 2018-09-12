###############################################################################
#		 Plot found YSOs on SCUBA-2 images
###############################################################################
#				jul 30, 2018
###############################################################################


import os
#import montage_wrapper as montage
#from astropy.io import fits

import aplpy
import matplotlib.pyplot as plt

#####################################
# Paths
#####################################
scuba_in = '/home/emma/gradu/data/All_data/SCUBA_850/'
clumps_in = '/home/emma/gradu/data/All_data/SCUBA_850/clumps/fits_clumps/'

list_in = '/home/emma/gradu/data/All_data/YSO/clump_and_YSO_list/no_duplicates/'

outDir = '/home/emma/gradu/data/All_data/YSO/clump_and_YSO_list/plotted_YSOs/'




#####################################
# Functions
#####################################

def outName(field):
	out = field+'_associated_YSO.png'
	return out

def outName_cl(field):
	out = field+'_clumps_and_YSO.png'
	return out


#####################################
# Run code
#####################################

for filename in os.listdir(list_in):
	if filename.endswith('.txt'):
		field = filename.strip('_no_duplicates.txt')
		e = field.split('.')
		field_s = e[0]+'_'+e[1]+'_'+e[2]
		#print field
		
		# For plotting on scuba maps
		scuba = field+'_850.fits'
		cl = field_s+'_850_clump.fits'

	#	if filename != 'G089.66-06.62.txt':	# several candidates, 1 YSO, 1 TT
	#	if field != 'G202.31+02.53':	# 18 YSOs + candidates
	#	if field != 'G215.87-17.50':	# 18 YSOs + candidates
	#		continue

	#	if field != 'G026.54+00.72':
	#		continue
		f = aplpy.FITSFigure(scuba_in+scuba)#,figure=fig)

		f_clump = aplpy.FITSFigure(clumps_in+cl)

		

		
		f.show_colorscale(cmap='magma',vmin=0) #'jet_r', vmin=0,vmax=(number)
		f_clump.show_colorscale(cmap='rainbow',vmin=0)	

		# 'color_r' flips the direction of the colorbar around
		#f.save('test.png')


		with open(list_in+filename) as coList:
			for line in coList:
				if line[0] == '#':
					continue
				p = line.split()


				name = p[0]

				
					
				
				
		
				ra = float(p[1])
				dec = float(p[2])
				clump = p[3]
				if clump == 'nan':
					continue
				#print ra
				#print dec

				print clump
				
				
			
				if field == 'G026.54+00.72':
					f_clump.show_circles(ra, dec, 0.005,facecolor='midnightblue',edgecolor='white')
					f.show_circles(ra, dec, 0.005,facecolor='white')
					f.show_circles(ra, dec, 0.005,color='black')
				else:
					f.show_circles(ra, dec, 0.001,color='black')
					f.show_circles(ra, dec, 0.0015,facecolor='white')
					f_clump.show_circles(ra, dec, 0.0015,facecolor='midnightblue',edgecolor='white')

				#f.add_label(ra,dec,'*',color='white',size=20)				
				#f.add_label(ra,dec,'',color='black',size=15)
				#f.show_circles(ra, dec, 0.0023,facecolor='white',alpha=0.5)#,alpha=0.5)
				#f_clump.add_label(ra,dec,'*',color='midnightblue',size=25)
				#if starType == 'cor': # Plot cores
				#	try:
				#f.show_circles(ra, dec, 0.005,color='white')
						#f.add_label(ra,dec,'',color='white',size=10)	
				#	except:
				#		print '\n'+'Unable to draw plot\n'+starType+'	'+str(ra)+'	'+str(dec)+'\n'+10*'*'
				

		name_out = outName(field)
		out = outName_cl(field)
		f.save(outDir+name_out)	
		f.close()
		f_clump.save(outDir+out)	
		f_clump.close()	
	#	plt.close()
		
"""		
"""
