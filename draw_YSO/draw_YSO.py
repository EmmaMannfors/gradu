##############################################################
#		Draw YSO sources onto maps
##############################################################
#		onto SCUBA maps? NH2 maps????
##############################################################
#			SIMBAD / OLD VERSION
##############################################################


OLD VERSION!!!!


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

n_tau_in = '/home/emma/gradu/data/All_data/combinations/A_plus_B_min_C/'

n_tau__all_in = '/home/emma/gradu/data/All_data/combinations/PACS_SPIRE/removedBadPixels/'
matches = '/home/emma/gradu/codes/matches/cropMatches.txt'
#####################################
# Functions
#####################################
def ra_to_deg(h,m,s):
	#g = ra.split(':')
	#h = float(g[0])
	#m = float(g[1])
	#s = float(g[2])
	h += (m/60.0)
	h += (s/3600.0)
	deg = h*(360.0/24.0)
	return deg

def dec_to_deg(deg,amin,asec):
	#y = dec.split(':')
	#deg,amin,asec = float(y[0]),float(y[1]),float(y[2])
	if deg < 0:
		deg -= (amin/60.0)
		deg -= (asec/3600.0)
		return round(deg,6)
	deg += (amin/60.0)
	deg += (asec/3600.0)
	return deg

def outName(field):
	out = field+'_n_tau_SPIRE.png'
	return out
def outName_2(field):
	out = field+'_n_tau_SPIRE_2.png'
	return out
#####################################
# Run code
#####################################

for filename in os.listdir(list_in):
	if filename.endswith('.txt'):
		field = filename.strip('.txt')
		print field
		
		# For plotting on scuba maps
		#scuba = field+'_850.fits'
		#f = aplpy.FITSFigure(scuba_in+scuba)#,figure=fig)

		# For plotting on N_tau maps
		scuba = field+'_850.fits'
		psw_now = 'hh'
		with open(matches) as m:
			next(m)
			for line in m: 
				e = line.split()
				now = e[0]
				if now == scuba:
					psw_now = e[6]
					break
		s = psw_now.split('_')
		if s[0] == 'NA' or s[0] == 'hh':
			continue
		print s
		#psw = s[0]+'_AE_zeros_1.fits'
		if s[3] == 'crop.fits':
			psw = s[0]+'_ABC_1.fits'
			name_out = outName(field)
		elif s[4] == '2.fits':
			psw = s[0]+'_ABC_2.fits'
			name_out = outName_2(field)
		else: 
			print 'AAAAAA\n',s[3],'	',s[4]
			continue
		print psw
		try: 
			f = aplpy.FITSFigure(n_tau_in+psw)	
		except: 
			print 'map does not exist'
			continue
	#	fig = plt.figure()
	#	if filename != 'G171.53-14.91.txt':	# 2 YSO, 0 candidates
	#	if filename != 'G202.31+02.53.txt':	# 18 YSOs + candidates
	#		continue


		
		f.show_colorscale(cmap='default') #'jet_r', vmin=0,vmax=(number)
		# 'color_r' flips the direction of the colorbar around
		#f.save('test.png')


		with open(list_in+filename) as coList:
			for i in range(0,9):
				next(coList)
			for line in coList:
				if line[0] == '=':
					break
				p = line.split()


				# Only collecting YSOs
				starType = p[4]
				if starType != 'Y*O' and starType != 'Y*?':
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
				if starType == 'Y*?': # Plot YSO candidates
					try:
						#f.show_circles(ra, dec, 0.002,color='k')
						f.show_markers(ra, dec,color='y')#'ro'

#add_label(ra,dec,'*',color='black',size=13)
						f.show_rectangles(ra,dec, 0.004, 0.004)
					except:
						print '\n'+'Unable to draw plot\n'+starType+'	'+str(ra)+'	'+str(dec)+'\n'+10*'*'

				else:  # Plot YSOs

					try:
						f.show_circles(ra, dec, 0.003)#, layer=False, zorder=None, **kwargs)
						#f.show_markers(ra, dec,'co')
					except:
						print '\n'+'Unable to draw plot\n'+str(starType)+'	'+str(ra)+'	'+str(dec)+'\n'
		f.save(outDir+name_out)	
		f.close()	
	#	plt.close()
		
		



# YSO = Y*O 
# YSO candidate = Y*?

# T tauri = TT*
# TT candidate = TT?

# pre-main sequence = pr*
# PMS candidate = pr?


#s_a = subplot_area(band)
#		loc = path(band)
#		f = aplpy.FITSFigure(loc+filename,figure=fig,subplot=s_a)
#		currentPath = path(band)
#		f.set_tick_labels_font(size='x-small')
#		f.set_axis_labels_font(size='small')
#		#show y-axis labels and ticks
#		f.axis_labels.hide_y()
#		f.tick_labels.hide_y()
#		f.axis_labels.hide_x()
#		f.tick_labels.hide_x()
#		f.tick_labels.set_xformat('ddd.d')#'hh:mm:ss.ss')
#	#	f.tick_labels.set_xposition('top')
#	#	f.tick_labels.set_xformat('hh:mm:ss')#'hh:mm:ss.ss')
#		if band == 70 or band == 250:
#			f.tick_labels.show_y()
#			f.axis_labels.show_y()
#		if band == 350 or band == 250 or band == 500:
#			f.tick_labels.show_x()
#			f.axis_labels.show_x()
#		#f.add_label(0.1, 0.9, '(a)', relative=True)
#		f.show_colorscale(cmap='gist_heat')
#	#	f.add_colorbar()
#		f.add_label(0.25, 0.93, str(band)+' um', relative=True)
#		f.add_label(0.45,0.05,filename.strip(currentPath).split('_')[0],relative=True)
#		fig.canvas.draw()
	


