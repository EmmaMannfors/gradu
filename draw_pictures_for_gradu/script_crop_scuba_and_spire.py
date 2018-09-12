#####################################################################################
#		Draw SCUBA location on large Herschel image. 
#####################################################################################
#				Jul 19, 2018
#####################################################################################

import os
import aplpy
import matplotlib.pyplot as plt


##############################################
#	Paths
##############################################
scuba_in = '/home/emma/gradu/data/All_data/SCUBA_850/'
psw_in = '/home/emma/gradu/data/All_data/SPIRE_PSW_250_Herschel/cropped_250/'
big_in = '/home/emma/gradu/data/All_data/SPIRE_PSW_250_Herschel/'


dir_out = '/home/emma/gradu/codes/draw_pictures_for_gradu/all_Herschel_crops/'


doublesList = ['G017.69-00.15_PSW_L2.5.fits','G070.07-01.60_PSW_L2.5.fits','G092.12+03.88_PSW_L2.0.fits']

##############################################
# Plot fields with one SCUBA match
##############################################
#with open('centers_250_850.txt') as c: 
#	next(c)
#	for line in c: 
#		p = line.split() 
#		psw = p[1]
#		if psw == 'NA':
#			continue
#		e = psw.split('_')
#		herschel = e[0]+'_'+e[1]+'_'+e[2]+'.fits'
#		ra = float(p[2])
#		dec = float(p[3])
#
#		if herschel in doublesList: 
#			continue
#
#		#print herschel
#		try: 
#			outName = e[0]+'_'+e[1]+'_'+e[2]+'_'+e[4].strip('.fits')+'.png'
#		except: 
#			outName = e[0]+'_'+e[1]+'_'+e[2]+'_1.png'
#
#		#print outName
#		fig = plt.figure()
#		fb = aplpy.FITSFigure(big_in+herschel,figure=fig,subplot=[0.0,0.0,0.8,0.8])
#		fb.show_colorscale(cmap='inferno',vmin=0)
#		fb.tick_labels.hide_y()
#		fb.tick_labels.hide_x()
#		fb.axis_labels.hide_x()
#		fb.axis_labels.hide_y()
#		fb.ticks.hide()		
#		fb.frame.set_color('w')
#
#		fb.show_rectangles(ra, dec, 0.2, 0.2,color='w')
#
#		f1 = aplpy.FITSFigure(psw_in+psw,figure=fig,subplot=[0.6,0.6,0.2,0.27])
#		f1.show_colorscale(cmap='inferno',vmin=0)
#		f1.tick_labels.hide_y()
#		f1.tick_labels.hide_x()
#		f1.axis_labels.hide_x()
#		f1.axis_labels.hide_y()
#		f1.ticks.hide()
#
#
#		fig.canvas.draw()
#		plt.savefig(dir_out+outName,bbox_inches='tight')
#		plt.close()

##############################################
# Plot fields with two/three SCUBA matches
##############################################
def plotMany(herschel,psw1,psw2,psw3,ra1,dec1,ra2,dec2,ra3,dec3,outName):
	fig = plt.figure()

	fb = aplpy.FITSFigure(big_in+herschel,figure=fig,subplot=[0.0,0.2,0.65,0.65])
	fb.show_colorscale(cmap='inferno',vmin=0)
	fb.tick_labels.hide_y()
	fb.tick_labels.hide_x()
	fb.axis_labels.hide_x()
	fb.axis_labels.hide_y()
	fb.ticks.hide()		
	fb.frame.set_color('w')
	
	fb.show_rectangles(ra1, dec1, 0.2, 0.2,color='w')
	fb.show_rectangles(ra2, dec2, 0.2, 0.2,color='w')
	if ra3 != 0:
		fb.show_rectangles(ra3, dec3, 0.2, 0.2,color='w')
	
	
	f1 = aplpy.FITSFigure(psw_in+psw1,figure=fig,subplot=[0.6,0.6,0.2,0.27])
	f1.show_colorscale(cmap='inferno',vmin=0)
	f1.tick_labels.hide_y()
	f1.tick_labels.hide_x()
	f1.axis_labels.hide_x()
	f1.axis_labels.hide_y()
	f1.ticks.hide()
	
	f2 = aplpy.FITSFigure(psw_in+psw2,figure=fig,subplot=[0.6,0.3,0.2,0.27])
	f2.show_colorscale(cmap='inferno',vmin=0)
	f2.tick_labels.hide_y()
	f2.tick_labels.hide_x()
	f2.axis_labels.hide_x()
	f2.axis_labels.hide_y()
	f2.ticks.hide()
	
	
	if psw3 != 'NA':
		f3 = aplpy.FITSFigure(psw_in+psw3,figure=fig,subplot=[0.6,0.02,0.2,0.27])
		f3.show_colorscale(cmap='inferno',vmin=0)
		f3.tick_labels.hide_y()
		f3.tick_labels.hide_x()
		f3.axis_labels.hide_x()
		f3.axis_labels.hide_y()
		f3.ticks.hide()
	
	
	fig.canvas.draw()
	plt.savefig(dir_out+outName,bbox_inches='tight')
	plt.close()



herschel = 'G017.69-00.15_PSW_L2.5.fits'
psw1 = 'G017.69-00.15_PSW_L2.5_crop.fits'
ra1,dec1 = 275.694470833, -14.9862636111
psw2 = 'G017.69-00.15_PSW_L2.5_crop_2.fits'
ra2,dec2 = 275.177591667, -14.0336327778
psw3 = 'G017.69-00.15_PSW_L2.5_crop_3.fits'
ra3,dec3 = 276.889741667,-14.6250719444
outName = 'G017.69-00.15_PSW_L2.5_all.png'
plotMany(herschel,psw1,psw2,psw3,ra1,dec1,ra2,dec2,ra3,dec3,outName)


herschel = 'G092.12+03.88_PSW_L2.0.fits'
psw1 = 'G092.12+03.88_PSW_L2.0_crop.fits'
ra1,dec1 = 316.033354167,52.5685483333
psw2 = 'G092.12+03.88_PSW_L2.0_crop_2.fits'
ra2,dec2 = 315.6164,52.4752347222
psw3 = 'G092.12+03.88_PSW_L2.0_crop_3.fits'
ra3,dec3 = 315.125583333,52.5216755556
outName = 'G092.12+03.88_PSW_L2.0_all.png'
plotMany(herschel,psw1,psw2,psw3,ra1,dec1,ra2,dec2,ra3,dec3,outName)


herschel = 'G070.07-01.60_PSW_L2.5.fits'
psw1 = 'G070.07-01.60_PSW_L2.5_crop.fits'
ra1,dec1 = 303.512608333,32.0154938889
psw2 = 'G070.07-01.60_PSW_L2.5_crop_2.fits'
ra2,dec2 = 303.393391667,31.3651077778
psw3 = 'NA'
ra3,dec3 = 0,0
outName = 'G070.07-01.60_PSW_L2.5_all.png'
plotMany(herschel,psw1,psw2,psw3,ra1,dec1,ra2,dec2,ra3,dec3,outName)















