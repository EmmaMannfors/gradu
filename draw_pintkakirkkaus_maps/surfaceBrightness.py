#####################################################################################
#			Program to draw surface brightness maps
#			of all Herschel channels
#####################################################################################
#				June 20, 2018
#####################################################################################
import os
import aplpy
import matplotlib.pyplot as plt


##############################################
#	Paths
##############################################
in_70 = '/home/emma/gradu/data/All_data/PACS_70_Herschel/MJy_cropped_70/'
in_100 = '/home/emma/gradu/data/All_data/PACS_100_Herschel/MJy_cropped_100/'
#in_160 = '/home/emma/gradu/data/All_data/PACS_160_Herschel/fixed_units_cropped_160/'
in_160 = '/home/emma/gradu/data/All_data/PACS_160_Herschel/cropped_160_MJy/'
in_250 = '/home/emma/gradu/data/All_data/SPIRE_PSW_250_Herschel/cropped_250/'
in_350 = '/home/emma/gradu/data/All_data/SPIRE_PMW_350_Herschel/cropped_350/'
in_500 = '/home/emma/gradu/data/All_data/SPIRE_PLW_500_Herschel/cropped_500/'
in_850 = '/home/emma/gradu/data/All_data/SCUBA_850/'

#out_test = '/home/emma/gradu/data/All_data/surface_brightness_maps/scriptTest/'
#outTest = '/home/emma/gradu/data/All_data/surface_brightness_maps/scriptTest_2/'
outTest = '/home/emma/gradu/data/All_data/PICTURES/surface_brightness_maps/newColor/'
#out_70 = '/'
#out_100 = '/'
#out_160 = '/'
#out_250 = '/'
#out_350 = '/'
#out_500 = '/'
#out_850 = ''

matchPath = '/home/emma/gradu/codes/matches/'
##############################################
#	Plot function
##############################################
#fig = plt.figure()

def subplot_area(band):
	if band == 70:
		subplot = [0.1,0.5,0.20,0.35]
	elif band == 100:
		subplot = [0.35,0.5,0.20,0.35]
	elif band == 160:
		subplot = [0.6,0.5,0.20,0.35]
	elif band == 250:
		subplot = [0.1,0.1,0.20,0.35]
	elif band == 350:
		subplot = [0.35,0.1,0.20,0.35]
	elif band == 500:
		subplot = [0.6,0.1,0.20,0.35]
	return subplot

def SCUBA_subplot_area(band):
	if band == 70:
		subplot = [0.05,0.5,0.2,0.35]
	elif band == 100:
		subplot = [0.3,0.5,0.2,0.35]
	elif band == 160:
		subplot = [0.55,0.5,0.2,0.35]
	elif band == 250:
		subplot = [0.05,0.1,0.2,0.35]
	elif band == 350:
		subplot = [0.3,0.1,0.2,0.35]
	elif band == 500:
		subplot = [0.55,0.1,0.2,0.35]
	elif band == 850:
		subplot = [0.8,0.3,0.2,0.35]
	return subplot


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

def outName(filename,band):
	path = path(filename,band)
	outname = out_test+filename.strip(path).strip('_crop.fits')+'_test.fits'
	return outname

def plotMap(filename,band):
	if filename == 'NA':
		print 'meow'
	else:	
		s_a = subplot_area(band)
		loc = path(band)
		f = aplpy.FITSFigure(loc+filename,figure=fig,subplot=s_a)
		currentPath = path(band)
		f.set_tick_labels_font(size='x-small')
		f.set_axis_labels_font(size='small')
		#show y-axis labels and ticks
		f.axis_labels.hide_y()
		f.tick_labels.hide_y()
		f.axis_labels.hide_x()
		f.tick_labels.hide_x()
		f.tick_labels.set_xformat('ddd.d')#'hh:mm:ss.ss')
	#	f.tick_labels.set_xposition('top')
	#	f.tick_labels.set_xformat('hh:mm:ss')#'hh:mm:ss.ss')
		if band == 70 or band == 250:
			f.tick_labels.show_y()
			f.axis_labels.show_y()
		if band == 350 or band == 250 or band == 500:
			f.tick_labels.show_x()
			f.axis_labels.show_x()
		#f.add_label(0.1, 0.9, '(a)', relative=True)
		f.show_colorscale(cmap='gist_heat')
	#	f.add_colorbar()
		f.add_label(0.25, 0.93, str(band)+' um', relative=True)
		f.add_label(0.45,0.05,filename.strip(currentPath).split('_')[0],relative=True)
		fig.canvas.draw()
	

def SCUBA_plotMap(filename,band):
	if filename == 'NA':
		print 'meow'

	else:	
		s_a = SCUBA_subplot_area(band)
		loc = path(band)
		f = aplpy.FITSFigure(loc+filename,figure=fig,subplot=s_a)
		currentPath = path(band)

		f.set_tick_labels_font(size='x-small')
		f.set_axis_labels_font(size='small')
		f.axis_labels.hide_y()
		f.tick_labels.hide_y()
		f.axis_labels.hide_x()
		f.tick_labels.hide_x()
		f.tick_labels.set_xformat('ddd.d')
		

		if band == 70 or band == 250:
			f.tick_labels.show_y()
			f.axis_labels.show_y()
		if band == 350 or band == 250 or band == 500:
			f.tick_labels.show_x()
			f.axis_labels.show_x()

	#	f.show_colorscale(cmap='gist_heat')

		f.show_colorscale(cmap='gist_earth')

		f.add_colorbar()
		f.colorbar.set_width(0.1)
		f.colorbar.set_location('top')

		f.colorbar.set_font(size='x-small')
		f.colorbar.set_axis_label_text('Flux (MJy/beam)')
	#	if band == 850:	
	#		f.add_colorbar()
	#		f.colorbar.set_location('right')

		f.add_label(0.25, 0.93, str(band)+' um', relative=True)
		if band == 850:	
			f.add_label(0.45,0.05,filename.strip(currentPath).split('_')[0],relative=True)
		#	f.ticks.show_y()
		#	f.ticks.show_x()
			f.ticks.set_color('black')
		fig.canvas.draw()	

############################
# Code
############################
with open(matchPath+'cropMatches.txt') as matches:
	next(matches)
	for line in matches:
		fig = plt.figure()
		p = line.split()
		
	#	print ' name = ',change(str(p[1]))
	#	print p[1]
		SCUBA_plotMap(p[0],850)
		SCUBA_plotMap(change(p[1]),70)
		SCUBA_plotMap(change(p[2]),100)
		SCUBA_plotMap(change(p[3]),160)
		SCUBA_plotMap(p[4],250)
		SCUBA_plotMap(p[5],350)
		SCUBA_plotMap(p[6],500)
	
		outName = p[0].split('_')[0]+'.png'
		print outName

		plt.savefig(outTest+outName,bbox_inches='tight')
		plt.close()



############################
# TEST
############################
#with open('testFiles.txt') as test:
#	for line in test:
#		fig = plt.figure()
#		p = line.split()
#
#		scuba = 'G159.23-20.09_850.fits'
#		SCUBA_plotMap(scuba,850)
#		SCUBA_plotMap(p[0],70)
#		SCUBA_plotMap(p[1],100)
#		SCUBA_plotMap(p[2],160)
#		SCUBA_plotMap(p[3],250)
#		SCUBA_plotMap(p[4],350)
#		SCUBA_plotMap(p[5],500)
#		
#		title1 = p[0].split('_')[0]
#		title2 = p[3].split('_')[0]
#		#plt.title(title)
#	#	plt.text(1, 145, title1)#120, title)
#	#	plt.text(1,65,title2)
#		plt.savefig(outTest+p[1].split('_')[0]+'.png',bbox_inches='tight')








