#####################################################################################
#		Draw the two crappy images that had to be removed from
#				calculate_N codes
#####################################################################################
#				Jul 16, 2018
#####################################################################################

import os
import aplpy
import matplotlib.pyplot as plt


##############################################
#	Paths
##############################################
dir_in = '/home/emma/gradu/data/All_data/SPIRE_PSW_250_Herschel/cropped_250/'
dir_out = '/home/emma/gradu/codes/draw_pictures_for_gradu/'

##############################################
#	Two files
##############################################
name1 = 'G013.90-00.51_PSW_L2.0_crop.fits'
name2 = 'G035.19-01.75_PSW_L2.0_crop.fits'

img1 = dir_in+name1
img2 = dir_in+name2

fig = plt.figure()


#[x1,y1,dx,dy]

f1 = aplpy.FITSFigure(img1,figure=fig,subplot=[0.1,0.3,0.35,0.4])
f2 = aplpy.FITSFigure(img2,figure=fig,subplot=[0.55,0.3,0.35,0.4])




f1.show_colorscale(cmap='gist_earth',vmin=0)
f2.show_colorscale(cmap='gist_earth',vmin=0)

f1.tick_labels.set_xformat('ddd.d')
f2.tick_labels.set_xformat('ddd.d')

f1.tick_labels.set_yformat('ddd.d')#'h')
f2.tick_labels.set_yformat('ddd.d')

f2.axis_labels.hide_y()

f1.ticks.set_color('black')
f2.ticks.set_color('black')


f1.add_colorbar()
f1.colorbar.set_width(0.1)
f1.colorbar.set_location('top')
f1.colorbar.set_font(size='xx-small')



f2.add_colorbar()
f2.colorbar.set_width(0.1)
f2.colorbar.set_location('top')
f2.colorbar.set_font(size='x-small')

fig.canvas.draw()

plt.savefig(dir_out+'crap_2.png',bbox_inches='tight')


#		f.set_tick_labels_font(size='x-small')
#		f.set_axis_labels_font(size='small')
#		f.axis_labels.hide_y()
#		f.tick_labels.hide_y()
#		f.axis_labels.hide_x()
#		f.tick_labels.hide_x()
#		f.tick_labels.set_xformat('ddd.d')
		

#		if band == 70 or band == 250:
#			f.tick_labels.show_y()
#			f.axis_labels.show_y()
#		if band == 350 or band == 250 or band == 500:
#			f.tick_labels.show_x()
#			f.axis_labels.show_x()

#	#	f.show_colorscale(cmap='gist_heat')
#
#		f.show_colorscale(cmap='default')
#
#		f.add_colorbar()
#		f.colorbar.set_width(0.1)
#		f.colorbar.set_location('top')
#
#		f.colorbar.set_font(size='x-small')
#		f.colorbar.set_axis_label_text('Flux (MJy/beam)')
	#	if band == 850:	
	#		f.add_colorbar()
	#		f.colorbar.set_location('right')

#		f.add_label(0.25, 0.93, str(band)+' um', relative=True)
#		if band == 850:	
#			f.add_label(0.45,0.05,filename.strip(currentPath).split('_')[0],relative=True)
		#	f.ticks.show_y()
		#	f.ticks.show_x()
#			f.ticks.set_color('black')
#		fig.canvas.draw()
