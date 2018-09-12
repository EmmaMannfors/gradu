#####################################################################################
#		Draw SCUBA-2 and PACS 70 um images next to each other
#####################################################################################
#				Jul 16, 2018
#####################################################################################


import os
import aplpy
import matplotlib.pyplot as plt

##############################################
#	Paths
##############################################
pacs_in = '/home/emma/gradu/data/All_data/PACS_70_Herschel/cropped_70/'
scuba_in = '/home/emma/gradu/data/All_data/SCUBA_850/'
dir_out = '/home/emma/gradu/codes/draw_pictures_for_gradu/'

##############################################
#	Code
##############################################
scubaName = 'G001.36+20.96_850.fits'
pacsName = 'G001.15+21.18_070_L2.5_J_crop.fits'

fig = plt.figure()

img1 = scuba_in+scubaName
img2 = pacs_in+pacsName

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


#f1.add_colorbar()
#f1.colorbar.set_width(0.1)
#f1.colorbar.set_location('top')
#f1.colorbar.set_font(size='xx-small')



#f2.add_colorbar()
#f2.colorbar.set_width(0.1)
#f2.colorbar.set_location('top')
#f2.colorbar.set_font(size='x-small')

fig.canvas.draw()

plt.savefig(dir_out+'70_850um.png',bbox_inches='tight')

