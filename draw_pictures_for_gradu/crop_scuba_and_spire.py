#####################################################################################
#		Draw the SCUBA and its corresponding SPIRE image
#		Along with large Herschel image if possible
#####################################################################################
#				Jul 16, 2018
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

dir_out = '/home/emma/gradu/codes/draw_pictures_for_gradu/'

scubaName = 'G006.01+36.74_850.fits'
pswName = 'G006.04+36.73_PSW_L2.0_crop.fits'
bigName = 'G006.04+36.73_PSW_L2.0.fits'

fig = plt.figure()

fb = aplpy.FITSFigure(big_in+bigName,figure=fig,subplot=[0.0,0.0,0.8,0.8])
fb.show_colorscale(cmap='inferno',vmin=0)
fb.tick_labels.hide_y()
fb.tick_labels.hide_x()
fb.axis_labels.hide_x()
fb.axis_labels.hide_y()
fb.ticks.hide()
fb.show_rectangles(238.536641667, -2.87842222222, 0.2, 0.2,color='w')
fb.frame.set_color('w')

fp = aplpy.FITSFigure(psw_in+pswName,figure=fig,subplot=[0.6,0.6,0.2,0.27])
fp.show_colorscale(cmap='inferno',vmin=0)
fp.tick_labels.hide_y()
fp.tick_labels.hide_x()
fp.axis_labels.hide_x()
fp.axis_labels.hide_y()
fp.ticks.hide()


fs = aplpy.FITSFigure(scuba_in+scubaName,figure=fig,subplot=[0.6,0.02,0.2,0.27])
fs.show_colorscale(cmap='inferno',vmin=0)
fs.tick_labels.hide_y()
fs.tick_labels.hide_x()
fs.axis_labels.hide_x()
fs.axis_labels.hide_y()
fs.ticks.hide()


fig.canvas.draw()
plt.savefig(dir_out+'crop.png',bbox_inches='tight')














