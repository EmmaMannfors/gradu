#################################################################################
#		Perform EXTRACTCLUMPS on FellWalker results
#################################################################################
#				Jul 31, 2018
#################################################################################
# 	Each instance of os.system() makes a new shell, and so Starlink-
#			commands must be run in each excecution
#################################################################################

import os
from astropy.io import fits
import numpy as np
#import aplpy



#############################
# paths
#############################
data_in = '/home/emma/gradu/data/All_data/SCUBA_850/sdfFIles/'
mask_in = '/home/emma/gradu/data/All_data/SCUBA_850/clumps/'

out = '/home/emma/gradu/data/All_data/SCUBA_850/clumps/extractResults/'
#############################
# Functions
#############################



#############################
# Init
#############################


cmd1 = "export STARLINK_DIR=/home/emma/star-2017A;"# echo $STARLINK_DIR;"
#cmd2 = "source $STARLINK_DIR/etc/profile;"
cmd2 = "source $STARLINK_DIR/etc/profile;"

convert = '/home/emma/star-2017A/bin/convert/convert.sh'
cupid = '/home/emma/star-2017A/bin/cupid/cupid.sh'

#init = cmd1+cmd2+cmd4

# Individual routines
fits2ndf = '/home/emma/star-2017A/bin/convert/fits2ndf'
ndf2fits = '/home/emma/star-2017A/bin/convert/fits2ndf'
#findclumps = '/home/emma/star-2017A/bin/cupid/findclumps'

smurf      =  "/home/emma/star-2017A/bin/smurf/smurf.csh"
#convert    =  "/home/emma/star-2017A/bin/convert/convert.csh"
#cupid      =  "/home/emma/star-2017A/bin/cupid/cupid.csh"
kappa      =  "/home/emma/star-2017A/bin/kappa/kappa.csh"
makemap       =  "/home/emma/star-2017A/bin/smurf/makemap"
makesnr       =  "/home/emma/star-2017A/bin/kappa/makesnr"
findclumps    =  "/home/emma/star-2017A/bin/cupid/findclumps"
extractclumps =  "/home/emma/star-2017A/bin/cupid/extractclumps"
thresh        =  "/home/emma/star-2017A/bin/kappa/thresh"
#fits2ndf      =  "/home/emma/star-2017A/bin/convert/fits2ndf"
#ndf2fits      =  "/home/emma/star-2017A/bin/convert/ndf2fits"
findback      =  "/home/emma/star-2017A/bin/cupid/findback"
stats         =  "/home/emma/star-2017A/bin/kappa/stats"
ndfcopy       =  "/home/emma/star-2017A/bin/kappa/ndfcopy"
cmd3 = convert
cmd4 = cupid
init = cmd1+cmd2+cmd3+';'+cmd4+';'




# Data: G001_36+20_96_850.sdf
# Mask: G001_36+20_96_850_clumps.sdf
#############################
# Extractclumps
#############################
for filename in os.listdir(mask_in):

	if filename.endswith('.sdf'):
#	if filename.endswith('G001_36+20_96_850_clumps.sdf'):
		data = filename.strip('_clumps.sdf')+'.sdf'
		print filename, data
		
		outName = filename.strip('_clumps.sdf')+'_out.sdf'
		catName = filename.strip('_clumps.sdf')+'_cat'

		dataIN = data_in+data
		maskIN = mask_in+filename
		fileOUT = out+outName
		cat = out+catName

#shape = ellipse
		
		command = init+extractclumps+' '+maskIN+' '+dataIN+' '+fileOUT+' outcat='+cat+' deconv=false shape=ellipse'
		os.system(command)



















