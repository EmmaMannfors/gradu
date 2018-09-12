#################################################################################
#		Take all FellWalker sdf files and makes them .fits files
#################################################################################
#				Jul 20, 2018
#################################################################################
# 	Each instance of os.system() makes a new shell, and so Starlink-
#			commands must be run in each excecution
#################################################################################
#	export STARLINK_DIR=/home/emma/star-2017A
#	source $STARLINK_DIR/etc/profile
#	convert
#	cupid

# /home/emma/gradu/codes/find_clumps





import os
#############################
# paths
#############################
in_850 = '/home/emma/gradu/data/All_data/SCUBA_850/clumps/'
out_850 = '/home/emma/gradu/data/All_data/SCUBA_850/clumps/fits_clumps/'

#############################
# Init
#############################
cmd1 = "export STARLINK_DIR=/home/emma/star-2017A;"# \n'
cmd2 = "source $STARLINK_DIR/etc/profile;"
#cmd3 = 'convert;'
fits2ndf = '/home/emma/star-2017A/bin/convert/fits2ndf'
ndf2fits = '/home/emma/star-2017A/bin/convert/ndf2fits'
convert = '/home/emma/star-2017A/bin/convert/convert.sh'
cmd3 = convert
#############################
# Run script
#############################
for filename in os.listdir(in_850):
	if filename.endswith('.sdf'):
		
		name_in = in_850+filename
	#	print name_in

		# CODE CANNOT HANDLE . IN FITS FILENAMES! Make sure changenames.py has been run
		# Removing _850.fits fucks up the naming
		name_out = out_850+filename.strip('.sdf')+'.fits'
		#print name_in, '	', name_out
		os.system(cmd1+cmd2+cmd3+';'+ndf2fits+' '+name_in+' '+name_out)




