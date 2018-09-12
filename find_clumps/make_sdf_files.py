#################################################################################
#		Take all .fits 850 um files and turn them into sdf files
#################################################################################
#				Jul 13, 2018
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
in_850 = '/home/emma/gradu/data/All_data/SCUBA_850/changedNames/'
out_850 = '/home/emma/gradu/data/All_data/SCUBA_850/sdfFIles/'

#############################
# Init
#############################
cmd1 = "export STARLINK_DIR=/home/emma/star-2017A;"# \n'
cmd2 = "source $STARLINK_DIR/etc/profile;"
#cmd3 = 'convert;'
fits2ndf = '/home/emma/star-2017A/bin/convert/fits2ndf'
convert = '/home/emma/star-2017A/bin/convert/convert.sh'
cmd3 = convert
#############################
# Run script
#############################
for filename in os.listdir(in_850):
	if filename.endswith('.fits'):
		
		name_in = in_850+filename
	#	print name_in

		# CODE CANNOT HANDLE . IN FITS FILENAMES! Make sure changenames.py has been run
		# Removing _850.fits fucks up the naming
		name_out = out_850+filename.strip('.fits')+'.sdf'
		#print name_in, '	', name_out
		os.system(cmd1+cmd2+cmd3+';'+fits2ndf+' '+name_in+' '+name_out)







#os.system(cmd1+cmd2+cmd3+' '+fits2ndf+' '+in_850_test+name+' '+in_850_test+name_out)


#os.system('cd /home/emma/gradu/codes/find_clumps')


# config = config_file.txt # fellWalker commands 



