#################################################################################
#		Change filenames so they don't have . in them 
#			except for the ending
#################################################################################
#			jul 17, 2018
#################################################################################
import os



##################################
# Paths
##################################
#in_850_test = '/home/emma/gradu/data/test/starlinkTest/'
in_850 = '/home/emma/gradu/data/All_data/SCUBA_850/'
out_850 = '/home/emma/gradu/data/All_data/SCUBA_850/changedNames/'


for filename in os.listdir(in_850):
	if filename.endswith('.fits'):
		f = filename.split('.')
		newName = f[0]+'_'+f[1]+'_'+f[2]+'.fits'
		if not newName.endswith('_850.fits'):
			print 'ERROR!!!! \n'+filename+'\n'+newName
		#print filename
		#print newName
		print 10*'*'
		os.system('cd '+in_850+';'+'cp '+filename+' '+out_850+newName)
		#print 'cd '+in_850+';'+'cp '+filename+' '+out_850+newName
