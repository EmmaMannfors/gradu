##############################################################################
#		Find number of sourcse for field
##############################################################################
#			Jul 9, 2018
##############################################################################

import os

path_in = '/home/emma/gradu/data/All_data/catalogs/Akari/'
file_out = '/home/emma/gradu/data/All_data/catalogs/Akari/number_of_sources'

num = open(file_out,'w')
num.write('field		FIS	PSC\n')
for filename in os.listdir(path_in):
	fis = 0
	psc = 0
	if filename.endswith('.txt'):
		f = filename.split('_')
		field = f[0]
		
		# FIS files
		with open(path_in+field+'_FIS.txt') as now:
			for line in now: 
				fis+=1
		
		# PSC catalog
		with open(path_in+field+'_PSC.txt') as current:
			for line in current: 
				psc+=1
			
		# There are 7 lines of info	
		fis -= 7
		psc -= 7
		
		print field+'	'+str(fis)+'	'+str(psc)+'\n'
		num.write(field+'	'+str(fis)+'	'+str(psc)+'\n')


