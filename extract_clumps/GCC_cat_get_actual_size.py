################################################################################
#		Using distance estimates from GCC catalogue
#			calculate actual size of clumps
################################################################################

import os
import math

#################################
# Paths
#################################
file_in = '/home/emma/gradu/data/All_data/SCUBA_850/distances_to_sources/distance_from_PCGG_cat.txt'

dir_in = '/home/emma/gradu/data/All_data/SCUBA_850/clumps/extractResults/results/resultTables/'

dir_out = '/home/emma/gradu/data/All_data/SCUBA_850/distances_to_sources/size_of_clumps/GCC_cat/'


#################################
# Code
#################################

with open(file_in) as f: 
	next(f)
	for line in f: 


		p = line.split()
		

		
		try: 
			d = float(p[1])
		except:
			print p[1]
			continue


		FIELD = p[0]

		FF = FIELD.split('.')

		field = FF[0]+'_'+FF[1]+'_'+FF[2]

		#if field != 'G006_01+36_74':
		#	continue

		scubaF = field+'_850_clumps.txt'


		outF = open(dir_out+field+'_size_of_clumps_GCC.txt','w')
		outF.write('# Distances from GCC catalogue\n')


		outF.write('# Distance = '+str(d)+' [kpc]'+'\n')
		
		outF.write('# Clump ID	size[pc]\n')
		#print d,'+',pl,'= ',dMax
		#print d,'-',m,'= ',dMin

		doc = open(dir_in+scubaF,'r')
		next(doc)
		for l in doc: 
			if l[0] == '#':
				break
			q = l.split()

			clumpID = q[0]
			outF.write(clumpID+'		')

			size = q[3].split('(')[1].strip(')')

			ra,dec = size.split(',')
			
			x,y = float(ra)/3600,float(dec)/3600	# In arcsec: x = ra, y = dec
		
			rx = 1000*d*math.sin(x)	# *1000 because d is in kpc and I want pc
			ry = 1000*d*math.sin(y)

			## max values
			#rxMax = 1000*dMax*math.sin(x)	# *1000 because d is in kpc and I want pc
			#ryMax = 1000*dMax*math.sin(y)

			## min values
			#rxMin = 1000*dMin*math.sin(x)	# *1000 because d is in kpc and I want pc
			#ryMin = 1000*dMin*math.sin(y)

			#rxplus = rxMax - rx
			#rxmi = rx - rxMin

			#ryplus = ryMax - ry
			#rymi = ry - ryMin



			#print clumpID
			
			#print rx,' + ',rxplus,' - ',rxmi

			
			RX = str(round(rx,3))
			RY = str(round(ry,3))
	
			#RXp = str(round(rxplus,3))
			#RXm = str(round(rxmi,3))

			#RYp = str(round(ryplus,3))
			#RYm = str(round(rymi,3))
			
			#if round(rxplus,12) == round(rxmi,12): 
			outF.write('('+RX+','+RY+')'+'\n')
			#else: 
			#	outF.write('('+RX+' + '+RXp+' - '+RXm+','+RY+' + '+RYp+' - '+RYm+')'+'\n')
				

			#if round(ryplus,12) == round(rymi,12): 
			#	outF.write()
			#else: 
		#		outF.write(str(ry)+' + '+str(ryplus)+' - '+str(rymi)+'\n')
				











