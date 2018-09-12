#######################################################################
#		Remove stars that are in the data twice. 
#		Stars with distace < 2 -3 arcsec should be merged
#######################################################################
#		Do this for all stars as well as clump-associated
#				(no_nan) files
#######################################################################
#				Jul 26, 2018
#######################################################################

import os
import math
from functions import distance_on_sphere_deg
#################################
# Paths
#################################
all_in = '/home/emma/gradu/data/All_data/YSO/clump_and_YSO_list/'
clumps_in = '/home/emma/gradu/data/All_data/YSO/clump_and_YSO_list/no_nan_clumps/'

list_out = '/home/emma/gradu/data/All_data/YSO/clump_and_YSO_list/no_duplicates/'
no_nan_out = '/home/emma/gradu/data/All_data/YSO/clump_and_YSO_list/no_duplicates/no_nan_or_duplicates/'
#################################
# Dictionary
#################################
#stars = dict()
#starList = []
#multiples = dict()
#
r_max = 0.000555555555556

#################################
# Functions
#################################
def multi(ra0,dec0,ra,dec):
	#r0 = math.sqrt(ra0**2 + dec0**2)
	#r1 = math.sqrt(ra**2 + dec**2)
	d = distance_on_sphere_deg(ra0,dec0,ra,dec)
	#print '	'+str(ra),str(dec)
	#print '	'+str(d)
	if d <= r_max: 
		return True
	else: 
		return False


#################################
# Code
#################################

doThis = 'n'



for filename in os.listdir(all_in):
	if filename.endswith('.txt'):
		
	#	if filename != 'G001.36+20.96.txt' and filename != 'G005.93-01.00.txt':
	#		continue
			
		#################################
		if doThis == 'n':
			print 'Remove doThis variable to write files'
			break
		#################################

		#################################
		# Dictionary
		#################################
		stars = dict()
		starList = []
		multiples = dict()
		matchList = []
		#################################

		i = 0 # For naming Akari
		
		f = open(all_in+filename,'r')
		outName = filename.strip('.txt')+'_no_duplicates.txt'
		#print outName
		fOut = open(list_out+outName,'w')

		typ = ''
		
		for line in f: 
			p = line.split()
			#print p[1]
			if p[1] == 'Vizier' and p[2] == 'I/II':#line == '# Vizier I/II ':
				print 'meow'
				typ = 'vizI_II'
			if p[1] == 'Vizier' and p[2] == 'III':
				typ = 'viz_III'
			if p[1] == 'Akari' and p[2] == 'FIS':
				typ = 'FIS'
			if p[1] == 'Akari' and p[2] == 'PSC':
				typ = 'PSC'
			if p[1] == 'Simbad':
				typ = 'sim'
				print 'meow'
			else: 
				typ = typ

			if line[0] == '#':
		#		fOut.write(line)
				continue

			
			if p[0] == '-':
				name = 'Akari_'+str(i)
			else: 
				name = p[0]
			i += 1

			ra = float(p[1])
			dec = float(p[2])
			clump = float(p[3])
			if typ == 'sim':
				typ += '_'+p[4]
			#try:
			#	typ = p[4]
			#except: 
			#	typ = ''
			#print name

			stars[name] = (ra,dec,clump,typ)
			starList.append(name)

			valList = stars.values()

			# ra = valList[0][0], dec = valList[0][1]
			#print valList[0][0],valList[0][1]

		# After running through file, check dict
		for key in stars: 
			if key in matchList: 
				continue
		#	print key
			valList = stars[key]
			ra0 = valList[0]
			dec0 = valList[1]

		#	print ra0,dec0

			##################################
			# Makes dict entry for every key
			# with value as a blank list
			##################################
			multiples[key] = []
			
			for star in stars: 
				if star != key: 
					ra,dec = stars[star][0],stars[star][1]
					#print '		'+star
					if multi(ra0,dec0,ra,dec):
						print '	match!'
						# List of already matched objects
						matchList.append(star)
						multiples[key].append(star)
		#	try: 
		#		print multiples[key]
		#	except: 
		#		print 'no match'

		#test = open('test.txt','w')
		fOut.write('# Name		ra		dec		clump	type\n')
		
		for key in multiples: 
			raList, decList,clumpList = [],[],[]
			fOut.write('##########\n')
			s = stars[key]
			# Make a list of associated ra/dec values
			raList.append(float(s[0]))
			decList.append(float(s[1]))
			clumpList.append(float(s[2]))
			vals = multiples[key]
			if len(vals) != 0:
				fOut.write('# '+key+'	')
			else: 
				fOut.write(key+'	')
			for n in range(0,len(s)): 
				fOut.write(str(s[n])+'	')
			fOut.write('\n')
			for star in vals: 
				fOut.write('# '+star+'	')
				q = stars[star]
				raList.append(float(q[0]))
				decList.append(float(q[1]))
				clumpList.append(float(q[2]))
				for n in range(0,len(q)):
					fOut.write(str(q[n])+'	')
				fOut.write('\n')

			############################
			# avg ra and dec for matches
			############################
			ra_tot,dec_tot,clump_tot = 0.0,0.0,0.0
			for n in range(0,len(raList)):
				ra_tot += raList[n]
			#print raList
			ra_avg = ra_tot/len(raList)

			for n in range(0,len(decList)):
				dec_tot += decList[n]
			#print decList
			dec_avg = dec_tot/len(decList)


			for n in range(0,len(clumpList)):
				clump_tot += clumpList[n]
			print clumpList
			clump_avg = clump_tot/len(clumpList)		
			#print clump_avg
			
			if len(vals) != 0:
				fOut.write('match		'+str(ra_avg)+'	'+str(dec_avg)+'	'+str(clump_avg)+'	- \n')
			#print 5*'*'





#################################
# Make files without nan clumps
#################################
#makeFiles = raw_input('Make list of only existing clumps? (y/n)')

makeFiles = 'y'
if makeFiles == 'y':
	for filename in os.listdir(list_out):
	#	if filename == 'G001.36+20.96_no_duplicates.txt':
		if filename.endswith('.txt'):

		#	if filename != 'G001.36+20.96_no_duplicates.txt':
		#		continue
		#	f = open(list_out+filename,'r')
			outName = filename.strip('.txt')+'_no_nan.txt'
			outF = open(no_nan_out+outName,'w')

		#	print list_out+filename
		#	print no_nan_out+outName

			print filename
			with open(list_out+filename) as f:
				for line in f: 
				#	print line
		#	for line in file_in: 
		#		print line
		#		#print 'meow'
		#		print line[0]
					p = line.split()
					if line == '##########':
						continue
					elif line[0] == '#' or line[0] == '*':

						###################
						# Remove nan-matches too
						###################
						try: 
							if p[4] == 'nan':
								continue
						except: 
							print ''
						###################
						outF.write(line)
						continue
					clump = p[3]
					if clump == 'nan':
						continue
					outF.write(line)
					#print clump
				
			
			
		
#no_nan_out








































				
		
		
