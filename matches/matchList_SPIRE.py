#############################################################################
#        Script for going thru Final_File_wavelength and 
# 	combining them into one matches.txt-style file
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#				May 28, 2018
#############################################################################

import os
from itertools import islice    #For printing random 4 items

#############################################################################
#			Go thru every single list
#			Find all SCUBA-2 files and make hashmap (?)
#			of all the files with corresponding 
#			Herschel files. 
#############################################################################
#			If H. file does not exist, value should
#			become NA
#############################################################################

# Create empty dictionary
d = dict()

#How to append shit to dictionary
#d['a'] = ['b','c']



with open('/home/emma/gradu/data/matches.txt') as matches:
	next(matches)
	next(matches)
	next(matches)
	for line in matches:
		p = line.split()
		if str(p[2]) not in d: #If dict does not contain key yet
			space = '	'*5
			d[str(p[2])] = ['NA'+space,'NA'+space,'NA'+space]

		
print len(d)






with open('Final_sources_250.txt') as f:
	next(f)
	for line in f:
		p = line.split()
		valList = d[p[0]]
		val350 = valList[1]
		val500 = valList[2]
		d[str(p[0])] = [str(p[1])+'	',val350,val500]

with open('Final_sources_350.txt') as f:
	next(f)
	for line in f:
		p = line.split()
		valList = d[p[0]]
		val250 = valList[0]
		val500 = valList[2]
		d[str(p[0])] = [val250,str(p[1])+'	',val500]

with open('Final_sources_500.txt') as f:
	next(f)
	for line in f:
		p = line.split()
		valList = d[p[0]]
		val250 = valList[0]
		val350 = valList[1]
		d[str(p[0])] = [val250,val350,str(p[1])+'	']





#############################################################################
#		Write dictionary to file
#############################################################################

final = open('SPIREMatches.txt','w')
space3 = '	'*3
space4 = '	'*4
final.write('SCUBA-2'+space3+'  SPIRE PSW (250)'+space4+'   SPIRE PMW (350)'+space4+'    SPIRE PLW (500)\n')
for key in d:
	valList = d[key]
	val250 = valList[0]
	val350 = valList[1]
	val500 = valList[2]
	space = ' '*4
	final.write(key + space + val250 + space + val350 + space + val500 + '\n')

final.close()














 
