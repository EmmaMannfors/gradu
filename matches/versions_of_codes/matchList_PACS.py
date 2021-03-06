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
		# SCUBA = p[2] (0 = line #, 1 = NRO name)
		if str(p[2]) not in d: #If dict does not contain key yet
			space = '	'*5
			d[str(p[2])] = ['NA'+space,'NA'+space,'NA'+space]

		
print len(d)







# Print 4 items from ditionary d
#n_items = list(islice(d.iteritems(), 4))
#print n_items



#for line in f: 
#		p = line.split()
#		SCUBAList.append(str(p[0]))
#		HerschelList.append(str(p[1]))




#for filename in os.listdir('/home/emma/gradu/data/All data'):
#	# if filename ends with .txt and contains the string '70'
#	if filename.endswith('.txt') and '70' in filename:
#		with open(filename) as f:
#			next(f)
#			for line in f:
#				p = line.split()
#				d[str(p[0])] = [str(p[1]),'NA']
#		#print filename
#		continue
#	elif filename.endswith('.txt') and '100' in filename:
#		with open(filename) as f:
#			next(f)
#			for line in f:
#				p = line.split()
#				valList = d[p[0]]
#				val70 = valList[0]
#				d[str(p[0])] = [val70,str(p[1])]
#		continue
#	else:
#		continue


def addToDict(filename,wavelength,i):
#70, i=0; 	100, i=1; 	160, i=2; 	250, i=3; 	350, i=4; 	500, i=5
	lambdaList = [70,100,160,250,350,500]
	index = 0
	with open(filename) as f:
		next(f)
		for line in f:
			p = line.split()
			valList = d[p[0]]
			values = []
			while index <= 5:
				if index == i:
					values.append(str(p[index]))
					continue
				else:
					values.append(valList[index])
			#print values
			d[str(p[0])] = [values[0],values[1],values[2]]





with open('Final_sources_70.txt') as f:
	next(f)
	for line in f:
		p = line.split()
		valList = d[p[0]]
		val100 = valList[1]
		val160 = valList[2]
		d[str(p[0])] = [str(p[1])+'	',val100,val160]

with open('Final_sources_100.txt') as f:
	next(f)
	for line in f:
		p = line.split()
		valList = d[p[0]]
		val70 = valList[0]
		val160 = valList[2]
		d[str(p[0])] = [val70,str(p[1])+'	',val160]

with open('Final_sources_160.txt') as f:
	next(f)
	for line in f:
		p = line.split()
		valList = d[p[0]]
		val70 = valList[0]
		val100 = valList[1]
		d[str(p[0])] = [val70,val100,str(p[1])+'	']





#############################################################################
#		Write dictionary to file
#############################################################################

final = open('PACSMatches.txt','w')
space3 = '	'*3
space4 = '	'*4
final.write('SCUBA-2'+space3+'  PACs 70'+space4+'    PACS 100'+space4+'    PACS 160\n')
for key in d:
	valList = d[key]
	val70 = valList[0]
	val100 = valList[1]
	val160 = valList[2]
	space = ' '*4
	final.write(key + space + val70 + space + val100 + space + val160 + '\n')

final.close()














 
