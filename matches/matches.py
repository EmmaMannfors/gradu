#############################################################################
#		Script for going thru Final_File_wavelength and 
#		combining them into one matches.txt-style file
#
#				For all six wavelengths
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#				May 29, 2018
#############################################################################

import os


# Create empty dictionary
d = dict()



# Add every SCUBA file as the key to a separate dictionary entry.
# Every entry takes 6 string which will be the wavelength filenames. 
# Currently every single string is 'NA'
with open('/home/emma/gradu/data/matches.txt') as matches:
	next(matches)
	next(matches)
	next(matches)
	for line in matches:
		p = line.split()
		# SCUBA = p[2] (0 = line #, 1 = NRO name)
		if str(p[2]) not in d: #If dict does not contain key yet
			space = '	'*5
			d[str(p[2])] = ['NA'+space,'NA'+space,'NA'+space,'NA'+space,'NA'+space,'NA'+space]

		

#############################################################################
#		Function which adds everything to dictionary. 
#############################################################################
def addToDict(filename,wavelength,i):
	with open(filename) as f:
		next(f)
		for line in f:
			p = line.split()
			valList = d[p[0]]
			index = 0
			values = []
			while index <= 5:
				if index == i:
					values.append(str(p[1]))
				else:
					values.append(str(valList[index]))
				index += 1
			d[str(p[0])] = [values[0],values[1],values[2],values[3],values[4],values[5]]


#############################################################################
#		Running function for every wavelength file
#############################################################################
addToDict('Final_sources_70.txt',70,0)
addToDict('Final_sources_100.txt',100,1)
addToDict('Final_sources_160.txt',160,2)
addToDict('Final_sources_250.txt',250,3)
addToDict('Final_sources_350.txt',350,4)
addToDict('Final_sources_500.txt',500,5)




#############################################################################
#		Write dictionary to file
#############################################################################

final = open('cropMatches.txt','w')
space3 = '	'*3
space4 = '	'*4
final.write('SCUBA-2'+space3+'  PACs 70'+space4+'    PACS 100'+space4+'    PACS 160')
final.write(space4+'SPIRE PSW (250)'+space4+'SPIRE PMW (350)'+space4+'SPIRE PLW (500)\n')
for key in d:
	valList = d[key]
	val70 = valList[0]
	val100 = valList[1]
	val160 = valList[2]
	val250 = valList[3]
	val350 = valList[4]
	val500 = valList[5]
	space = ' '*4
	final.write(key + space + val70 + space + val100 + space + val160 )
	final.write(space + val250 + space + val350 + space + val500+'\n')

final.close()



#############################################################################
#		Makes LaTex table compatible .txt file
#############################################################################
#tex = open('matchesforLatex.txt','w')

#space = ' & '
#tex.write('SCUBA-2'+space+'PACs 70'+space+'PACS 100'+space+'PACS 160')
#tex.write(space+'SPIRE PSW (250)'+space+'SPIRE PMW (350)'+space+'SPIRE PLW (500)'+'\n')
#tex.write('\hline')
#for key in d:
#	valList = d[key]
#	val70 = valList[0]
#	val100 = valList[1]
#	val160 = valList[2]
#	val250 = valList[3]
#	val350 = valList[4]
#	val500 = valList[5]
#	tex.write(key + space + val70 + space + val100 + space + val160 )
#	tex.write(space + val250 + space + val350 + space + val500+'\n')
#	tex.write('\hline'+'\n')

#tex.close()







 
