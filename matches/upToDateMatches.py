#############################################################################
#		Script for making an up-to-date version of
#				matches.txt
#############################################################################
#############################################################################
#		Take the SCUBA, 70 - 160 um matches from
#		matches.txt and the correct matches from
#		the SPIRE images from SPIREmatches.txt
#############################################################################
#				Jun 18, 2018
#############################################################################
import os

#################################
#	Paths
#################################
match_in = '/home/emma/gradu/codes/matches/cropMatches.txt'
spire_in = '/home/emma/gradu/codes/matches/SPIREMatches.txt'


#new = open('newMatches.txt','w')
space3 = 3*'	'
space4 = 4*'	'
spaces = 5*'	'
#new.write('SCUBA-2'+space3+'PACs 70'+spaces+'PACS 100'+space4+'PACS 160'+space4+'SPIRE PSW (250)'+spaces+'SPIRE PMW (350)'+spaces+'SPIRE PLW (500)'+'\n')

spire = open(spire_in,'r')

none = 'NA				'
spc = '	'
#with open(match_in) as f:
#	next(f)
#	next(f)
#	next(f)
#	for line in f:
#		p = line.split()
#		if p[1] == 'NA':
#			p[1] = 'NA				'
#		if p[2] == 'NA':
#			p[2] = 'NA				'
#		if p[3] == 'NA':
#			p[3] = 'NA				'
#
#		for line in spire:
#			q = line.split()
#			if p[0] == q[0]:
#				new.write(p[0]+spc+p[1]+spc+p[2]+spc+p[3]+spc+q[1]+spc+q[2]+spc+q[3]+'\n')
#				continue
#			else:
#				continue
#		
#		new.write(p[0]+spc+p[1]+spc+p[2]+spc+p[3])
#		continue
		






# Create empty dictionary
d = dict()



# Add every SCUBA file as the key to a separate dictionary entry.
# Every entry takes 6 string which will be the wavelength filenames. 
# Currently every single string is 'NA'
with open(match_in) as matches:
	next(matches)
	next(matches)
	next(matches)
	for line in matches:
		p = line.split()
		if p[1] == 'NA':
			p[1] = 'NA				'
		if p[2] == 'NA':
			p[2] = 'NA				'
		if p[3] == 'NA':
			p[3] = 'NA				'
		if str(p[0]) not in d: #If dict does not contain key yet
			space = '	'*5
			d[str(p[2])] = [p[1],p[2],p[3],none,none,none]

		

#############################################################################
#		Function which adds everything to dictionary. 
#############################################################################
def addToDict(filename,wavelength,i):
	with open(filename) as f:
#		next(f)
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
f250 = open('sources_250.txt','w')
f350 = open('sources_350.txt','w')
f500 = open('sources_500.txt','w')
#f250.write('SCUBA			SPIRE PSW')
#f350.write('SCUBA			SPIRE PMW')
#f500.write('SCUBA			SPIRE PLW')


with open(spire_in) as f:
	next(f)
	for line in f:
		t = line.split()
		f250.write(t[0]+'	'+t[1]+'\n')
		f350.write(t[0]+'	'+t[2]+'\n')
		f500.write(t[0]+'	'+t[3]+'\n')


addToDict('sources_250.txt',250,3)
addToDict('sources_350.txt',350,4)
addToDict('sources_500.txt',500,5)




#############################################################################
#		Write dictionary to file
#############################################################################

final = open('newMatches.txt','w')
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

		


