##########################################################################
#		Find and store info about YSOs and others 
#				from Simbad	
##########################################################################
#			partly from python course final
##########################################################################
#			Jul 9, 2018
##########################################################################


#########################################################
#			YSO dictionary: 
#########################################################
#		Key to dictionary is source name
#	Value of dictionary  contains info about the YSO: 
#		- name (key)
#		- RA (J2000)
#		- dec (J2000)
#		- Units of coordinate
#		- Type of source
#		- Flux
#		- Link to Simbad page
#########################################################
#			Jul 7, 2018
#########################################################
#			Emma Mannfors
#	 		013744560
#		Python ohjelmointia aloittelijoille
#########################################################

import os.path

###########################
# Set up YSO class
###########################
class newStar:
	def __init__(self,name,co,unit,typ,f,link):
		self.name = name
		self.coordinate = co
		self.unit = unit
		self.typ = typ
		self.flux = f
		self.link = link





###########################
# Set up YSO dictionary
# and list of star-class
# objects
###########################
starDic = dict()
starsList = []

		

###########################
# Create new YSO object
###########################
def newYSO():
	name = input('Name of source\n')
	co = input('Coordinates\n')
	unit = input('units of coordinates (deg/sexagesimal(/s))\n')
	typ = input('Type of source? (YSO/YSO candidate/other)\n')
	f = input('flux\n')
	link = input('Link to Simbad page\n')
	star = newStar(name,co,unit,typ,f,link)
	starDic[name] = star
	starsList.append(star)

###########################
# Functions to write
# to document
###########################

# The actual text that is written to document
def w(document):
	doc = open(document,'w')
	for star in starsList:
		doc.write(star.name+': '+star.typ+'\n')
		doc.write('(ra,dec) ['+star.unit+']: '+star.coordinate+'\n')
		doc.write('flux: '+star.flux+'\n')
		doc.write(star.link+'\n')

# The function to open document and do the writing
def writeToDocument():
	docName = input('Input document\n (including path!)\n')
	if os.path.isfile(docName):
		yn = input('WARNING! Running write command will overwrite previous document! Proceed? (y/n)\n')
		if yn == 'y':
			w(docName)
			print('Directory written to '+docName)
		else:
			print('File not written')	
	else:	
		w(docName)
		print('Directory written to '+docName)

###########################
# Prints all current
# YSO objects to terminal
###########################
def printYSO():
	for star in starsList:
		print(5*'*')
		print(star.name,': ',star.typ)
		print('(ra,dec) [',star.unit,']: ',star.coordinate)
		print('flux: ',star.flux)
		print(star.link)
	print(5*'*')
		
###########################
# Converting ra,dec to 
# degrees. Used for 
# list(item=coordinate)
###########################
def ra_to_deg(ra):
	g = ra.split(':')
	h = float(g[0])
	m = float(g[1])
	s = float(g[2])
	h += (m/60.0)
	h += (s/3600.0)
	deg = h*(360.0/24.0)
	return deg

def dec_to_deg(dec):
	y = dec.split(':')
	deg,amin,asec = float(y[0]),float(y[1]),float(y[2])
	if deg < 0:
		deg -= (amin/60.0)
		deg -= (asec/3600.0)
		return round(deg,6)
	deg += (amin/60.0)
	deg += (asec/3600.0)
	return deg
	

###########################
# Lists features of 
# YSO objects
###########################
def list(item):
	# Listing coordinates: all are changed to degrees first!
	if item == 'coordinate' or item == 'co':
		raList = []
		decList = []
		ra_source_dic = dict()
		dec_source_dic = dict()
		print('	')
		for star in starsList:
			if star.unit == 'sexagesimal' or star.unit == 's':
				co = star.coordinate
				a = co.strip('(').strip(')')
				b = a.split(',')
				ra_s,dec_s = b[0],b[1]
				ra,dec = ra_to_deg(ra_s),dec_to_deg(dec_s)
				
			elif star.unit == 'deg':	
				co = star.coordinate
				a = co.strip('(').strip(')')
				b = a.split(',')
				ra,dec = float(b[0]),float(b[1])
			else: 
				print(star.name+': units must be in deg or s\n')
				continue

			if not ra in raList:
				raList.append(ra)
				ra_source_dic[ra] = star.name
			if not dec in decList:
				decList.append(dec)
				dec_source_dic[dec] = star.name
					
		ra_sort = sorted(raList)
		print(5*'*'+'\n'+'RA in ascending order\n'+5*'*')
		for i in range(0,len(ra_sort)):	
			r = ra_sort[i]
			print(ra_source_dic[r],'	',r)
			
		print(' ')

		d_sort = sorted(decList)
		print(5*'*'+'\n'+'DEC in ascending order\n'+5*'*')
		for i in range(0,len(d_sort)):	
			d = d_sort[i]
			print(dec_source_dic[d],'	',d)

	# Listing flux. Assumes all are in same units.
	elif item == 'flux':
		fluxList = []
		f_s_dic = dict()
		print('	')
		for star in starsList:
			if star.flux != 'NA':
				fluxList.append(float(star.flux))
				f_s_dic[float(star.flux)] = star.name
			else: 
				print(star.name+': flux unknown\n')

		f_sort = sorted(fluxList)
		print(5*'*'+'\n'+'Flux in ascending order\n'+5*'*')
		for i in range(0,len(f_sort)):
			f = f_sort[i]
			print(f_s_dic[f],'	',f)

	# Listing sources in alphabetical order
	elif item == 'source':
		print('	')
		print(5*'*'+'\n'+'Sources\n'+5*'*')
		sourceList = []
		for star in starsList:
			sourceList.append(star.name)
		s_sort = sorted(sourceList)
		for line in s_sort:
			print(line)
	
	# Listing sources by type
	elif item == 'type':
		ysoList = []
		yso_cList = []
		otherList = []
		typeList = []
		print('	')
		print(5*'*'+'\n'+'Objects by type:\n'+5*'*')
		for star in starsList:
			if not star.typ in typeList:
				typeList.append(star.typ)
		for item in typeList:
			print(item+':')
			for star in starsList:
				if star.typ == item:
					print(star.name)
			print('	')
	
	# Listing links, though in no particular order
	elif item == 'link':
		print('	')
		print(5*'*'+'\n'+'SIMBAD links to YSOs\n'+5*'*')
		for star in starsList:
			print('	')
			print(star.name)
			print(star.link)


	
###########################
# Removes YSO object
###########################
def removeYSO(key):
	try:
		n = starDic[key]
		del starDic[key]
		starsList.remove(n)

	# In case user tries to input nonexisting source
	except:
		print('\n **Source does not exist**')
		list('source')	

###########################
# Changes values of 
# YSO objects
###########################
def overWrite(name):
	# Changes star's info
	try:
		star = starDic[name]
	
		print('leave answer blank if unchanged\n')
		newName = input('New name\n')
		if newName == '':
			newName = name
		newCo = input('New coordinates\n')
		if newCo == '':
			newCo = star.coordinate
			newUnits = star.unit
		else: 
			newUnits = input('Units of new coordinates\n')
			if newUnits == '':
				newUnits = star.unit
		newType = input('New type?\n')
		if newType == '':
			newType = star.typ
		newFlux = input('New flux?\n')
		if newFlux == '':
			newFlux = star.flux
		newLink = input('New link?\n')
		if newLink == '':
			newLink = star.link
		# Remove references to old object
		removeYSO(name)
		
		# Make new object and add it to lists
		star_new = newStar(newName,newCo,newUnits,newType,newFlux,newLink)
		starDic[newName] = star_new
		starsList.append(star_new)

	# In case user makes a typo or tries to change source which does not exist. A list of existing sources is printed out. 
	except:
		print('\n **Source does not exist**')
		list('source')




###########################
# For testing purposes
###########################
def testYSO_1():
	name = '2Mass 0123.4'
	co = '(01:34:56.6,+12:25:46.00)'
	unit = 's'
	typ = 'YSO'
	f = '666'
	link = 'www.cheese.com'
	star = newStar(name,co,unit,typ,f,link)
	starDic[name] = star
	starsList.append(star)

def testYSO_3():
	name = 'MJR2015 1446'
	co = '(324.69876,-34.212)'
	unit = 'deg'
	typ = 'YSO candidate'
	f = '27.5'
	link = 'www.cats.com'
	star = newStar(name,co,unit,typ,f,link)
	starDic[name] = star
	starsList.append(star)

def testYSO_2():
	name = 'Cthulhu'
	co = '(01:34:56.6,+12:25:46.00)'#'(b,a)'
	unit = 's'
	typ = 'YSO candidate'
	f = '100'
	link = 'www.cats.com'
	star = newStar(name,co,unit,typ,f,link)
	starDic[name] = star
	starsList.append(star)
###########################

##################
testYSO_1()
testYSO_2()
testYSO_3()
##################


while True:
	print('\n')
	do = input('What to do?\n(new/remove/print/stop/write/change/list)\n')
	if do == 'new':
		newYSO()
	elif do == 'remove':
		key = input('What YSO to remove? (give name of YSO)\n')
		removeYSO(key)
	elif do == 'stop':
		break
	elif do == 'print':
		printYSO()
	elif do == 'write':
		writeToDocument()
	elif do == 'change':
		name = input('Name of source\n')
		overWrite(name)
	elif do == 'list':
		thing = input('What to list?\n(source/coordinate/type/flux/link)\n')
		list(thing)

















