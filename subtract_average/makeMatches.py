###################################################################
#	makes list of matching 500 and 160 um files
###################################################################
matchPath = '/home/emma/gradu/codes/matches/'
coordinateFile = 'emptySpaces.txt'

f = open('500_160_matches.txt','w')

space = 4*'	'
f.write('PACS 160um'+space+'SPIRE PSW'+space+'SPIRE PMW'+space+'SPIRE PLW\n')

########################################
#	Add matches to dictionary
########################################
d = dict()

with open(matchPath+'cropMatches.txt') as matches:
	next(matches)
	for line in matches:
		t = line.split()
	#	print t[0]
		d[str(t[0])] = [t[1],t[2],t[3],t[4],t[5],t[6]]




########################################
#	Run code
########################################
locList = []

with open(coordinateFile) as co:
	next(co)
	for line in co:
		# Getting all the matches
		p = line.split()
		scuba = p[0]+'_850.fits'
		matchList = d[scuba]
		if matchList[2] != 'NA':
			f.write(matchList[2]+'	'+matchList[3]+'	'+matchList[4]+'	'+matchList[5]+'\n')



