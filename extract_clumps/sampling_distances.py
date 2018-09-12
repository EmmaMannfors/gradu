############################################################################
#			Are there matches in sampling survey
############################################################################
#				Aug 2, 2018
############################################################################
from functions import distance_on_sphere_deg

IN = '/home/emma/Desktop/gal_coord.txt'


print 'field,l,b,DR1_(l,b),distance'
with open(IN) as g: 
	for line in g: 
		cat = open('./DR1_source_info.txt','r')
		p = line.split()
		if line[0] == ':':
			field = p[1]
			continue
		else: 
			field = field
			l = float(p[0])
			b = float(p[1])


		#print field, l, b
		#if field != 'G202.31+02.53':
		#	continue
		#if field != 'G105.44+09.88':
		#	continue
		#print l,type(l)
		#print field,l,b

		for li in cat: 
			if li[0] == '#': 
				continue
			e = li.split()
			#try: 
			lCat = float(e[0])
			bCat = float(e[1])
			#except: 
			#	#print e[0]
			#	print 'ahh'

			#if field == 'G202.31+02.53':
			#print lCat, bCat

			#print lCat
			#L = float(l)

			dist = distance_on_sphere_deg(l, b, lCat, bCat)
			
			if dist <= 0.17:#10 arcmin
				print field#,l,b,lCat,bCat,dist
				print li
				continue
			

#DD = distance_on_sphere_deg(105.44, 9.88, 105.415, 0.87943)
#print DD
