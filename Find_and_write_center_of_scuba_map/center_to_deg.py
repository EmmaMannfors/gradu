f = open('center_in_deg.txt','w')
f.write('field			RA		DEC\n')

with open('converted_centers.txt') as c:
	next(c)
	next(c)
	for line in c:
		p = line.split()
		name = p[0]
		ra_i = p[1]+':'+p[3]+':'+p[5]	
		dec_i = p[7]
		
		# RA to deg
		# 	First, convert all to hours
		
		h = float(p[1])+(float(p[3])/60.0)+(float(p[5])/3600.0)
		ra = h*(360.0/24.0)

		d = dec_i.split(':')
		if float(d[0]) < 0:
			dec = float(d[0])-(float(d[1])/60.0)-(float(d[2])/3600.0)
		else: 
			dec = float(d[0])+(float(d[1])/60.0)+(float(d[2])/3600.0)
		
		f.write(name+'	'+str(ra)+'	'+str(dec)+'\n')
