
f = open('centers_for_vizier.txt','w')

with open('center_in_deg.txt') as c: 
	next(c)
	for line in c:
		p = line.split()
		f.write(': '+p[0]+'\n')
		f.write(p[1]+'	'+p[2]+'\n')
		
