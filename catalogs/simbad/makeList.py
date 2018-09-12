################################################################
#		Make list of sources that will make 
#			searching Simbad easier
################################################################
#			Jul 9, 2018
################################################################

in_File = '/home/emma/gradu/codes/Find_and_write_center_of_scuba_map/converted_centers.txt'

out_File = 'scuba_centers_simbad.txt'

f = open(out_File,'w')
f.write('SCUBA			ra		dec\n')

sDic = dict()
scubaList = []

with open(in_File) as matches: 
	next(matches)
	next(matches)
	for line in matches:
		p = line.split()
		scuba = p[0]
		h,m,s = p[1],p[3],p[5]
		dec = p[7]
		d = dec.split(':')
		if int(d[0]) >= 0:
			#print d[0]
			d[0] = '+'+d[0]
			
		#else: 
			
			#print d[0]
		ra = h+' '+m+' '+s
		dec = d[0]+' '+d[1]+' '+d[2]
		
		scubaList.append(scuba)
		sDic[scuba] = ra+'	'+dec
		#print scuba+' '+h+' '+m+' '+s+' '+d[0]+' '+d[1]+' '+d[2]+'\n'
		#f.write(scuba+'	'+h+' '+m+' '+s+'	'+d[0]+' '+d[1]+' '+d[2]+'\n')

scuba_sort = sorted(scubaList)
for i in range(0,len(scuba_sort)):
	key = scuba_sort[i]
	co = sDic[key]
	print key+'	'+co
	f.write(key+'	'+co+'\n')



















