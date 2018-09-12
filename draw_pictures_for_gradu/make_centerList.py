#####################################################################################
#		Writes file which has scuba centers corresponding with Herschel files
#####################################################################################
#				Jul 20, 2018
#####################################################################################

import os



# Centers, in deg, corresponds to SCUBA files
center_deg = '/home/emma/gradu/codes/Find_and_write_center_of_scuba_map/center_in_deg.txt'

# Matches, scuba = p[0], psw = p[4]


out = open('centers_250_850.txt','w')

out.write('SCUBA'+3*'	'+'PSW(250 um)'+4*'	'+'RA		DEC\n')

with open(center_deg) as c: 
	next(c)
	for line in c: 
		match = open('/home/emma/gradu/codes/matches/cropMatches.txt','r')
		p = line.split()
		scuba = p[0]
		ra = p[1]
		dec = p[2]	
		
		print scuba
		
		for line in match: 
			e = line.split()
			print e[0]
			if e[0] == scuba: 
				psw = e[4]
				if psw == 'NA':		
					psw = 'NA'+4*'	'
				
			else: 
				continue
		print psw
		out.write(scuba+'	'+psw+'	'+ra+'	'+dec+'\n')
		print 10*'*'
		psw = 'NA'
