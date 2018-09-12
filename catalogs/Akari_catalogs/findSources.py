##############################################################################
#		Find the matches of SCUBA maps in AKARI
#			catalogs. 
##############################################################################
#			Jul 6, 2018
##############################################################################

import os
########################################
# Paths
########################################
matchDir = '/home/emma/gradu/codes/Find_and_write_center_of_scuba_map/'
SCUBA_in = '/home/emma/gradu/data/All_data/SCUBA_850/'


Akari_in = '/home/emma/gradu/Akari_catalogs/'
Akari_FIS = Akari_in+'Akari_FIS_bright_source.txt'
Akari_PSC = Akari_in+'AKARI-IRC_PSC_b2.1_rc6_public.txt'

########################################
# Functions: unit conversions
########################################
def min_to_arcmin(m):
	amin = m*(360.0/24.0)
	return amin

def asec_to_sec(asec):
	s = asec*(24.0/360.0)
	return s

########################################
# radius of map
########################################
# radius of map, 15/2 arcmin
#r_dec = (15.0/2.0)		# in arcmin
#r_ra = 	amin_to_min(r_dec)	# in minutes

r_dec_amin = 7				# 7.5 arcmin
r_dec_asec = 30
r_dec = (15.0/2.0)*(60.0)		# in arcsec
r_ra = asec_to_sec(r_dec)		# in seconds
########################################
# Functions: max and min values
########################################
def ra_max(h,m,s):
	r = r_ra	# in sec now
	if s+r < 60.0:
		seconds = s+r
		minutes = m
		hours = h
	else:
		seconds = r-(60.0-s)
		if m == 59.0:
			minutes = 0.0
			hours = h+1
			
		else:	
		# If there were sources with h = 23, you'd need to deal with it
			minutes = m+1
			hours = h
	return hours,minutes,seconds
	

def ra_min(h,m,s):
	r = r_ra	# sec
	if s-r > 0.0:
		seconds = s-r
		minutes = m
		hours = h
	else:
		seconds = 60.0-(r-s)
		if m == 0.0:
			minutes = 59.0
			hours = h-1
		else:
			minutes = m-1
			hours = h
		
		#minutes = 60.0-(r-m)
	return hours,minutes,seconds
		

def dec_max_positive(deg,amin,asec):
### ONLY FOR POSITIVE DECLINATION!
	r_amin = 7.0
	r_asec = 30.0 # = 7.5 arcmin
	if amin + r_amin < 60.0:
		arcmin = amin + r_amin
		degrees = deg
		if asec+r_asec < 60.0:
			arcsec = asec+r_asec
		else: 
			arcsec = r_asec - (60-asec)
			arcmin +=1
	else: 
		arcmin = (60.0-amin) 
		degrees = deg + 1
		if asec+r_asec < 60.0:
			arcsec = asec+r_asec
		else: 
			arcsec = r_asec - (60-asec)
			arcmin +=1

	return degrees,arcmin,arcsec
		

def dec_min_positive(deg,amin,asec):
	r_amin = 7.0
	r_asec = 30.0 # = 7.5 arcmin
	if amin-r_amin > 0.0:
		arcmin = amin-r_amin
		degrees = deg
		if asec - r_asec > 0.0:
			arcsec = asec - r_asec
		else:
			arcsec = 60.0-(r_asec-asec)
			arcmin -= 1
	else:
		arcmin = 60.0-(r_amin-amin)
		degrees = deg - 1
		if asec - r_asec > 0.0:
			arcsec = asec - r_asec
		else:
			arcsec = 60.0-(r_asec-asec)
			arcmin -= 1
	return degrees,arcmin,arcsec




def dec_min_negative(deg,amin,asec):
# Amin and asec are positive, deg is negative!
	r_amin = 7.0
	r_asec = 30.0 # = 7.5 arcmin
	if asec + r_asec < 60.0:
		arcsec = asec + r_asec
	else: 
		arcsec = (asec + r_asec) - 60.0
		amin += 1
	if amin + r_amin < 60.0: 
		arcmin = amin + r_amin
		degrees = deg
	else: 
		arcmin = (amin + r_amin) - 60.0
		degrees = deg - 1
		
	return degrees,arcmin,arcsec

# OLD
#	if amin-r_amin > 0.0:
#		arcmin = amin-r_amin
#		degrees = deg
#		if asec - r_asec > 0.0:
#			arcsec = asec - r_asec
#		else:
#			arcsec = 60.0-(r_asec-asec)
#			arcmin -= 1
#	else:
#		arcmin = 60.0-(r_amin-amin)
#		degrees = deg - 1
#		if asec - r_asec > 0.0:
#			arcsec = asec - r_asec
#		else:
#			arcsec = 60.0-(r_asec-asec)
#			arcmin -= 1



def dec_max_negative(deg,amin,asec):
# ARCMIN AND ARCSEC ARE POSITIVE, AND RETURNED POSITIVE!!!
	r_amin = 7.0
	r_asec = 30.0 # = 7.5 arcmin
	if asec - r_asec >= 0.0:
		arcsec = asec - r_asec
	else: 
		arcsec = asec + r_asec
		amin -= 1
	if amin - r_amin > 0.0:
		arcmin = amin - r_amin
		degrees = deg
	else:
		arcmin = 60.0 - r_amin + amin
		# degrees/deg are negative! 
		degrees = deg + 1
	return degrees,arcmin,arcsec

# OLD
#	if amin + r_amin < 60.0:
#		arcmin = amin + r_amin
#		degrees = deg
#		if asec+r_asec < 60.0:
#			arcsec = asec+r_asec
#		else: 
#			arcsec = r_asec - (60-asec)
#			arcmin +=1
#	else: 
#		arcmin = r_amin - (60.0-amin) 
#		degrees = deg + 1
#		if asec+r_asec < 60.0:
#			arcsec = asec+r_asec
#		else: 
#			arcsec = r_asec - (60-asec)
#			arcmin +=1


########################################
# Functions: RA and DEC to degrees
########################################
def ra_to_deg(h,m,s):
	#g = ra.split(':')
	#h = float(g[0])
	#m = float(g[1])
	#s = float(g[2])
	h += (m/60.0)
	h += (s/3600.0)
	deg = h*(360.0/24.0)
	return deg

def dec_to_deg(deg,amin,asec):
	y = dec.split(':')
	#deg,amin,asec = float(y[0]),float(y[1]),float(y[2])
	if deg < 0:
		deg -= (amin/60.0)
		deg -= (asec/3600.0)
		return round(deg,6)
	deg += (amin/60.0)
	deg += (asec/3600.0)
	return deg


########################################
# Path to outputs
########################################
outPath = '/home/emma/gradu/data/All_data/catalogs/Akari/'

def outFile_FIS(filename):
	o = filename.split('_')
	outName = o[0]+'_FIS.txt'
	return outName

def outFile_PSC(filename):
	o = filename.split('_')
	outName = o[0]+'_PSC.txt'
	return outName
########################################
# Run code
########################################


############################
# Find min and max 
# coordinates for sources
############################
with open(matchDir+'converted_centers.txt') as matches:
	next(matches)
	next(matches)
	for line in matches: 
		m = line.split()
		
		# Name 
		scuba = m[0]

		# ra and dec
		ra = m[1]+':'+m[3]+':'+m[5]
		dec = m[7]
		d = dec.split(':')

		# Individual coordinates
		h,m,s = float(m[1]),float(m[3]),float(m[5])
		deg,amin,asec = float(d[0]),float(d[1]),float(d[2])

		h_max,m_max,s_max = ra_max(h,m,s)
		h_min,m_min,s_min = ra_min(h,m,s)

		if deg >= 0:
			deg_max,amin_max,asec_max = dec_max_positive(deg,amin,asec)
			deg_min,amin_min,asec_min = dec_min_positive(deg,amin,asec)
		else: 
			deg_max,amin_max,asec_max = dec_max_negative(deg,amin,asec)
			deg_min,amin_min,asec_min = dec_min_negative(deg,amin,asec)
		
	#	print 5*'*'
	#	print deg,':',amin,':',asec,'\n'
	#	print deg_max,':',amin_max,':',asec_max
	#	print deg_min,':',amin_min,':',asec_min
	#	print 5*'*'

		# Convert these to degrees

		max_ra = ra_to_deg(h_max,m_max,s_max)
		min_ra = ra_to_deg(h_min,m_min,s_min)
		max_dec = dec_to_deg(deg_max,amin_max,asec_max)
		min_dec = dec_to_deg(deg_min,amin_min,asec_min)

		print 5*'*'
		print max_ra,'	',min_ra
		print max_dec,'	',min_dec

		ra_deg = ra_to_deg(h,m,s)
		
		dec_deg = dec_to_deg(deg,amin,asec)
		############################
		# Open files
		############################	
		file_out_FIS = outFile_FIS(scuba)
		outPutFile_FIS = open(outPath+file_out_FIS,'w')
		
		file_out_PSC = outFile_PSC(scuba)
		outPutFile_PSC = open(outPath+file_out_PSC,'w')

		############################
		# Write preliminary info
		# to files
		############################	
		outPutFile_FIS.write(scuba+'\n'+'ra: '+str(ra)+' = '+str(ra_deg)+'\n')
		outPutFile_FIS.write('dec: '+str(dec)+' = '+str(dec_deg)+'\n')
		outPutFile_FIS.write('max ra: '+str(max_ra)+'\n')
		outPutFile_FIS.write('min ra: '+str(min_ra)+'\n'+'max dec: '+str(max_dec)+'\n')
		outPutFile_FIS.write('min dec: '+str(min_dec)+'\n')
		
		outPutFile_PSC.write(scuba+'\n'+'ra: '+str(ra)+' = '+str(ra_deg)+'\n')
		outPutFile_PSC.write('dec: '+str(dec)+' = '+str(dec_deg)+'\n')
		outPutFile_PSC.write('max ra: '+str(max_ra)+'\n')
		outPutFile_PSC.write('min ra: '+str(min_ra)+'\n'+'max dec: '+str(max_dec)+'\n')
		outPutFile_PSC.write('min dec: '+str(min_dec)+'\n')
		############################
		# Go through Akari catalog
		# and find matches
		############################

		#if scuba != 'G202.31+02.53_850.fits': 
		#	continue	


		with open(Akari_FIS) as tbl:
			next(tbl)
		#	print max_ra,'	',min_ra
			for line in tbl:
				e = line.split()
				RA = float(e[2])
				D = float(e[3])
			#	print max_ra,'	',min_ra,'	',RA
				if RA < max_ra and RA > min_ra:
					if D < max_dec and D > min_dec:
						outPutFile_FIS.write(line)#e[0]+'\n')

	
		with open(Akari_PSC) as psc:
			next(psc)
		#	print max_ra,'	',min_ra
			for line in psc:
				e = line.split()
				RA = float(e[2])
				D = float(e[3])
			#	print max_ra,'	',min_ra,'	',RA
				if RA < max_ra and RA > min_ra:
					if D < max_dec and D > min_dec:
						outPutFile_PSC.write(line)#e[0]+'\n')



		continue


























	

		
