################################################################################
#		Take results from extractcllumps and make them readable
################################################################################

import os

#################################
# Paths
#################################
dir_in = '/home/emma/gradu/data/All_data/SCUBA_850/clumps/extractResults/results/'
dir_out = '/home/emma/gradu/data/All_data/SCUBA_850/clumps/extractResults/results/resultTables/'

#################################
# Functions
#################################
def deg_to_ra(deg): 
	H = deg*(24.0/360)
	h = int(H)
	re = H-h
	M = re*60.0
	m = int(M)
	re = M-m
	S = re*60.0
	s = round(S,2)
	string = str(h)+':'+str(m)+':'+str(s)
	return h,m,s,string
	

def deg_to_dec(deg):
	d = int(deg)
	re = deg-d
	AMIN = re*60.0
	amin = int(AMIN)
	re = AMIN-amin
	ASEC = re*60.0
	asec = round(ASEC,2)
	string = str(d)+':'+str(abs(amin))+':'+str(abs(asec))
	return d,abs(amin),abs(asec),string




#################################
# Code
#################################
for filename in os.listdir(dir_in):
	if filename.endswith('clump_info.txt'):
	#if filename.endswith('G202_31+02_53_850_clump_info.txt'):
		outName = filename.strip('_clump_info.txt')+'_clumps.txt'
		out = open(dir_out+outName,'w')
		out.write('Clump ID	peak[deg]			center[deg]			size[arcsec]			sum[mJy/arcsec**2]	peak[mJy/arcsec**2]	vol[arcsec.arcsec]\n')

		
		with open(dir_in+filename) as t: 
			next(t)
			for line in t: 
				p = line.split()
				clumpNum = p[0].split('(')[1].strip(',')

		
				ra_peak = float(p[1].strip(','))
				dec_peak = float(p[2].strip(','))
		
				ra_cen = float(p[3].strip(','))
				dec_cen = float(p[4].strip(','))
		
				ra_size = float(p[5].strip(','))
				dec_size = float(p[6].strip(','))

		
				SUM = float(p[7].strip(','))
			
				peak = float(p[8].strip(','))
		
				volume = float(p[9].strip(','))

				shape = p[10].split("'")[1]+'_'+p[11]+'_'+p[12]+'_'+p[13]+'_'+p[14]+'_'+p[15]+'_'+p[16]+'_'+p[17].strip("')")
		
		
				#h_p,m_p,s_p,ra_p = deg_to_ra(ra_peak)
				#d_p,amin_p,asec_p,dec_p = deg_to_dec(dec_peak)
		
				#h_c,m_c,s_c,ra_c = deg_to_ra(ra_cen)
				#d_c,amin_c,asec_c,dec_c = deg_to_dec(dec_cen)

				#h_p,m_p,s_p,ra_p = deg_to_ra(ra_peak)
				#d_p,amin_p,asec_p,dec_p = deg_to_dec(dec_peak)

				print ra_peak,dec_peak
		
				print 10*'*'
		
		
				out.write(clumpNum+'		'+'('+str(ra_peak)+','+str(dec_peak)+')')
				out.write('	'+'('+str(ra_cen)+','+str(dec_cen)+')')
				out.write('	'+'('+str(ra_size)+','+str(dec_size)+')')
				out.write('	'+str(SUM)+'		'+str(peak)+'		'+str(volume)+'\n')


		with open(dir_in+filename) as t: 
			next(t)
			out.write(60*'#'+'\n')
			out.write('Clump ID	ra[deg]		dec[deg]	ra size[deg]	dec size[deg]	position angle	type\n')
			for line in t: 
				p = line.split()
				clumpNum = p[0].split('(')[1].strip(',')

				ra_c = p[13]
				dec_c = p[14]
		
				d_ra = p[15]
				d_dec = p[16]
		
				angle = p[17].strip("')")
		
				out.write(clumpNum+'		'+ra_c+'	'+dec_c+'	'+d_ra+'	'+d_dec+'	'+angle)
				out.write('	'+p[10].split("'")[1]+'	'+p[11]+'	'+p[12]+'\n')


		
#################################
# OLD
#################################


		














		
