####################################################################
# For calculating center of SCUBA-2 maps
####################################################################

from astropy.io import fits
import montage_wrapper as montage
import os

####################################################################
#	Path
####################################################################

scuba_in = '/home/emma/gradu/data/All_data/SCUBA_850/'



####################################################################
#	Convert to arcsec and stuff
####################################################################
#Create the units to convert to
hours = 0.0
minutes = 0.0
seconds = 0.0



def truncate(f, n): #from the internetz
    '''Truncates/pads a float f to n decimal places without rounding'''
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return '.'.join([i, (d+'0'*n)[:n]])


i = 0

def ra_to_h(deg):
	hours = deg*(24.0/360)
	current = hours
	hours = int(hours)
	current -= hours
	minutes = current * 60.0
	current = minutes
	minutes = int(minutes)
	current -= minutes
	seconds = current * 60.0
	seconds = truncate(seconds,3)
	ans = str(hours)+" h, "+str(minutes)+" m, "+str(seconds)+" s"
	return ans

deg = 0.0
arcmin = 0.0
arcsec = 0.0


def deg_to_arcsec(x):
	deg = int(x)  	#x is already in degrees
	current = x - deg
	current = abs(current)   #To get rid of negative values
	arcmin = 60*current
	current = arcmin
	arcmin = int(arcmin)
	current -= arcmin
	arcsec = 60*current
	arcsec = truncate(arcsec,3)
	ans = str(deg)+":"+str(arcmin)+":"+str(arcsec)
	return ans
	


####################################################################
#From SCUBA file, calculate center of map
####################################################################
i = 0
ra_center_list = []
dec_center_list = []


deg = open('centers.txt','w')
deg.write('Centers of SCUBA-2 maps in degrees \n')
deg.write('RA                       dec \n')

conv = open('converted_centers.txt','w')
conv.write('Centers of SCUBA-2 maps in degrees \n')
conv.write('Field			RA                       dec \n')


ans = raw_input('LaTex-compatible table? ')
tex = open('latex_scuba_center.txt','w')


for filename in os.listdir(scuba_in):
	if filename.endswith('.fits'):	
		SCUBA_fits = fits.open(scuba_in+filename)
		
		b_left_ra = SCUBA_fits[0].header['OBSRABL']
		b_left_dec = SCUBA_fits[0].header['OBSDECBL']
		t_left_ra = SCUBA_fits[0].header['OBSRATL']
		t_left_dec = SCUBA_fits[0].header['OBSDECTL']
		b_right_ra = SCUBA_fits[0].header['OBSRABR']
		b_right_dec = SCUBA_fits[0].header['OBSDECBR']
		t_right_ra = SCUBA_fits[0].header['OBSRATR']	
		t_right_dec = SCUBA_fits[0].header['OBSDECTR']

		left_ra_avg = (b_left_ra + t_left_ra)/2.0
		right_ra_avg = (b_right_ra + t_right_ra)/2.0
	
		t_dec_avg = (t_right_dec + t_left_dec)/2.0
		b_dec_avg = (b_right_dec + b_left_dec)/2.0
	
		if left_ra_avg > right_ra_avg:
			ra_center = ((left_ra_avg - right_ra_avg)/2.0) + right_ra_avg
		else:
			ra_center = ((right_ra_avg - left_ra_avg)/2.0) + left_ra_avg
	
		if t_dec_avg > b_dec_avg:
			dec_center = ((t_dec_avg - b_dec_avg)/2.0) + b_dec_avg
		else: 
			dec_center = ((b_dec_avg - t_dec_avg)/2.0) + t_dec_avg
	
#		ra_center_list.append(ra_center)
#		dec_center_list.append(dec_center)


		print 'meow'

		ra_conv = ra_to_h(ra_center)
		dec_conv = deg_to_arcsec(dec_center)

		conv.write(filename+'	'+ra_conv+'	'+dec_conv+'\n')
	

		if ans == 'y':
			#do Latex compatible table
			e = filename.split('_')
			tex.write('		\hline'+'\n')
			tex.write('		'+e[0]+' & '+ra_conv+' & '+dec_conv+'\\\\'+'\n')

	
		continue
	else: 
		continue



if ans == 'y':
	tex.write('		\hline')





#while n < len(ra_center_list):
#	ra = str(ra_center_list[n])
#	dec = str(dec_center_list[n])
#	center.write(ra + '     ' + dec + '\n')
#	n += 1





	





