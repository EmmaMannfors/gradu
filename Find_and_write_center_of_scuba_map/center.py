####################################################################
# For calculating center of SCUBA-2 maps
####################################################################

from astropy.io import fits
import montage_wrapper as montage
import os

####################################################################
#open files
####################################################################
SCUBAList = []
#HerschelList = []
SCUBA = open('SCUBA.txt','r')
#Her = open('Herschel.txt','r')
for line in SCUBA:
	line = line.strip('\n')
	SCUBAList.append(line)

#for line in Her:
#	line = line.strip('\n')
#	HerschelList.append(line)

####################################################################
#From SCUBA file, calculate center of map
####################################################################
i = 0
ra_center_list = []
dec_center_list = []

while i < len(SCUBAList):
	SCUBA_fits = fits.open(SCUBAList[i])
	
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

	ra_center_list.append(ra_center)
	dec_center_list.append(dec_center)
	
	i +=1


n = 0

center = open('centers.txt','w')
center.write('Centers of SCUBA-2 maps in degrees \n')
center.write('RA                       dec \n')
while n < len(ra_center_list):
	ra = str(ra_center_list[n])
	dec = str(dec_center_list[n])
	center.write(ra + '     ' + dec + '\n')
	n += 1








