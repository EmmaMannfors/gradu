####################################################################
#Script to resize one Herschel image to the size of one SCUBA image
####################################################################

from astropy.io import fits
import montage_wrapper as montage
import os

####################################################################
#open files
####################################################################
SCUBAList = []
HerschelList = []

# SCUBA 850 um and Herschel 250 um
#SCUBA = open('SCUBA.txt','r')
#Her = open('Herschel.txt','r')
#for line in SCUBA:
#	line = line.strip('\n')
#	SCUBAList.append(line)

#for line in Her:
#	line = line.strip('\n')
#	HerschelList.append(line)

#SCUBA 850 um and Herschel 160 um
#f = open('/home/emma/gradu/data/All data/pub/PACS_160_Herschel/SCUBA850_Herschel160.txt','r')

with open('/home/emma/gradu/data/All data/pub/PACS_160_Herschel/SCUBA850_Herschel160.txt') as f:
	next(f)
	for line in f: 
		p = line.split()
		SCUBAList.append(str(p[0]))
		HerschelList.append(str(p[1]))

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


####################################################################
#make sure NAXIS >= 2 in Herschel file
####################################################################

i = 0

NAXISList = []

while i < len(HerschelList):
	outPutFile = HerschelList[i].strip('.fits')+'_NAXIS_2.fits'
	data,header = fits.getdata(HerschelList[i],header=True)
	header['NAXIS'] = 2
	if not os.path.exists(outPutFile):
		fits.writeto(outPutFile,data,header)
	NAXISList.append(outPutFile)	
	i+=1



####################################################################
#Crop image
####################################################################
i = 0
cropList = []

while i < len(NAXISList):
	#print NAXISList[i]
	currentFile = str(NAXISList[i])
#These are for Herschel SPIRE 250 um!
#	if i == 4 or i == 9 or i == 17:
#		outputFile = currentFile.strip('_NAXIS_2.fits')+'_crop_2.fits'
#	if i == 8:
#		outputFile = currentFile.strip('_NAXIS_2.fits')+'_crop_3.fits'
#	else:
#		outputFile = currentFile.strip('_NAXIS_2.fits')+'_crop.fits'
	outputFile = currentFile.strip('_NAXIS_2.fits')+'_crop.fits'
	montage.mSubimage(currentFile,outputFile,ra_center_list[i],dec_center_list[i],0.2)
	cropList.append(outputFile)
	i += 1
	



i = 0
#finalFile = open('Final_sources.txt','w')
#finalFile.write('SCUBA                     Herschel\n')

#850 & 160 um
finalFile = open('Final_sources_160.txt','w')
finalFile.write('SCUBA                     Herschel 160um\n')

while i < len(cropList):
	scubaNow = SCUBAList[i]
	herNow = cropList[i]
	finalFile.write(scubaNow + '    ' + herNow + '\n')
	#finalFile.write('meow')
	i+=1
	

finalFile.close()








