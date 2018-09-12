##################################################################################
#			Reproject the A(500-160) and B-C maps
#			to the resolution of the 250 um map
##################################################################################
#				Jul 4, 2018
##################################################################################

from astropy.io import fits
import os
import montage_wrapper as montage
###
from astropy.utils.data import get_pkg_data_filename
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy.convolution import Gaussian1DKernel


########################################
# Paths
########################################
# A
A_in = '/home/emma/gradu/data/All_data/n_tau_all/500_160_n_tau/'
A_out = '/home/emma/gradu/data/All_data/n_tau_all/500_160_n_tau/reproject_to_250/'

# B
BC_in = '/home/emma/gradu/data/All_data/n_tau_all/350_160_n_tau/B_minus_C/'
BC_out = '/home/emma/gradu/data/All_data/n_tau_all/350_160_n_tau/B_minus_C/reproject_to_250/'

# headers
hdr_in = '/home/emma/gradu/data/All_data/SPIRE_PSW_250_Herschel/headers_250/'


########################################
# Functions
########################################
def reproject(filename,path):
	print 'meow'

def check(path,filename):
	f = fits.open(path+filename)[0]
	if f.data.shape != (121,121):
		print filename, '	',f.data.shape
	else:
		print filename


########################################
# Convolve
########################################
def sigma(img,band):
	if band == 500:
		sigmaU = 35.4		# A FWHM
	elif band == 350:
		sigmaU = 24.2		# BC FWHM
	sigma250 = 17.9		# Result FWHM
	
	minus = sigmaU**2 - sigma250**2
	sq = math.sqrt(minus)
	
	CDELT1_deg = img.header['CDELT1']
	CDELT1 = abs(CDELT1_deg*3600.0)
	print 'C= ',CDELT1
	new_FWHM = sq / CDELT1
	sigma = new_FWHM/(math.sqrt(CDELT1*math.log(2)))
#	sigma = new_FWHM/(math.sqrt(CDELT1*8*math.log(2)))
	return sigma


#def convolve(filename):
#	current = get_pkg_data_filename(filename)
#	img_2 = fits.open(current)[0]
#	img = img_2.data
#	sKernel2 = sigma(img_2)
#	# Create kernel and convolve
#	kernel = Gaussian2DKernel(sKernel2)
#	conv = convolve(img,kernel)
#	# Write to file
#	img_2.data = conv
#	return img_2



print 5*'MEOW'

########################################
# Reproject
########################################
L25list = ['G017.69-00.15','G070.07-01.60','G094.11-04.86','G158.55-20.92','G178.31+00.28']
for filename in os.listdir(BC_in):
	if filename.endswith('.fits'):
		inFile = BC_in+filename
		p = filename.split('_')
		l = 'L2.0'
		for i in range(0,len(L25list)):
			now = L25list[i]
			if p[0] == now:
				l = 'L2.5'
				break
		if filename.endswith('1.fits'):
			hdr = hdr_in+p[0]+'_PSW_'+l+'.txt'
			outFile = BC_out+filename.strip('_B_C_1.fits')+('_rep_BC_1.fits')
			AFile = A_in+p[0]+'_test_NH2.fits'
			A_rep = A_out+p[0]+'_rep_1.fits'
		elif filename.endswith('2.fits'):
			hdr = hdr_in+p[0]+'_PSW_'+l+'_2.txt'
			outFile = BC_out+filename.strip('_B_C_2.fits')+('_rep_BC_2.fits')
			AFile = A_in+p[0]+'_2_test_NH2.fits'
			A_rep = A_out+p[0]+'_rep_2.fits'
		else: 
			print '!!!!!!!!!!!!!!!!!!!!!!!!!'
			print filename
			print 'AAAAA PANIC'
			print '!!!!!!!!!!!!!!!!!!!!!!!!!'

	#	print filename,'\n',hdr
	#	print AFile,'\n',A_rep
	#	print 10*'*'
		montage.reproject(inFile,outFile,header=hdr,exact_size=True)
		montage.reproject(AFile,A_rep,header=hdr,exact_size=True)


for filename in os.listdir(BC_out):
	if filename.endswith('.fits'):
		check(BC_out,filename)

for filename in os.listdir(A_out):
	if filename.endswith('.fits'):
		check(A_out,filename)







