########################################################################
#		Same as subtract.py but or I250, T and tau maps
########################################################################
#				Calculating B - C
#			Where B = 250 + 350 um map
#		And C = B convolved to resolution of 500 um map
########################################################################
#			FIRST reproject A maps to 350 um
########################################################################
#			Then calculating A + (B - C)
#				Aug 7, 2018
########################################################################
import numpy as np
from astropy.io import fits
import math
import os
import montage_wrapper as montage
import scipy
import scipy.ndimage as nd

from astropy.utils.data import get_pkg_data_filename
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy.convolution import Gaussian1DKernel
######################################
#	Paths
######################################


Adir = '/home/emma/gradu/data/All_data/n_tau_all/n_tau_maps/reprojected_to_350/'
Bdir = '/home/emma/gradu/data/All_data/n_tau_all/two_map_N_tau/'
Cdir = '/home/emma/gradu/data/All_data/n_tau_all/two_map_N_tau/convolved_to_500/'


Apath = Adir
Bpath = Bdir
Cpath = Cdir

A_B_C_out = '/home/emma/gradu/data/All_data/combinations/A_plus_B_min_C/'
B_C_out = '/home/emma/gradu/data/All_data/n_tau_all/two_map_N_tau/B_minus_C/'


######################################
#	Functions
######################################
def sigma(img,band):
	sigma350 = 35.4
	if band == 350:
		sigmaU = 24.2
	minus = sigma350**2 - sigmaU**2
	sq = math.sqrt(minus)
	
	CDELT1_deg = img.header['CDELT1']
	CDELT1 = abs(CDELT1_deg*3600.0)
	print 'C= ',CDELT1
	new_FWHM = sq / CDELT1
	print 'new FWHM= ', new_FWHM
	sigma = new_FWHM/(math.sqrt(CDELT1*8*math.log(2)))
	print sigma
	return sigma

def convolve_map(A, fwhm, mode='reflect', radius=3.6):
    # fwhm = fwhm of the beam,  in pixels
    # make the kernel with up to ?? fwhm
    n      = int(2+fwhm*radius)
    n = min(n, A.shape[0]-1)
    if (n%2==0):
        n += 1   # make sure it is odd (?)
    kernel = np.zeros((n,n), np.float32)
    K      = 4.0*math.log(2.0)/(fwhm*fwhm) 
    D2     = max(1.7, (1.5*fwhm))**2.0   # outer radius in pixels
    for i in range(n):                   # loop over kernel pixels
        dx  = (n-1.0)/2.0 - i
        for j in range(n):
            dy  = (n-1.0)/2.0 - j
            d2  = dx*dx+dy*dy
            if (d2<D2):
                kernel[i,j] =  np.exp(-K*d2)
    # normalize: integral over kernel equals one
    kernel /= sum(kernel.flat)
    tmp     = np.asarray(nd.convolve(np.asarray(A,np.float32), kernel, mode=mode), np.float32)
    return tmp
######################################
#	Convolve B to 500 um
######################################
for filename in os.listdir(Bdir):
	if filename.endswith('.fits'):

		break

		#if filename != 'G202.16+02.64_test_T.fits':
		#	continue

		p = filename.split('_')
		#print p[2]
		
		field = p[0]
		
		if p[1] == 'test':
			num = '1'
			typ=p[2].strip('.fits')
		else: 
			num = p[1]
			typ=p[3].strip('.fits')
		
		
		if typ == 'BETA':
			continue
		#elif typ == 'NH2':
		#	continue
		elif typ == 'au':
			typ = 'tau'


		outF = field+'_'+typ+'_conv_'+num+'.fits'

		print outF



		f = fits.open(Bdir+filename)[0]

		fData = f.data

		sKernel = sigma(f,350)

		kernel = Gaussian2DKernel(sKernel)
		
		conv = convolve(fData,kernel,normalize_kernel=True)#,boundary='wrap'
		#conv2 = convolve_map(fData,0.346956011926)

		f.data = conv

		f.writeto(Cdir+outF)
		#f.writeto('test.fits',clobber=True)

#C=  10.0
#new FWHM=  2.58364084191
#0.346956011926

######################################
#	B - C
######################################
missedList = []



for filename in os.listdir(Cdir):
	if filename.endswith('.fits'):


		break

		c = fits.open(Cpath+str(filename))[0]
		cData = c.data
		p = filename.split('_')

		field = p[0]

		if filename.endswith('conv_1.fits'):
			typ=p[1]
			bFile = p[0] + '_test_'+p[1]+'.fits'
			outPutFile = B_C_out + p[0] + '_'+typ+'_B_C_1.fits'
		else: 
			typ=p[1]
			num=p[3].split('.')[0]
			bFile = p[0] + '_'+num+'_test_'+typ+'.fits'
			outPutFile = B_C_out + p[0] + '_'+typ+ '_B_C_'+num+'.fits'

			
		
		b = fits.open(Bpath+bFile)[0]
		bData = b.data
		minus = bData - cData
		c.data = minus
		c.writeto(outPutFile,clobber=True)
		





########################################################################
#			Calculating A + (B - C)
#		where A = the 250, 350, and 500 um map
#			And (B - C) is calculated above
########################################################################
resList = []
#for filename in os.listdir(B_C_out):
#	if filename.endswith('res_1.fits'):
#		resList.append(filename)
#	elif filename.endswith('res_2.fits'):
#		resList.append(filename)
#	elif filename.endswith('res_3.fits'):
#		resList.append(filename)
#	elif filename.endswith('res_4.fits'):
#		resList.append(filename)
#	elif filename.endswith('res_5.fits'):
#		resList.append(filename)
#	else: 
#		print 'uhh'

i = 0

#G006.04+36.73_I250_B_C_1.fits	BCfile
# G006.04+36.73_I250_rep_1.fits	AFile

for filename in os.listdir(B_C_out):
	if filename.endswith('.fits'):
		#break
		bc = fits.open(B_C_out + str(filename))[0]
		bcData = bc.data
		p = filename.split('_')

		field=p[0]
		typ=p[1]
		num=p[4].split('.')[0]

		print filename,num

		if typ != 'NH2':
			aFile = field+'_'+typ+'_rep_'+num+'.fits'
			outPut = field+'_'+typ+'_ABC_'+num+'.fits'
		else: 
			aFile = 'NH2/'+field+'_rep_'+num+'.fits'
			outPut = field+'_'+typ+'_ABC_'+num+'.fits'
		
		#if filename.endswith('_1.fits'):
		#	aFile = p[0] + '_rep_1.fits'
		#	outPut = p[0] + '_ABC_1.fits'
		#elif filename.endswith('_2.fits'):
		#	aFile = p[0] + '_rep_2.fits'
		#	outPut = p[0] +	'_ABC_2.fits'		
		#elif filename.endswith('_3.fits'):
		#	aFile = p[0] + '_rep_3.fits'
		#	outPut = p[0] +	'_ABC_3.fits'		
		#elif filename.endswith('_4.fits'):
		#	aFile = p[0] + '_rep_4.fits'
		#	outPut = p[0] +	'_ABC_4.fits'		
		#elif filename.endswith('_5.fits'):
		#	aFile = p[0] + '_rep_5.fits'
		#	outPut = p[0] +	'_ABC_5.fits'		
		#else:
		#	print 'meow',filename
		#	i += 1
		#	continue

		a = fits.open(Apath+aFile)[0]
		aData = a.data
		plus = aData + bcData
		bc.data = plus
		bc.writeto(A_B_C_out + outPut,clobber=True)
		continue



#print i # i counts number of missed images
if i == 0:
	print 'jee'






