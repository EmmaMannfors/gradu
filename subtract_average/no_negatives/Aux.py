import numpy, matplotlib
from   matplotlib.pylab import *

from   scipy.optimize import fsolve, fmin, bisect, fmin_powell, leastsq
from   scipy.interpolate import interp1d, bisplrep, bisplev, splrep, splev, splder
from   scipy.integrate import quad
from   numpy import array, NaN, nonzero, log10, concatenate, exp, dot
from   scipy.ndimage import rotate, zoom
# from   MJ.Aux.ColourCorrect import *
# from   MJ.Aux.FITS import CopyEmptyFits
import math
import pylab
import os
import matplotlib
import time
import multiprocessing



#########################################################################################
##############################  From MJ/mjDefs.py  ######################################
#########################################################################################

try:
    from astropy.io import fits as pyfits
except:
    import pyfits


# the following (mostly) in cgs system !!
C_LIGHT          = 2.99792458E10
C_LIGHT_SI       = 2.99792458E8
AMU              = 1.6605E-24 
H_K              = 4.7995074E-11 
#BOLZMANN        = 1.3806E-16
BOLZMANN         = 1.3806488e-16
BOLTZMANN        = 1.3806488e-16
# BOLTZMANN_SI     = 1.3806E-23  
BOLTZMANN_SI     = 1.3806488e-23
STEFAN_BOLTZMANN = 5.670373e-5
SB_SI            = 5.670373e-8
CGS_TO_JY_SR     = 1e23          # erg/cm2/sr/Hz = CGS_TO_JY_SR * Jy/sr
# PLANCK         = 6.6262E-27
PLANCK           = 6.62606957e-27 
# PLANCK_SI      = 6.6262E-34
PLANCK_SI        = 6.62606957e-34
M0               = 1.99e33
MJupiter         = 1.9e30        # [g]
GRAV             = 6.67e-8
GRAV_SI          = 6.673e-11
PARSEC           = 3.0857E18
ELECTRONVOLT     = 1.6022e-12
AU               = 149.597871e11
RSUN             = 6.955e10
RSUN_SI          = 6.955e8
DSUN             = 1.496e13  # cm
DSUN_SI          = 1.496e11  # 1.496e8 km
MSUN             = 1.9891e33
MSUN_SI          = 1.9891e30
M_EARTH          = 5.972e27
LSUN             = 3.839e33
LSUN_SI          = 3.839e26
TSUN             = 5778.0
MJUPITER         = 1.9e30
H_C2             = PLANCK/(C_LIGHT*C_LIGHT)
H_C2_GHz         = PLANCK/(C_LIGHT*C_LIGHT)*1.0e27

ARCSEC_TO_DEGREE =  (1.0/3600.0)
DEGREE_TO_RADIAN =  0.0174532925199432958
ARCMIN_TO_RADIAN =  (2.9088820e-4)
ARCSEC_TO_RADIAN =  (4.8481368e-6)
HOUR_TO_RADIAN   =  (0.261799387)
MINUTE_TO_RADIAN =  (4.3633231e-3)
SECOND_TO_RADIAN =  (7.2722052e-5)

RADIAN_TO_DEGREE =  57.2957795130823208768
RADIAN_TO_ARCMIN =  3437.746771
RADIAN_TO_ARCSEC =  206264.8063
RADIAN_TO_HOUR   =  3.819718634
RADIAN_TO_MINUTE =  229.1831181
RADIAN_TO_SECOND =  13750.98708

ARCMIN_TO_DEGREE =   (1.0/60.0)
DEGREE_TO_ARCMIN =   60.0
DEGREE_TO_ARCSEC =   3600.0
BSD              =   0.94e21    #  N(H2)=0.94e21*Av

kB = 1024.0
MB = 1024.0*kB
GB = 1024.0*MB

# Rayleigh = 1e6/(4*pi)  photons / cm2 / s / sr = 2.41e-7 erg/cm2/sr/s 
# H-alpha 6563 Angstrom => 4.5679e+14                                  
# 1 km/s = 1.5237e+9 Hz                                                
# R/(km/s) = 2.41e-7 / 1.52e9 = 1.58e-16  erg/cm2/s/sr/Hz              
RAYLEIGH_KMS_TO_JY_SR = 1.5807943647067913e-16

def HMS2RAD(h, m, s):
    if (h>=0): 
        return  (h+(1.0*m)/60.0+(1.0*s)/3600.0)*HOUR_TO_RADIAN
    else:    # someone uses negative hours ?
        return  -(-(1.0*h)+(1.0*m)/60.0+(1.0*s)/3600.0)*HOUR_TO_RADIAN
    
def DMS2RAD(d, m, s):
    # If  -1<DEC<0.0, d must be float so that we have it as -0.0
    # ok as long as -0.0 has STRING REPRESENTATION WITH MINUS SIGN INCLUDED
    if (d==0.0):
        sss = '%.1f' % d
        if (sss.find('-')>=0):  # d was -0.0 !!
            return -(-(1.0*d)+(1.0*m)/60.0+(1.0*s)/3600.0)*DEGREE_TO_RADIAN            
    if (d>=0):
        return ((1.0*d)+(1.0*m)/60.0+(1.0*s)/3600.0)*DEGREE_TO_RADIAN
    else:
        return -(-(1.0*d)+(1.0*m)/60.0+(1.0*s)/3600.0)*DEGREE_TO_RADIAN

def HMS2DEG(h,m,s):
    return RADIAN_TO_DEGREE*HMS2RAD(h,m,s)

def DMS2DEG(d, m, s):
    # Note -- if -1.0<DEC<0.0, d must be ***float*** -0.0 so that we get the sign right
    return RADIAN_TO_DEGREE*DMS2RAD(d,m,s)
    
def RAD2HMS(x):
    xx  =  abs(x*RADIAN_TO_SECOND)
    h   =  int(xx/3600)
    xx -=  h*3600.0 
    h   =  h % 24
    m   =  int(xx/60) 
    s   =  xx-m*60.0 
    if (abs(s-60.0)<0.1):
        s  = 0.0 
        m += 1 
        if (m==60):
            m  = 0 
            h += 1 
        if (h==24):
            h  = 0
    return (h, m, s)
    
def RAD2DMS(x):
    xx  =  abs(x*RADIAN_TO_ARCSEC) 
    d   =  int(xx/3600)
    xx -=  d*3600.0 
    m   =  int(xx/60) 
    s   =  xx-m*60.0 
    if (abs(s-60.0)<0.01):
        s  = 0.0 
        m += 1 
        if (m==60):
            m  = 0 
            d += 1
    if (x<0.0):
        d = -d 
    return (d, m, s)

def DEG2HMS(x):
    return RAD2HMS(x*DEGREE_TO_RADIAN)

def DEG2DMS(x):
    return RAD2DMS(x*DEGREE_TO_RADIAN)

def STR_HMS2RAD(s):
    """
    Convert string 'h m s' into radians
    """
    ss = s.lower().replace(':', ' ')
    ss = ss.replace('h',' ').replace('m', ' ').replace('s', ' ').replace(':', ' ').split()
    h, m, s = float(ss[0]), float(ss[1]), float(ss[2])
    return HMS2RAD(h,m,s)

def STR_DMS2RAD(s):
    """
    Convert string 'd am as' into radians    
    """
    ss = s.lower().replace(':', ' ')
    ss = ss.replace('d',' ').replace('m', ' ').replace('s', ' ').replace('\'',' ').replace('"', ' ').replace(':', ' ').split()
    d, am, ars = float(ss[0]), float(ss[1]), float(ss[2])
    # if -1.0<DEC<0.0, d=-0.0 => sign should be correctly interpreted in DMS2RAD
    res = DMS2RAD(d,am,ars)
    if (ss[0].find('-')>=0): 
        if (res>0.0):
            print('*'*80)
            print('ERROR IN STR_DMS2RAD !!!!   %s --> %.7f' % (s, res))
            print('*'*80)
        res = -abs(res)
    return res


def  STR_HMS2DEG(s):
    return RADIAN_TO_DEGREE*STR_HMS2RAD(s)

def STR_DMS2DEG(s):
    return RADIAN_TO_DEGREE*STR_DMS2RAD(s)


def RAD2STRS(x, y, digits=[1,1], colons=False):
    # digits = decimal places in RA and DEC strings
    h, m, s     = RAD2HMS(x)
    dd, mm, ss  = RAD2DMS(y)
    if (digits[0]==0):   
        str_s = '%02d' % rint(s)
        if (str_s=='60'):
            str_s = '00'
            m += 1
            if (m==60):
                h, m = (h+1)%24, 0
    else:
        fmt    = '%0' + '%d.%df' % (digits[0]+3, digits[0])
        str_s = fmt % s
        if (str_s[0:2]=='60'):
            str_s = fmt % 0.0
            m += 1
            if (m==60):
                h, m = (h+1)%24, 0
    if (digits[1]==0):   
        str_ss = '%02d' % rint(ss)
        if (str_ss=='60'):
            str_ss = '00'
            mm += 1
            if (mm==60):
                dd, mm = (dd+1)%360, 0
    else:
        fmt    = '%0' + '%d.%df' % (digits[1]+3, digits[1])
        str_ss = fmt % ss
        if (str_ss[0:2]=='60'):
            str_ss = fmt % 0.0
            mm += 1
            if (mm==60):
                dd, mm = (dd+1)%360, 0
    #
    if (colons):
        return '%02d:%02d:%s' % (h, m, str_s), '%+03d:%02d:%s' % (dd, mm, str_ss)
    return '%02d %02d %s' % (h, m, str_s), '%+03d %02d %s' % (dd, mm, str_ss)



def DEG2STRS(x, y, digits=[1,1], colons=False):
    return RAD2STRS(x*DEGREE_TO_RADIAN, y*DEGREE_TO_RADIAN, digits=[1,1], colons=colons)


def coo_delta(h1, m1, s1, d1, am1, as1, h2, m2, s2, d2, am2, as2):
    ra1 =  HMS2RAD(h1, m1, s1)
    ra2 =  HMS2RAD(h2, m2, s2)
    dec1=  DMS2RAD(d1, am1, as1)
    dec2=  DMS2RAD(d2, am2, as2)
    cosy=  cos(mean([dec1,dec2]))
    print('CENTRE %8.4f %8.4f degrees\n' % ( 0.5*(ra1+ra2)*RADIAN_TO_DEGREE, 0.5*(dec1+dec2)*RADIAN_TO_DEGREE ))
    print('SIZE   %8.4f %8.4f degrees\n' % ( cosy*(ra1-ra2)*RADIAN_TO_DEGREE, (dec1-dec2)*RADIAN_TO_DEGREE ))
    (h, m, s)   = RAD2HMS(0.5*(ra1+ra2))
    (d, am, ars) = RAD2DMS(0.5*(dec1+dec2))
    print('CENTRE %3d %2d %4.1f  %3d %2d %4.1f\n' % ( h, m, s, d, am, ars ))
    
    
    
def planck_function(T, freq, derivatives=False):
    """
    Given a temperature T [K] and frequency freq [Hz], return the corresponding
    blackbody intensity in cgs units [erg/s/cm2/sr/Hz]. Optionally, also return derivatives
    wrt temperature. The precision set explicitly to float32.
    Uses _lots_ of memory (use pyx_planck_function from pyx.pyxPSM if memory is an issue)
    """
    # make T and freq into arrays of similar dimension
    if (0):
        T    = asarray(T, ndmin=1)
        freq = asarray(freq, ndmin=1)
    else:
        if (isscalar(T)):     T    = asarray([T,], float32)
        if (isscalar(freq)):  freq = asarray([freq,], float32)
    if    (len(T)>len(freq)):  freq = freq*ones(T.shape, float32) # freq must have been scalar !
    elif  (len(T)<len(freq)):  T    = T*ones(freq.shape, float32) # T must have been scalar !
    tmp  =  H_K*freq/T                        # argument of the exponential function
    S    =  zeros(len(tmp), float32)
    if (derivatives): dS_dT = S.copy()
    m    =  nonzero((tmp<100.0)&(tmp>1.0e-4)) # _outside_ Rayleigh-Jeans-regime (tmp>100 => result 0.0)
    if (len(m)>0):
        S[m] =  (2.0*H_C2_GHz*(freq[m]/1.0e9)**3.0/(exp(tmp[m])-1.0)).astype(float32)
        if (derivatives):
            dS_dT[m] =  (S[m] * (tmp[m]/T[m]) / (1.0-exp(-tmp[m]))).astype(float32)
    m    =  nonzero(tmp<=1.0e-4)              # _inside_ Rayleigh-Jeans-regime
    if (len(m)>0):
        S[m] = (2.0*freq[m]**2.0*BOLTZMANN*T[m]/(C_LIGHT**2.0)).astype(float32)
        if (derivatives):
            dS_dT[m] =  (S[m]/T[m]).astype(float32)
    if (derivatives):  return S, dS_dT
    return S



def um2f(um):
    """
    Convert wavelength [um] to frequency [Hz]
    """
    return C_LIGHT/(um*1.0e-4)

def f2um(f):
    """
    Convert frequency [Hz] to wavelenth [um]
    """
    return 1.0e4*C_LIGHT/f


def intensity_conversion(unit, UNIT, freq=0.0):
    """
    Return conversion coefficient for transformation of surface brightness 
    from units 'unit' to units 'UNIT'. 
    The recognized unit names are listed in the left hand column below.
    'SI'      -   [W/m2/sr/Hz]
    'cgs'     -   [erg/cm2/sr/s/Hz]
    'MJy/sr'  -   [MJy/sr]
    'Jy/sr'   -   [Jy/sr]
    'KCMB'    -   thermodynamic temperature [K]
    'mKCMB'   -   thermodynamic temperature [mK]
    'KRJ'     -   [K]  in Rayleigh-Jeans approximation
    'mKRJ'    -   [mK] in Rayleigh-Jeans approximation
    In the last four cases the optional argument is the frequency [Hz] (scalar).
    """
    k = 1.0
    if (unit!='cgs'):    # Conversion from 'unit' to cgs units
        if   (unit=='Jy/sr'):    k  = 1.0e-23
        elif (unit=='MJy/sr'):   k  = 1.0e-17        
        elif (unit=='SI'):       k  = 1.0e+3
        elif (unit=='mKCMB'):
            B, dB_dT = planck_function(2.726, asarray([freq,]), True)  # arg = frequency [Hz]
            k        = 1.0e-3*dB_dT[0]
        elif (unit=='KCMB'):
            B, dB_dT = planck_function(2.726, asarray([freq,]), True)  # arg = frequency [Hz]
            k        = dB_dT[0]
        elif (unit=='mKRJ'):  
            k        = 2.0e-3*BOLTZMANN*(freq/C_LIGHT)**2.0  # arg = frequency [Hz]
        elif (unit=='KRJ'):  
            k        = 2.0*BOLTZMANN*(freq/C_LIGHT)**2.0     # arg = frequency [Hz]
        else:
            return 0.0 
    if (UNIT!='cgs'):    # Conversion from cgs units to 'UNIT'
        if   (UNIT=='Jy/sr'):    k *= 1.0e+23
        elif (UNIT=='MJy/sr'):   k *= 1.0e+17        
        elif (UNIT=='SI'):       k *= 1.0e-3
        elif (UNIT=='mKCMB'):
            B, dB_dT = planck_function(2.726, asarray([freq,]), True)  # arg = frequency [Hz]
            k       /=  (1.0e-3*dB_dT[0])
        elif (UNIT=='KCMB'):
            B, dB_dT = planck_function(2.726, asarray([freq,]), True)  # arg = frequency [Hz]
            k       /=  dB_dT[0]
        elif (UNIT=='mKRJ'):
            k       /= 2.0e-3*BOLTZMANN*(freq/C_LIGHT)**2.0  # arg = frequency [Hz]
        elif (UNIT=='KRJ'):
            k       /= 2.0*BOLTZMANN*(freq/C_LIGHT)**2.0     # arg = frequency [Hz]
        else:
            return 0.0
    return k


#########################################################################################
################################  From MJ/FITS.py  ######################################
#########################################################################################

def ExtractHDU(F, hdu=0):
    """
    Make a new pyfits object by extracting one hdu from a multi-hdu image.
    Input:
        hdu  =   the number of the header unit to be extracted
    """
    FF            =  pyfits.HDUList(pyfits.PrimaryHDU())    
    FF[0].header  =  F[hdu].header    
    FF[0].data    =  F[hdu].data
    try:
        FF[0].header.remove('XTENSION', ignore_missing=True)
    except:
        print('No XTENSION keyword...')
    FF.verify('silentfix')
    return FF


def CopyFits(F, input_hdu=0):
    hdu      = pyfits.PrimaryHDU(asarray(F[input_hdu].data.copy(), float32))
    hdulist  = pyfits.HDUList([hdu])
    hdulist[0].header = F[input_hdu].header
    hdulist[0].data   = 1.0*F[input_hdu].data
    return hdulist


def CopyEmptyFits(F, input_hdu=0, data=None, plane=-1):
    """
    Make an empty FITS object using F[input_hdu] as model.
    Input:
        F          = pyfits object
        input_hdu  = select the header unit, default is 0
        data       = if not None, copy this to the newly created Fits object
        plane>=0   = assume 3D FITS, taking one plane
    """
    hdu      = pyfits.PrimaryHDU(asarray(F[input_hdu].data.copy(), float32))
    hdulist  = pyfits.HDUList([hdu])
    A        = F[input_hdu].header
    B        = hdulist[0].header
    M, N     = int(A['NAXIS1']), int(A['NAXIS2'])
    B.update(NAXIS1 = M)
    B.update(NAXIS2 = N)
    B.update(CRVAL1 = A['CRVAL1'])
    B.update(CRVAL2 = A['CRVAL2'])
    try:
        ### B.update(CDELT1 = -abs(A['CDELT1']))  # make sure raw image already has correct orientation
        B.update(CDELT1 = A['CDELT1'])   # make sure raw image already has correct orientation
        B.update(CDELT2 = A['CDELT2'])
    except:
        pass
    for tag in ['CD1_1', 'CD2_2', 'CD3_3']:
        try:
            B.update(tag = A[tag])
        except:
            pass
    B.update(CRPIX1 = A['CRPIX1'])
    B.update(CRPIX2 = A['CRPIX2'])
    B.update(CTYPE1 = A['CTYPE1'])
    B.update(CTYPE2 = A['CTYPE2'])
    try:
        B.update(CROTA2 = A['CROTA2'])
    except:
        B.update(CROTA2 = 0.0)
    hdulist[0].data *= 0.0
    if (plane>=0):
        hdulist[0].data = zeros((N,M), float32)
        hdulist[0].data = 1.0*F[0].data[plane,:,:].reshape(N,M)
    if (data!=None):
        hdulist[0].data = data.copy()
    return hdulist




def CopyEmptyFits2(F, input_hdu=0, data=None, plane=-1):
    """
    Make an empty FITS object using F[input_hdu] as model.
    This one copies the entire header.
    Input:
        F          = pyfits object
        input_hdu  = select the header unit, default is 0
        data       = if not None, copy this to the newly created Fits object
        plane>=0   = assume 3D FITS, taking one plane
    """
    hdu      = pyfits.PrimaryHDU(asarray(F[input_hdu].data.copy(), float32))
    hdulist  = pyfits.HDUList([hdu])
    B        = hdulist[0].header
    B        = F[0].header.copy()
    hdulist[0].data *= 0.0
    if (plane>=0):
        hdulist[0].data = zeros((N,M), float32)
        hdulist[0].data = 1.0*F[0].data[plane,:,:].reshape(N,M)
    if (data!=None):
        hdulist[0].data = data.copy()
    return hdulist



def MakeEmptyFits(lon, lat, radius, pix, sys_req):
    """
    Make an empty fits object.
    Inputs:
        lon, lat  = centre coordinates of the field [radians]
        radius    = map radius [radians]
        pix       = pixel size [radians]
        sys_req   = coordinate system, WCS_GALACTIC or WCS_J2000
    """
    npix      = int(2.0*radius/pix)+1
    A         = zeros((npix, npix), float32)
    hdu       = pyfits.PrimaryHDU(A)
    F         = pyfits.HDUList([hdu])
    if (0):
        F[0].header.update('CRVAL1',  lon*RADIAN_TO_DEGREE)
        F[0].header.update('CRVAL2',  lat*RADIAN_TO_DEGREE)
        F[0].header.update('CDELT1', -pix*RADIAN_TO_DEGREE)
        F[0].header.update('CDELT2',  pix*RADIAN_TO_DEGREE)
        F[0].header.update('CRPIX1',  0.5*(npix+1))
        F[0].header.update('CRPIX2',  0.5*(npix+1))
        if (sys_req==WCS_GALACTIC):
            F[0].header.update('CTYPE1', 'GLON-TAN')
            F[0].header.update('CTYPE2', 'GLAT-TAN')
            F[0].header.update('COORDSYS', 'GALACTIC')
        else:
            F[0].header.update('CTYPE1', 'RA---TAN')
            F[0].header.update('CTYPE2', 'DEC--TAN')
            F[0].header.update('COORDSYS', 'EQUATORIAL')
            F[0].header.update('EQUINOX', 2000.0)
    else:
        F[0].header.update(CRVAL1 =  lon*RADIAN_TO_DEGREE)
        F[0].header.update(CRVAL2 =  lat*RADIAN_TO_DEGREE)
        F[0].header.update(CDELT1 = -pix*RADIAN_TO_DEGREE)
        F[0].header.update(CDELT2 =  pix*RADIAN_TO_DEGREE)
        F[0].header.update(CRPIX1 =  0.5*(npix+1))
        F[0].header.update(CRPIX2 =  0.5*(npix+1))
        if (sys_req==WCS_GALACTIC):
            F[0].header.update(CTYPE1   = 'GLON-TAN')
            F[0].header.update(CTYPE2   = 'GLAT-TAN')
            F[0].header.update(COORDSYS = 'GALACTIC')
        else:
            F[0].header.update(CTYPE1   = 'RA---TAN')
            F[0].header.update(CTYPE2   = 'DEC--TAN')
            F[0].header.update(COORDSYS = 'EQUATORIAL')
            F[0].header.update(EQUINOX  = 2000.0)
    return F



def MakeEmptyFitsDim(lon, lat, pix, m, n, sys_req):
    """
    Make an empty fits object.
    Inputs:
        lon, lat  = centre coordinates of the field [radians]
        pix       = pixel size [radians]
        m, n      = width and height in pixels
        sys_req   = coordinate system, WCS_GALACTIC or WCS_J2000
    """
    A         = zeros((n, m), float32)
    hdu       = pyfits.PrimaryHDU(A)
    F         = pyfits.HDUList([hdu])
    if (0):
        F[0].header.update('CRVAL1',  lon*RADIAN_TO_DEGREE)
        F[0].header.update('CRVAL2',  lat*RADIAN_TO_DEGREE)
        F[0].header.update('CDELT1', -pix*RADIAN_TO_DEGREE)
        F[0].header.update('CDELT2',  pix*RADIAN_TO_DEGREE)
        F[0].header.update('CRPIX1',  0.5*(m+1))
        F[0].header.update('CRPIX2',  0.5*(n+1))
        if (sys_req==WCS_GALACTIC):
            F[0].header.update('CTYPE1', 'GLON-TAN')
            F[0].header.update('CTYPE2', 'GLAT-TAN')
            F[0].header.update('COORDSYS', 'GALACTIC')
        else:
            F[0].header.update('CTYPE1', 'RA---TAN')
            F[0].header.update('CTYPE2', 'DEC--TAN')
            F[0].header.update('COORDSYS', 'EQUATORIAL')
            F[0].header.update('EQUINOX', 2000.0)
    else:
        F[0].header.update(CRVAL1 =  lon*RADIAN_TO_DEGREE)
        F[0].header.update(CRVAL2 =  lat*RADIAN_TO_DEGREE)
        F[0].header.update(CDELT1 = -pix*RADIAN_TO_DEGREE)
        F[0].header.update(CDELT2 =  pix*RADIAN_TO_DEGREE)
        F[0].header.update(CRPIX1 =  0.5*(m+1))
        F[0].header.update(CRPIX2 = 0.5*(n+1))
        if (sys_req==WCS_GALACTIC):
            F[0].header.update(CTYPE1   = 'GLON-TAN')
            F[0].header.update(CTYPE2   = 'GLAT-TAN')
            F[0].header.update(COORDSYS = 'GALACTIC')
        else:
            F[0].header.update(CTYPE1   = 'RA---TAN')
            F[0].header.update(CTYPE2   = 'DEC--TAN')
            F[0].header.update(COORDSYS = 'EQUATORIAL')
            F[0].header.update(EQUINOX  = 2000.0)
    return F



def ResampleFits(A, HDR, fwhm=0.0, fwhm_ori=None, hdu=0, scanamorphos=-1):
    """
    Resample FITS image A onto the pixels defined by HDR, optionally smoothing the
    data first with a gaussian beam with the given fwhm [radians].
    Input:
        A         =   origin of the data, a pyfits object (no default)
        HDR       =   pyfits header object defining the target pixelisation (no default)
        fwhm      =   if >0, convolve data to this resolution [radians] when on B pixels (0.0)
        fwhm_ori  =   fwhm [radians] of the original data in A  (None)
        hdu       =   the hdu number of A to use (default 0)
        scanamorphos = if >=0, assume A is a scanamorphos map = 3D cube
                       from which plane scanamorphos is extracted (default -1)
    Output:
        new pyfits object with the resampled (and convolved) image.
    Note:
        If fwhm is given, also fwhm_ori must be specified.
        However, as a back up routine looks for keyword MJ_FWHM [degrees].
        This is also written to the header of the returned pyfits object.
    """
    # print fwhm, fwhm_ori, hdu, scanamorphos
    if ((fwhm>0.0)&(fwhm_ori==None)):
        try:
            fwhm_ori = HDR['MJ_FWHM']*DEGREE_TO_RADIAN
        except:
            print('ResampleFits: fwhm_ori must be specified!')
        if (fwhm_ori==None): 
            return None        
    if (scanamorphos>=0):
        print('Input map is scanamorphos ??', A)
        A.verify('silentfix')
        A[hdu].data = A[hdu].data[scanamorphos,:,:].copy()
        A[hdu].header['NAXIS'] = 2                        
    if (len(A)>1):
        print('Extract hdu %d' % hdu)
        A = ExtractHDU(A, hdu=hdu)
    A.verify('silentfix')
    if (fwhm>0.0):
        print('ConvolveFits')
        A = ConvolveFits(A, fwhm, fwhm_orig=fwhm_ori, hdu=0)       
    A.verify('silentfix')
    A.writeto('tmp003.fits', clobber=True, output_verify='silentfix')
    HDR.toTxtFile('tmp.header', clobber=True)
    montage.reproject('tmp003.fits', 'tmp004.fits', header='tmp.header', exact_size=True)
    F = pyfits.open('tmp004.fits')
    F[0].header.update('MJ_FWHM', fwhm*RADIAN_TO_DEGREE)
    return F


    
def ReprojectFits(A, B, hdu=[0,0]):
    """
    Reproject Fits A to pixels of fits B using montage.
    """
    B[hdu[1]].header.toTxtFile('tmp.header', clobber=True)
    A.writeto('tmp000.fits', clobber=True, output_verify='silentfix')
    montage.reproject('tmp000.fits', 'tmp001.fits', header='tmp.header', exact_size=True)
    C = pyfits.open('tmp001.fits')
    return C

    
def GetFitsMask(A, fwhm, hdu=0, mini=-99.0, maxi=1e12, ignore=0.0, bad_limit=0.95):
    """
    Return a 2D  mask by checking valid pixels in the FITS image and extending 
    masked area by convolution with given fwhm.
    Input:
        A          =  input pyfits object
        fwhm       =  FWHM [radians] of the convolving kernel
        bad_limit  =  ok pixels are set to 1, bad 0 => set to zero also pixels that
                      in the convolution of the mask would fall below bad_limit
        hdu        =  number of the header unit (default=0)
        mini, maxi =  limits for valid pixels
                      nonzero((X>mini)&(X<maxi)&(X!=ignore)) == valid pixels
        ignore     =  ignore pixels with this value (e.g., missing value)
    """
    pix    = get_fits_pixel_size(A, hdu) * DEGREE_TO_RADIAN
    FWHM   = fwhm/pix
    X      = A[hdu].data.copy()
    m      = nonzero((X>mini)&(X<maxi)&(X!=ignore))  # ok pixels
    X      = zeros(A[0].data.shape, float32)
    X[m]   = 1.0      # ok pixels
    X      = convolve_map_fast(X.copy(), FWHM)
    bad    = nonzero((X<bad_limit)|(X==0.0))
    X      = ones(A[0].data.shape, float32)
    X[bad] = 0.0
    return X



#########################################################################################
############################  From MJ/Aux/Dust.py #######################################
#########################################################################################



def PlanckFunction(freq, T):
    """
    Return value of the Planck function. 
    Input parameters frequency (freq, [Hz]) and temperature (T, [K]).
    Returns the intensity in cgs units ([erg/s/cm2/sr/Hz]).
    """
    return (2.0*PLANCK*(freq/C_LIGHT)**2.0*freq) / (exp(H_K*freq/T)-1.0)


def MBB_function(f, T, B):
    # Modified Blackbody function
    return  PlanckFunction(f, T)/PlanckFunction(um2f(200.0), T)  * (f/um2f(200.0))**B
                


def make_cc_spline2d(filter, um, Tmin, Tmax, dT, beta_min, beta_max, dbeta):
    """
    Prepare a 2d spline object that can be used for calculation of colour corrections for
    (T, beta) pairs. In case of large maps this should be faster than re-calculating the
    colour correctin for each pixel separately.
    Input:
        filter  =  name of the detector filter
        um      =  nominal wavelength [um] of the filter
        Tmin, Tmax, dT = requested temperature grid
        beta_min, beta_max, dbeta = requested beta grid
    Returns
        tck     =  coefficients for the 2d spline
                   bisplev(T, beta, tck) returns the colour correction coefficient
    """
    # print 'make_cc_spline2d'
    nbeta    = int((beta_max-beta_min)/float(dbeta)+1)
    nT       = int((Tmax-Tmin)/float(dT)+1.0)
    B        = numpy.linspace(beta_min-0.1*dbeta, beta_max+0.1*dbeta, nbeta)
    T        = numpy.linspace(Tmin-0.1*dT, Tmax+0.1*dT, nT)
    B2d, T2d = meshgrid(B, T)  # T = rows    
    N, M     = T2d.shape
    Y        = zeros((N,M), float32)
    for i in range(N):
        for j in range(M):
            t       = T2d[i,j]
            b       = B2d[i,j]            
            Y[i, j] = CalculateColourCorrection(filter, um2f(um), b, t)
            #print 'Y  ', Y[i,j]
    #
    tck = bisplrep(T2d, B2d, Y, s=0.0001, full_output=0, quiet=1)
    # print 'make_cc_spline2d ---- DONE!'
    return tck

        
        
def HenyeyGreenstein(theta, g):
    """
    Return value of the Henyey-Greenstein function for angle
    theta (rad) and asymmetry parameter g)
    """
    return  ((1-g*g)/(4.0*pi))*pow(1.0+g*g-2.0*g*cos(theta), -1.5)


def ModifiedBlackbody(f, T, beta):
    """
    Return value of modified black body curve B(T)*f^beta at given frequency f [Hz].
    Result will be normalized with value at 200um.
    """
    ### f0 = C_LIGHT/200.0e-4
    return ((f/1.4989623e+12)**(3.0E0+beta)) * (exp(H_K*1.4989623e+12/T)-1.0) / (exp(H_K*f/T)-1.0)


def ModifiedBlackbody250(f, T, beta):
    """
    Return value of modified black body curve B(T)*f^beta at given frequency f [Hz].
    Result will be normalized with value at 200um.
    """
    return ((f/1.19917e+12)**(3.0E0+beta)) * (exp(H_K*1.19917e+12/T)-1.0) / (exp(H_K*f/T)-1.0)


def ModifiedBlackbody2B(f, fx, t, beta1, beta2):
    if (isscalar(f)):
        if (f>=f0):
            return ModifiedBlackbody(f, t, beta1)
        else:
            return ModifiedBlackbody(f, t, beta2) * (fx/1.4989623e+12)**(beta1-beta2)
    else:
        res = zeros(len(f), float32)
        m = nonzero(f>=fx)
        res[m] = ModifiedBlackbody(f[m], t, beta1)
        m = nonzero(f<fx)
        res[m] = ModifiedBlackbody(f[m], t, beta2) * (fx/1.4989623e+12)**(beta1-beta2)
    return res


def BBChi2(p, f, S, dS):
    """
    Given a modified blackbody curve with parameters p = [ I200, T, beta ],
    observed frequencies f and intensities S+-dS, return chi2 error.
    """
    #### I  = p[0]*ModifiedBlackbody(f, p[1], p[2])   # model predictions
    I  = p[0]*((f/1.4989623e+12)**(3.0E0+p[2])) * (exp(H_K*1.4989623e+12/p[1])-1.0) / (exp(H_K*f/p[1])-1.0)
    I  = (I-S)/dS
    #print f2um(f)
    #print p[0], p[1], p[2]
    #print sum(I*I)
    # add penalty for odd beta ??
    #### return dot(I,I)/len(I)
    return sum(I*I)/len(I)



def BBChi2_250(p, f, S, dS):
    """
    Given a modified blackbody curve with parameters p = [ I250, T, beta ],
    observed frequencies f and intensities S+-dS, return chi2 error.
    """
    I  = p[0]*((f/1.19917e+12)**(3.0E0+p[2])) * (exp(H_K*1.19917e+12/p[1])-1.0) / (exp(H_K*f/p[1])-1.0)
    I  = (I-S)/dS
    return sum(I*I)/len(I)



def BBChi2Cov(p, f, S, dS, COVI):
    """
    Given a modified blackbody curve with parameters p = [ I200, T, beta ],
    observed frequencies f and intensities S+-dS, return chi2 error.
    COVI = inverse of the covariance matrix for values between bands
    """
    I  = p[0]*((f/1.4989623e+12)**(3.0E0+p[2])) * (exp(H_K*1.4989623e+12/p[1])-1.0) / (exp(H_K*f/p[1])-1.0)
    #  instead of    [(x-u)/dx]^2
    #     I   = (I-S)/dS
    #     res = sum(I*I)/len(I)
    #  we have       (x-u)*Inv(cov)*(x-u)
    d   =  I-S
    res =  dot(d, dot(COVI, d))
    return sum(res)


def BBChi2Cov_250(p, f, S, dS, COVI):
    """
    Given a modified blackbody curve with parameters p = [ I250, T, beta ],
    observed frequencies f and intensities S+-dS, return chi2 error.
    COVI = inverse of the covariance matrix for values between bands
    """
    I  = p[0]*((f/1.19917e+12)**(3.0E0+p[2])) * (exp(H_K*1.19917e+12/p[1])-1.0) / (exp(H_K*f/p[1])-1.0)
    #  instead of    [(x-u)/dx]^2
    #     I   = (I-S)/dS
    #     res = sum(I*I)/len(I)
    #  we have       (x-u)*Inv(cov)*(x-u)
    d   =  I-S
    res =  dot(d, dot(COVI, d))
    return sum(res)


def BBChi2_2B(p, f1, f2, S1, dS1, S2, dS2, fx):
    """
    Given a modified blackbody curve with parameters p = [ I200, T, beta1, beta2 ],
    observed frequencies f and intensities S+-dS, return chi2 error.
       f1, S1, dS1  fitted with    p[0]*MBB(p[1], p[2])
       f2, S2, dS2  fitted with    p[0]*MBB(p[1], p[3]) *  (fx/f0)**(beta1-beta2)
    """
    # f0 = um2f(200.0) = 1498962290000.0
    I  = p[0]*((f1/1.4989623e+12)**(3.0E0+p[2])) * (exp(H_K*1.4989623e+12/p[1])-1.0) / (exp(H_K*f1/p[1])-1.0)
    I  = (I-S1)/dS1
    chi2 = sum(I*I)/len(I)
    # the MBB with beta2
    I  = p[0]*((f2/1.4989623e+12)**(3.0E0+p[3])) * (exp(H_K*1.4989623e+12/p[1])-1.0) / (exp(H_K*f2/p[1])-1.0)
    I *= (fx/1.4989623e+12)**(p[2]-p[3])
    I  = (I-S2)/dS2
    chi2 += sum(I*I)/len(I)    
    return chi2


def BBChi2_fix_T(p, f, S, dS, T):
    """
    Given a modified blackbody curve with parameters p = [ I200, beta ],
    observed frequencies f and intensities S+-dS, return chi2 error.
    p = [ I, beta], temperature T is kept constant.
    """
    I  = p[0]*ModifiedBlackbody(f, T, p[1])   # model predictions
    I  = (I-S)/dS
    return sum(I*I)/len(I)

    
def BBChi2_fix_T_250(p, f, S, dS, T):
    """
    Given a modified blackbody curve with parameters p = [ I250, beta ],
    observed frequencies f and intensities S+-dS, return chi2 error.
    p = [ I, beta], temperature T is kept constant.
    """
    I  = p[0]*ModifiedBlackbody250(f, T, p[1])   # model predictions
    I  = (I-S)/dS
    return sum(I*I)/len(I)

    
def BBChi2_fix_beta(p, f, S, dS, beta):
    """
    Given a modified blackbody curve with parameters p = [ I200, T ],
    observed frequencies f and intensities S+-dS, return chi2 error.
    p = [ I, T], emissivity index beta is kept constant.
    """
    # print 'p f S dS beta', p, f, S, dS, beta
    I  = p[0]*ModifiedBlackbody(f, p[1], beta)   # model predictions
    I  = (I-S)/dS
    return sum(I*I)/len(I)


def BBChi2_fix_beta_250(p, f, S, dS, beta):
    """
    Given a modified blackbody curve with parameters p = [ I250, T ],
    observed frequencies f and intensities S+-dS, return chi2 error.
    p = [ I, T], emissivity index beta is kept constant.
    """
    I  = p[0]*ModifiedBlackbody250(f, p[1], beta)   # model predictions
    I  = (I-S)/dS
    return sum(I*I)/len(I)



def Delta_fix_beta(p, f, S, dS, beta):
    """
    Given a modified blackbody curve with parameters p = [ I200, T ],
    observed frequencies f and intensities S+-dS, return normalised residuals.
    p = [ I, T], emissivity index beta is kept constant.
    """
    # print 'p f S dS beta', p, f, S, dS, beta
    I  = p[0]*ModifiedBlackbody(f, p[1], beta)   # model predictions
    I  = (I-S)/dS
    return I/len(I)


def Deriv_fix_beta(p, f, S, dS, beta):
    """
    Return derivatives - NOT TESTED
    """
    D  = zeros((len(S), 2), float32)
    for i in range(len(S)):
        D[i,0] = ModifiedBlackbody(f, p[1], beta)/dS[i]
        D[i,1] = D[i,0] * H_K/(p[1]*p[1]) / (1.0-exp(-H_K*f[i]/p[1]))
    return D



def Delta_fix_T(p, f, S, dS, fixed_T):
    """
    Given a modified blackbody curve with parameters p = [ I200, beta ],
    observed frequencies f and intensities S+-dS, return chi2 error.
    p = [ I, T], emissivity index beta is kept constant.
    """
    I  = p[0]*ModifiedBlackbody(f, fixed_T, p[1])
    I  = (I-S)/dS
    return I/len(I)


def Deriv_fix_T(p, f, S, dS, fixed_T):
    """
    Return derivatives - NOT TESTED
    """
    f0 = um2f(200.0)
    D  = zeros((len(S), 2), float32)
    for i in range(len(S)):
        D[i,0] = ModifiedBlackbody(f[i], fixed_T, p[1])/dS[i]   #  d/dI
        D[i,1] = D[i,0] *  p[1] * (f0/f[i])                     #  d/dbeta
    return D



def FitModifiedBlackbody(um, S, dS, I200, T, beta, fix_T=False, fix_beta=False, MC=1,
                        Tmin=3.0, Tmax=35.0, beta_min=0.5, beta_max=4.0, 
                        powell=False, xtol=1.0e-5, ftol=1.0e-9, return_samples=False,
                        COV=[]):
    """
    Usage:
        p     = FitModifiedBlackbody(um, S, dS, I200, T, beta, fix_T=False, fix_beta=False, MC=1)
        p, dP, [samples] = FitModifiedBlackbody(um, S, dS, I200, T, beta, fix_T=False, fix_beta=False, MC=100)
    Given vector of wavelengths and intensities with uncertainties,
    fit a modified black body curve to the data. Initial values should be
    given for intensity at 200um, temperature, and emissivity index.
    If fix_T or fix_beta is true, the corresponding parameter is kept constant.
    Input:
        um       = vector of wavelenths [um]
        S        = intensity values
        dS       = uncertainties of the intensities
        I200     = initial value for 200um intensity (input)
        T        = initial value for temperature [K]
        beta     = initial value for spectral index
        fix_T    = if True, do fit keeping temperature fixed
        fix_beta = if True, do fit keeping spectral index fixed
        MC       = if >1, do this many Monte Carlo samples and return also the
                   estimated errors of the parameters (returns quartile interval scaled to std)
        Tmin, Tmax         = temperature values cut to this interval
        beta_min, beta_max = spectral indices cut to this interval
        return_samples     = if True, return array of samples MC x [1-3]
        COV      = optional full covariance matrix between S in different bands
    Returns:
        p, dp    = vector of fit results; if MC>1, dp contains the estimated errors
                   p depends on the input parameters
                      if fix_T     =>  [I, beta]
                      if fix_beta  =>  [I, T]
                      else         =>  [I, T, beta]
    """
    f   = C_LIGHT/(1.0e-4*asarray(um))
    NF  = len(f)
    tmp = None 
    if (MC>1):
        if ((fix_T)|(fix_beta)):
            tmp = zeros((MC, 2), float32)
        else:
            tmp = zeros((MC, 3), float32)
    else:
        MC = 1
    ###
    COVI, CHOL = [], []
    if (len(COV)>0):
        COVI = inv(COV)
        CHOL = cholesky(COV)
    ###
    for iii in range(MC):
        if (iii==(MC-1)):  # no Monte Carlo = the last one
            SS = S
        else:
            if (len(COV)<1):
                SS = S + dS*numpy.random.standard_normal(NF)
            else:
                # generate photometric errors using the covariance matrix
                SS = S +  dot(CHOL, randn(NF))
                #########################################################
        if (fix_T):        #  T fixed
            p0  = [ I200, max(beta_min, min(beta_max, beta)) ]
            if (1):
                p1  = fmin(BBChi2_fix_T, p0, args=(f, SS, dS, T), disp=0)
            else:
                # *** NOT OK ***
                # p1, p1Con = leastsq(Delta_fix_T, p0, args=(f, SS, dS, beta), Dfun=Deriv_fix_T)
                p1, p1Con = leastsq(Delta_fix_T, p0, args=(f, SS, dS, beta), Dfun=None)
        elif (fix_beta):   # BETA fixed
            p0  = [ I200, max(Tmin, min(Tmax, T)) ]
            if (1):
                p1  = fmin(BBChi2_fix_beta, p0, args=(f, SS, dS, beta), disp=0)
            else:          # NOT TESTED !!!!!!!!!
                p1, p1Con = leastsq(Delta_fix_beta, p0, args=(f, SS, dS, beta), Dfun=Deriv_fix_beta, disp=0)
        else:              # T and BETA both free
            p0  = [ I200, max(Tmin, min(Tmax, T)), max(beta_min, min(beta_max, beta)) ]
            if (powell):
                p1  = fmin_powell(BBChi2, p0, args=(f, SS, dS), disp=0, xtol=xtol, ftol=ftol)
            else:
                if (len(COV)>0):
                    p1   = fmin(BBChi2Cov, p0, args=(f, SS, dS, COVI), disp=0, xtol=xtol, ftol=ftol)
                else:
                    p1  =  fmin(BBChi2,    p0, args=(f, SS, dS),       disp=0, xtol=xtol, ftol=ftol)
        if (MC>1):
            tmp[iii,:] = p1
    if (MC<=1):
        return p1
    else:
        dp1 = []
        if (fix_T):
            m = nonzero((tmp[:,1]>beta_min)&(tmp[:,1]<beta_max))
        elif (fix_beta):
            m = nonzero((tmp[:,1]>Tmin)&(tmp[:,1]<Tmax))
        else:
            m = nonzero((tmp[:,1]>Tmin)&(tmp[:,1]<Tmax)&(tmp[:,2]>beta_min)&(tmp[:,2]<beta_max))
        for i in range(tmp.shape[1]):
            print('PRCTILE OF DATA ... ', len(tmp[:,i][m]))
            if (len(tmp[:,i][m])>3):                
                a, b = matplotlib.mlab.prctile(tmp[:,i][m], (25.0, 75.0))
            else:
                a, b = 0.0, 0.0
            dp1.append( (b-a)*0.74130 )  #  IQ -> sigma
        if (return_samples):
            return p1, asarray(dp1), tmp
        else:
            return p1, asarray(dp1) # p1 = always on the original data



        
def MBB_fit_chi2(UM, S, dS=None, fixed_T=None, fixed_beta=None, filters=None, 
            Tmin=5.0, Tmax=29.0, beta_min=0.5, beta_max=5.0,
            T_bins=200, beta_bins=80,
            result_I=None, result_T=None, result_B=None, result_ind0=None,
            spline_smooth=2.0e-7):
    """
    Fit a modified black body curve to the given observations.
    Input:
        UM            =  wavelengths [um]
        S, dS         =  intensity and its uncertainty 
                         S = [  vector(um1), vector(um2), ... ] => S[i] are observations of one pixel
        fixed_T       =  if given, the temperature will be fixed to this value
        fixed_beta    =  if given, the spectral index will be fixed to this value
        filters       =  if given, should contain an array with filter names (or None) 
                         for each of the bands so that the colour corrections can be
                         taken into account (S should in this case be the in-band values)
        Tmin, Tmax          =  limits for the allowed temperature
        beta_min, beta_max  =  limits for the allowed spectral index
        T_bins, beta_bins   =  size of internal tables (defaults 200 and 50)
        result_I, result_T, result_B, result_ind0
                            = if given, copy results to these arrays starting with position result_ind0
    Return:
        I, T, B        =   arrays of 200um surface brightness and estimated temperature and spectral index
                           three arrays returned even with fixed_T, fixed_beta
                           these are values of the MBB that, after de-colour correction, fits the
                           observed in-band values S
    Notes:
        even a small amount of smoothing introduces noticeable errors at low beta
        it may be better to keep beta_min below the expected range, still not small
        enough to cause oscillations
    """
    # print 'MBB_fit_chi2()'
    from MJ.mjDefs import um2f, PlanckFunction, PLANCK, HOMEDIR
    from numpy import linspace, float64, zeros, float32, meshgrid, ravel, ones
    from scipy.interpolate import bisplrep, bisplev
    from scipy.optimize import fmin
    # pre-calculate SEDs for a range of temperatures and spectral indices
    NT, NB, NW = T_bins, beta_bins, len(UM)
    if (  beta_min<0.1  ):
        T   =  linspace(Tmin, Tmax, NT)
        B   =  linspace(beta_min, beta_max, NB)
    else:
        T   =  numpy.logspace(log10(Tmin), log10(Tmax), NT)
        B   =  numpy.logspace(log10(beta_min), log10(beta_max), NB)        
    
    # make 2d spline interpolation functions for each wavelength
    print('Make 2d spline interpolation...')
    TCK = []        
    if (fixed_T!=None):
        Y = zeros(len(B), float32)
        for iw in range(NW):  # loop over UM
            for i in range(len(B)): # loop over beta bins
                um   =  UM[iw]
                y    =  PlanckFunction(um2f(um), fixed_T) / PlanckFunction(um2f(200.0), fixed_T)
                y   *=  (um2f(um)/um2f(200.0))**B[i]
                if (filters!=None):
                    if (filters[iw]!=None):
                        y *= CalculateColourCorrection(filters[iw], um2f(um), B[i], fixed_T)
                Y[i] = y   #      Y[T, beta, wavelength]
            tck = splrep(B, Y, s=spline_smooth, full_output=0, quiet=0)
            TCK.append(tck)
    elif (fixed_beta!=None):
        Y = zeros(len(T), float32)
        for iw in range(NW):            
            for i in range(len(T)):
                um   =  UM[iw]
                y    =  PlanckFunction(um2f(um), T[i]) / PlanckFunction(um2f(200.0), T[i])
                y   *=  (um2f(um)/um2f(200.0))**fixed_beta
                if (filters!=None):
                    if (filters[iw]!=None):
                        y *= CalculateColourCorrection(filters[iw], um2f(um), fixed_beta, T[i])
                Y[i] = y   #      Y[T, beta, wavelength]
            tck = splrep(T, Y, s=spline_smooth, full_output=0, quiet=0)
            TCK.append(tck)
            if (0):
                clf()
                plot(T, Y, 'x')
                plot(T, splev(T, tck), 'r.')
                draw()
                SHOW()
    else: # T and beta both free parameters
        Y    =  zeros((NT, NB, NW), float32)
        B2d, T2d = meshgrid(B, T)  # T = rows
        N, M = T2d.shape
        for i in range(N):      # over the rows = different T
            for j in range(M):  # over columns = different beta
                t = T2d[i,j]
                b = B2d[i,j]
                for iw in range(NW):
                    um = UM[iw]
                    # Y = intensity normalized so that value at 200um is 1.0
                    y   =  PlanckFunction(um2f(um), t) / PlanckFunction(um2f(200.0), t)
                    y  *=  (um2f(um)/um2f(200.0))**b
                    # if filter is known, de-colour correct this value so that it can be
                    # compared with un-colour corrected values
                    if (filters!=None):
                        if (filters[iw]!=None):
                            y *= CalculateColourCorrection(filters[iw], um2f(um), b, t)
                    Y[i, j, iw] = y   #      Y[T, beta, wavelength]
        # Note!!  Y are IN-BAND values corresponding to MONOCHROMATIC values of
        #         S200, beta and temperature
        # when in-band observations are fitted, results are (S,T,beta) for monochromatic data
        for iw in range(NW):
            tck = bisplrep(T2d, B2d, Y[:,:,iw].copy(), s=spline_smooth, full_output=0, quiet=0)
            # tck = bisplrep(T2d, B2d, Y[:,:,iw].copy(), s=0.0001, full_output=0, quiet=0, kx=1, ky=1)
            print('Add tck to TCK')
            TCK.append(tck)
            print(' ... added done')
            
        if (0):
            iw = 0
            clf()
            subplot(221)
            imshow(log10(Y[:,:,iw]), interpolation='nearest', aspect='auto')
            colorbar()
            subplot(222)
            yyy = 0.0*T2d
            for i in range(T2d.shape[0]):
                for j in range(T2d.shape[1]):
                    yyy[i,j] = bisplev(T2d[i,j], B2d[i,j], TCK[iw])                    
            imshow((yyy-Y[:,:,iw])/Y[:,:,iw], interpolation='nearest', aspect='auto')
            colorbar()
            subplot(223)
            plot(Y[10,:,iw], yyy[10,:], 'k.')
            xlabel('Input')
            ylabel('Spline')
            subplot(224)
            plot(Y[:,10,iw], yyy[:,10], 'k.')
            xlabel('Input')
            ylabel('Spline')
            SHOW()        

            
    print('Make 2d spline interpolation... DONE')
    
    if (S==None): return TCK
        
    # now (T, beta) can be determined by finding Y[T,beta,:] that gives the smallest chi2
    shape = S[0].shape
    II = zeros(len(ravel(S[0])), float32)
    TT = zeros(len(ravel(S[0])), float32)
    BB = zeros(len(ravel(S[0])), float32)
    s, ds = zeros(NW, float32), ones(NW, float32)
    # the initial estimate for S(200um) will be based on S[1], calculate S(200)/S(UM[1])
    ratio_S200_S1 = (PlanckFunction(um2f(200.0), 16.0)/PlanckFunction(um2f(UM[1]), 16.0)) * (um2f(200.0)/um2f(UM[1]))**2.0
    print('Start loop over pixels...')
    for ipix in range(len(ravel(S[0]))):  # loop over pixels
        # extract a vector of intensities in this pixel
        for iw in range(NW):
            s[iw] = ravel(S[iw])[ipix]
            if (dS[iw]!=None):
                ds[iw] = ravel(dS[iw])[ipix]
            else:
                ds[iw] = 0.2*s[iw]  # ???
        # find the scaling and the (T,beta) that minimize chi2
        p0 = asarray([ clip(s[1]*ratio_S200_S1, 1.0e-3, 1.0e12), 16.0, 2.0 ], float32)
        #
        p1 = fmin(MBB_fit_chi2_fun, p0, 
             args=(TCK, s, ds, Tmin, Tmax, beta_min, beta_max, fixed_T, fixed_beta),
             xtol=1.0e-6, ftol=1.0e-6, disp=0)
        #
        # @@@  MBB_fit_chi2
        # delta = zeros(len(s), float32)
        # p1, p1Con, info, msg, ierr = leastsq(Delta_fun, p0, 
        #     args=(delta, TCK, s, ds, Tmin, Tmax, beta_min, beta_max, fixed_T, fixed_beta), 
        #     Dfun=Deriv_fun,
        #     xtol=1.0e-4, ftol=1.0e-4, full_output=True)
        # if ((ierr<1)|(ierr>2)):
        #    print 'IERR ', ierr, '%5.2f' % p1[1], msg
        #
        
        # print p1
        II[ipix] = p1[0]
        if (fixed_T!=None):        # fixed T
            TT[ipix] = fixed_T
            BB[ipix] = p1[1]
        elif (fixed_beta!=None):   # fixed beta
            TT[ipix] = p1[1]
            BB[ipix] = fixed_beta
        else:                      # (T,beta) both free
            TT[ipix] = p1[1]
            BB[ipix] = p1[2]
    ###
    if (result_ind0!=None):
        result_I[result_ind0:(result_ind0+len(TT))] = II
        result_T[result_ind0:(result_ind0+len(TT))] = TT        
        result_B[result_ind0:(result_ind0+len(TT))] = BB
    else:
        II.shape = shape
        TT.shape = shape
        BB.shape = shape
        return II, TT, BB



    

def MBB_fit_LSQ(UM, S, dS=[], fixed_T=None, fixed_beta=None, filters=None, cc_routine=None,
            Tmin=5.0, Tmax=29.0, beta_min=0.5, beta_max=5.0,
            T_bins=200, beta_bins=80,
            result_I=None, result_T=None, result_B=None, result_ind0=None, TCK=[],
            spline_smooth=1e-6):
    """
    Fit a modified black body curve to the given observations.
    Input:
        UM            =  wavelengths [um]
        S, dS         =  intensity and its uncertainty 
                         S = [  vector(um1), vector(um2), ... ] => S[i] are observations of one pixel
        fixed_T       =  if given, the temperature will be fixed to this value
        fixed_beta    =  if given, the spectral index will be fixed to this value
        filters       =  if given, should contain an array with filter names (or None) 
                         for each of the bands so that the colour corrections can be
                         taken into account (S should in this case be the in-band values)
        cc_routine    =  if given,  cc_routine(filter, T, beta) returns the colour correction
                         overrides the normal behaviour where corrections calculated based on filter profiles
        Tmin, Tmax          =  limits for the allowed temperature
        beta_min, beta_max  =  limits for the allowed spectral index
        T_bins, beta_bins   =  size of internal tables (defaults 200 and 50)
        result_I, result_T, result_B, result_ind0
                            = if given, copy results to these arrays starting with position result_ind0
    Return:
        I, T, B        =   arrays of 200um surface brightness and estimated temperature and spectral index
                           three arrays returned even with fixed_T, fixed_beta
                           these are values of the MBB that, after de-colour correction, fits the
                           observed in-band values S
    Notes:
        Even a small amount of smoothing introduces noticeable errors at low beta.
        It may be better to keep beta_min below the expected range, still not small
        enough to cause oscillations.
        2017-03-21  T=5.0-35K, beta=0.5-5.0 is about ok,
                    extend beta to 0.0 and you are in trouble (even without colour corrections!)
    """
    print('MBB_fit_LSQ() .... spline_smooth %.3e' % spline_smooth)
    from Aux import um2f, PlanckFunction, PLANCK
    #from MJ.mjDefs import um2f, PlanckFunction, PLANCK, HOMEDIR
    from numpy import linspace, float64, zeros, float32, meshgrid, ravel, ones
    from scipy.interpolate import bisplrep, bisplev
    from scipy.optimize import fmin, leastsq
    
    # special case -- if fixed_beta used and that is an array, we need still
    # TCK as function of both T and beta; beta value set inside the loop from fixed_array[]
    FIXED_BETA_ARRAY = False
    if (fixed_beta!=None):
        if (isscalar(fixed_beta)):
            FIXED_BETA_ARRAY = False
        else:
            FIXED_BETA_ARRAY = True
            
            
    # pre-calculate SEDs for a range of temperatures and spectral indices
    NT, NB, NW = T_bins, beta_bins, len(UM)
    if (  beta_min<0.1  ):
        T   =  linspace(Tmin, Tmax, NT)
        B   =  linspace(beta_min, beta_max, NB)
    else:
        T   =  numpy.logspace(log10(Tmin), log10(Tmax), NT)
        B   =  numpy.logspace(log10(beta_min), log10(beta_max), NB)        
    
    # make 2d spline interpolation functions for each wavelength
    if (TCK==[]):
        print('Make 2d spline interpolation...')
        TCK = []        
        if (fixed_T!=None):
            Y = zeros(len(B), float32)
            for iw in range(NW):  # loop over UM
                for i in range(len(B)): # loop over beta bins
                    um   =  UM[iw]
                    y    =  PlanckFunction(um2f(um), fixed_T) / PlanckFunction(um2f(200.0), fixed_T)
                    y   *=  (um2f(um)/um2f(200.0))**B[i]
                    if (filters!=None):
                        if (filters[iw]!=None):
                            if (cc_routine==None):
                                y *= CalculateColourCorrection(filters[iw], um2f(um), B[i], fixed_T)
                            else:
                                y *= cc_routine(filters[iw], B[i], fixed_T)
                    Y[i] = y   #      Y[T, beta, wavelength]
                tck = splrep(B, Y, s=spline_smooth, full_output=0, quiet=0)
                TCK.append(tck) #  fixed T, free beta
            print('****** TCK FIXED T *******')
        elif ((fixed_beta!=None)&(not(FIXED_BETA_ARRAY))): # no beta dependence, only one beta
            Y = zeros(len(T), float32)
            for iw in range(NW):            
                for i in range(len(T)):
                    um   =  UM[iw]
                    y    =  PlanckFunction(um2f(um), T[i]) / PlanckFunction(um2f(200.0), T[i])
                    y   *=  (um2f(um)/um2f(200.0))**fixed_beta
                    if (filters!=None):
                        if (filters[iw]!=None):
                            if (cc_routine==None):
                                y *= CalculateColourCorrection(filters[iw], um2f(um), fixed_beta, T[i])
                            else:
                                y *= cc_routine(filters[iw], fixed_beta, T[i])
                    Y[i] = y   #      Y[T, beta, wavelength]
                tck = splrep(T, Y, s=spline_smooth, full_output=0, quiet=0)
                TCK.append(tck)  # fixed beta, T free
                if (0):
                    clf()
                    plot(T, Y, 'x')
                    plot(T, splev(T, tck), 'r.')
                    draw()
                    SHOW()
                    sys.exit()
            print('****** TCK FIXED BETA *******')
        else: # T and beta both free parameters -- or FIXED_BETA_ARRAY !!
            Y    =  zeros((NT, NB, NW), float32)
            B2d, T2d = meshgrid(B, T)  # T = rows
            N, M = T2d.shape
            for i in range(N):      # over the rows = different T
                for j in range(M):  # over columns = different beta
                    t = T2d[i,j]
                    b = B2d[i,j]
                    for iw in range(NW):
                        um = UM[iw]
                        # Y = intensity normalized so that value at 200um is 1.0
                        y   =  PlanckFunction(um2f(um), t) / PlanckFunction(um2f(200.0), t)
                        y  *=  (um2f(um)/um2f(200.0))**b
                        # if filter is known, de-colour correct this value so that it can be
                        # compared with un-colour corrected values
                        if (filters!=None):
                            if (filters[iw]!=None):
                                if (cc_routine==None):
                                    y *= CalculateColourCorrection(filters[iw], um2f(um), b, t)
                                else:
                                    y *= cc_routine(filters[iw], b, t)
                        Y[i, j, iw] = y   #      Y[T, beta, wavelength]
            # Note!!  Y are IN-BAND values corresponding to MONOCHROMATIC values of
            #         S200, beta and temperature
            # when in-band observations are fitted, results are (S,T,beta) for monochromatic data
            for iw in range(NW):
                tck = bisplrep(T2d, B2d, Y[:,:,iw].copy(), s=spline_smooth, full_output=0, quiet=0)
                # tck = bisplrep(T2d, B2d, Y[:,:,iw].copy(), s=0.0001, full_output=0, quiet=0, kx=1, ky=1)
                print('Add tck to TCK')
                TCK.append(tck)   # (T,beta) free or T free and we have array for fixed_beta
                print(' ... added done')
            print('****** TCK FREE BETA AND T ---- OR FIXED_BETA_ARRAY *******')
                            
            if (0):
                iw = 0
                clf()
                subplot(221)
                imshow(log10(Y[:,:,iw]), interpolation='nearest', aspect='auto')
                colorbar()
                subplot(222)
                yyy = 0.0*T2d
                for i in range(T2d.shape[0]):
                    for j in range(T2d.shape[1]):
                        yyy[i,j] = bisplev(T2d[i,j], B2d[i,j], TCK[iw])                    
                imshow((yyy-Y[:,:,iw])/Y[:,:,iw], interpolation='nearest', aspect='auto')
                colorbar()
                subplot(223)
                plot(Y[10,:,iw], yyy[10,:], 'k.')
                xlabel('Input')
                ylabel('Spline')
                subplot(224)
                plot(Y[:,10,iw], yyy[:,10], 'k.')
                xlabel('Input')
                ylabel('Spline')
                SHOW()        
    
                
        print('Make 2d spline interpolation... DONE')
        
    if (S==None): return TCK
        
    # now (T, beta) can be determined by finding Y[T,beta,:] that gives the smallest chi2
    shape = S[0].shape
    II = zeros(len(ravel(S[0])), float32)
    TT = zeros(len(ravel(S[0])), float32)
    BB = zeros(len(ravel(S[0])), float32)
    s, ds = zeros(NW, float32), ones(NW, float32)
    # the initial estimate for S(200um) will be based on S[1], calculate S(200)/S(UM[1])
    ratio_S200_S1 = (PlanckFunction(um2f(200.0), 16.0)/PlanckFunction(um2f(UM[1]), 16.0)) * (um2f(200.0)/um2f(UM[1]))**2.0
    print('Start loop over pixels...')
    delta = zeros(NW, float32)
    FIXED_BETA = fixed_beta


    for ipix in range(len(ravel(S[0]))):  # loop over pixels
        # extract a vector of intensities in this pixel
        for iw in range(NW):
            s[iw] = ravel(S[iw])[ipix]
            if (dS!=[]):
                if (dS[iw]!=[]):
                    ds[iw] = ravel(dS[iw])[ipix]
                else:
                    ds[iw] = 0.1*s[iw]  # ???
            else:
                ds[iw] = 0.1*s[iw]  # ???
                
        # find the scaling and the (T,beta) that minimize chi2
        if (FIXED_BETA_ARRAY):
            p0 = [ clip(s[1]*ratio_S200_S1, 1.0e-3, 1.0e12), 16.0 ]  # beta not free parameter
            FIXED_BETA = fixed_beta[ipix] # different for each pixel
        elif (fixed_beta!=None):
            p0 = asarray([ clip(s[1]*ratio_S200_S1, 1.0e-3, 1.0e12), 16.0 ], float32)
        elif (fixed_T!=None):
            p0 = [ clip(s[1]*ratio_S200_S1, 1.0e-3, 1.0e12), 2.0 ]
        else:
            p0 = [ clip(s[1]*ratio_S200_S1, 1.0e-3, 1.0e12), 16.0, 2.0 ]
        
            
            
        # print('FIXED_BETA_ARRAY', FIXED_BETA_ARRAY)
            
        #  @@@  MBB_fit_chi2_LSQ  --- FIXED_BETA is scalar, also for FIXED_BETA_ARRAY
        p1, pC1 = leastsq(Delta_fun, p0, 
                          args=(delta, TCK, s, ds, Tmin, Tmax, beta_min, beta_max, 
                          fixed_T, FIXED_BETA, FIXED_BETA_ARRAY),
                          Dfun=Deriv_fun,
                          xtol=1.0e-4, ftol=1.0e-4)
        
        if ((fixed_T==None)&(fixed_beta==None)):
            if (p1[1]==16.0): # failed ???   retry ???
                p0 = [ clip(s[1]*ratio_S200_S1, 1.0e-3, 1.0e12), 13.0, 2.5 ]                    
                p1, pC1 = leastsq(Delta_fun, p0, 
                                  args=(delta, TCK, s, ds, Tmin, Tmax, beta_min, beta_max, 
                                  fixed_T, FIXED_BETA, FIXED_BETA_ARRAY),
                                  Dfun=Deriv_fun,
                                  xtol=1.0e-4, ftol=1.0e-4)
                if ((fixed_T==None)&(p1[1]==13.0)):  # still failed
                    p1[1] = 15.0  # Warning -- this can hide bad fits in (T, beta) plots!
                    p1[2] =  2.0
                  
        II[ipix] = p1[0]
        if (fixed_T!=None):        # fixed T
            TT[ipix] = fixed_T
            BB[ipix] = p1[1]
        elif (fixed_beta!=None):   # fixed beta
            TT[ipix] = p1[1]
            BB[ipix] = FIXED_BETA  # covers scalar and array fixed_beta cases
        else:                      # (T,beta) both free
            TT[ipix] = p1[1]
            BB[ipix] = p1[2]
    ###
    if (result_ind0!=None):
        result_I[result_ind0:(result_ind0+len(TT))] = II
        result_T[result_ind0:(result_ind0+len(TT))] = TT        
        result_B[result_ind0:(result_ind0+len(TT))] = BB
    else:
        II.shape = shape
        TT.shape = shape
        BB.shape = shape
        return II, TT, BB


        
def Delta_fun(p, delta, TCK, s, ds, Tmin, Tmax, beta_min, beta_max, fixed_T=None, fixed_beta=None,
              FIXED_BETA_ARRAY=False):
    """
    Return fit errors for comparison of theoretical modified black body curve and observations.
    Input:
        delta =  work space, array[number of wavelengths]
        p     =  [ intensity, T, beta ], parameters of the modified black body curve
                 intensity is the value at 200um
        TCK   =  TCK[i] is a 2d spline representation of the i:th wavelength intensity
                 as function of temperature and spectral index (relative to 200um value)
        s, ds =  observed intensity values and their uncertainties
                 s[i] is the i:th wavelength
        fixed_T     =>  p = [ intensity, beta ]
        fixed_beta  =>  p = [ intensity, T    ]
        FIXED_BETA_ARRAY =  if true, use 2d interpolator even if fixed_beta set
    Return:
        vector of normalized deviations from the model
    """
    # print('--------------', FIXED_BETA_ARRAY)
    from   scipy.interpolate import bisplev, splev
    # delta    = zeros(len(s), float32)   # wavelengths
    penalty  = 1.0
    if (fixed_T!=None):  # fixed temperature
        for i in range(len(s)):
            delta[i]   =  (s[i] - p[0]*splev(p[1],TCK[i]))/ds[i]   # p[1] = beta
        if (p[1]<beta_min):
            penalty    *= (1.0+(beta_min-p[1]))
        if (p[1]>beta_max):
            penalty    *= (1.0+3.0*(p[1]-beta_max))
    elif ((fixed_beta!=None)&(not(FIXED_BETA_ARRAY))): # fixed spectral index
        for i in range(len(s)):
            delta[i]    = (s[i] - p[0]*splev(p[1], TCK[i]))/ds[i]  # p[1] = T
        if (p[1]<Tmin):
            penalty    *= (1.0+(Tmin-p[1]))
        if (p[1]>Tmax):
            penalty    *= (1.0+0.5*(p[1]-Tmax))
    elif (FIXED_BETA_ARRAY):
        for i in range(len(s)):
            delta[i] = (s[i] - p[0]*bisplev(clip(p[1],Tmin,Tmax), clip(fixed_beta,beta_min, beta_max), TCK[i]))/ds[i]
        # penalties
        if (p[1]<Tmin):
            penalty *= (1.0+3.0*(Tmin-p[1]))
        if (p[1]>Tmax):
            penalty *= (1.0+3.0*(p[1]-Tmax))        
    else:  # T and beta free --- or using FIXED_BETA_ARRAY
        for i in range(len(s)):
            delta[i] = (s[i] - p[0]*bisplev(clip(p[1],Tmin,Tmax), clip(p[2],beta_min, beta_max), TCK[i]))/ds[i]                
        # penalties
        if (p[1]<Tmin):
            penalty *= (1.0+3.0*(Tmin-p[1]))
        if (p[1]>Tmax):
            penalty *= (1.0+3.0*(p[1]-Tmax))
        if (p[2]<beta_min):
            penalty *= (1.0+3.0*(beta_min-p[2]))
        if (p[2]>beta_max):
            penalty *= (1.0+3.0*(p[2]-beta_max))
            
    # print 'MBB_fit_chi2_fun --- DONE!'
    return penalty*delta




def Deriv_fun(p, delta, TCK, s, ds, tmin, tmax, beta1, beta2, fixed_T, fixed_beta, FIXED_BETA_ARRAY=False):
    # Works with three free parameters and with T or spectral index as fixed
    if (FIXED_BETA_ARRAY):
        # we have two free parameters but 2d TCK; fixed_beta = fixed beta in this pixel
        D   =  zeros((len(s), 2), float32)
        for i in range(len(s)):
            D[i,0]  =  -bisplev(p[1], fixed_beta, TCK[i])/ds[i]                # d/dI
            D[i,1]  =  -p[0] * bisplev(p[1], fixed_beta, TCK[i], dx=1)/ds[i]   # d/dT
    elif (len(p)==2):   # one parameter is fixed, T or beta
        D   =  zeros((len(s), 2), float32)
        for i in range(len(s)):
            D[i,0]  =  -splev(p[1],TCK[i])/ds[i]
            D[i,1]  =  -p[0] * splev(p[1],TCK[i], der=1)/ds[i]
    else:             # three free parameters
        D   =  zeros((len(s), 3), float32)
        for i in range(len(s)):
            D[i,0]  =         -bisplev(p[1], p[2], TCK[i])      /ds[i]
            D[i,1]  =  -p[0] * bisplev(p[1], p[2], TCK[i], dx=1)/ds[i]
            D[i,2]  =  -p[0] * bisplev(p[1], p[2], TCK[i], dy=1)/ds[i]
    return D
    

    

def MBB_fit_chi2_mp(UM, S, dS=None, fixed_T=None, fixed_beta=None, filters=None, 
                           Tmin=5.0, Tmax=29.0, beta_min=0.5, beta_max=5.0,
                          T_bins=200, beta_bins=80, ncpus=2, spline_smooth=2.0e-7):
    """
    Wrapper for MBB_fit_chi2 to run the fit with multiprocessing.
    Uses the standard routine MBB_fit_chi2.
    """
    # print '--- MBB_fit_chi2_parallel_2 ---'
    shape      = S[0].shape
    nfreq      = len(S)                 # number of frequencies
    npix       = len(ravel(S[0]))       # pixels in each vector/map
    n          = int(npix/ncpus)        # pixels per job
    if (n==0):                          # less pixels than workers
        return MBB_fit_chi2(UM, S, dS, fixed_T, fixed_beta, filters,
                            Tmin, Tmax, beta_min, beta_max, T_bins, beta_bins,
                            spline_smooth=spline_smooth)
    # set up jobs
    manager = multiprocessing.Manager()
    II      = multiprocessing.Array('f', npix)
    TT      = multiprocessing.Array('f', npix)
    BB      = multiprocessing.Array('f', npix)
    PROC = []
    print('MULTIPROCESSING WITH NCPUS %d' % ncpus)
    for i in range(ncpus):
        a, b = i*n, (i+1)*n+1
        if (i==(ncpus-1)): b = npix
        SS, dSS  = [], []
        for iw in range(nfreq):
            SS.append( ravel(S[iw])[a:b].copy() )
        if (dS==None):
            dSS = None
        else:
            for iw in range(nfreq):
                dSS.append( ravel(dS[iw])[a:b].copy() )
        # MBB_fit_chi2(UM, S, dS=None, fixed_T=None, fixed_beta=None, filters=None, 
        #              Tmin=5.0, Tmax=29.0, beta_min=0.5, beta_max=5.0,
        #              T_bins=200, beta_bins=80,
        #              result_T=None, result_B=None, result_ind=None)
        p =  multiprocessing.Process(target = MBB_fit_chi2,
                         args   = (UM, SS, dSS, fixed_T, fixed_beta, filters, 
                         Tmin, Tmax, beta_min, beta_max, 
                         T_bins, beta_bins,
                         II, TT, BB, a, spline_smooth))
            
        PROC.append(p)
        print('.... start %d ... %d-%d = %d pixels' % (i, a, b, b-a), SS[0].shape)
        p.start()
        
    # gather results
    # print 'gather'
    # T, B = zeros(npix, float32), zeros(npix, float32)
    for i in range(ncpus):
        PROC[i].join()
        print('gather ', i)
        
    II = asarray(II, float32)
    TT = asarray(TT, float32)
    BB = asarray(BB, float32)
    TT.shape = shape
    BB.shape = shape
    return II, TT, BB





def MBB_fit_LSQ_mp(UM, S, dS=None, fixed_T=None, fixed_beta=None, filters=None, cc_routine=None,
                           Tmin=6.9, Tmax=33.0, beta_min=0.59, beta_max=3.9,
                           T_bins=171, beta_bins=171, ncpus=2, spline_smooth=1.0e-6, TCK=[]):
    """
    Wrapper for MBB_fit_chi2 to run the fit with multiprocessing.
    Uses the standard routine MBB_fit_chi2 and MBB_fit_LSQ
    """
    FIXED_BETA_ARRAY = False
    if (fixed_beta!=None): 
        if (not(isscalar(fixed_beta))):
            FIXED_BETA_ARRAY = True  # beta fixed but using an array in fixed_beta
        
    # print '--- MBB_fit_chi2_parallel_2 ---'
    shape      = S[0].shape
    nfreq      = len(S)                 # number of frequencies
    npix       = len(ravel(S[0]))       # pixels in each vector/map
    n          = int(npix/ncpus)        # pixels per job
    if (n==0):                          # less pixels than workers
        return MBB_fit_LSQ(UM, S, dS, fixed_T, fixed_beta, filters, cc_routine,
                           Tmin, Tmax, beta_min, beta_max, T_bins, beta_bins,
                           spline_smooth=spline_smooth, TCK=TCK)
    # set up jobs
    manager = multiprocessing.Manager()
    II      = multiprocessing.Array('f', npix)
    TT      = multiprocessing.Array('f', npix)
    BB      = multiprocessing.Array('f', npix)
    PROC = []
    print('MULTIPROCESSING WITH NCPUS %d' % ncpus)
    for i in range(ncpus):
        a, b = i*n, (i+1)*n+1
        if (i==(ncpus-1)): b = npix
        SS, dSS  = [], []
        for iw in range(nfreq):
            SS.append( ravel(S[iw])[a:b].copy() )
        if (dS==None):
            dSS = None
        else:
            for iw in range(nfreq):
                dSS.append( ravel(dS[iw])[a:b].copy() )
        if (FIXED_BETA_ARRAY):
            FIXED_BETA = fixed_beta[a:b].copy()
            p =  multiprocessing.Process(target = MBB_fit_LSQ,
                         args   = (UM, SS, dSS, fixed_T, FIXED_BETA, filters, cc_routine,
                         Tmin, Tmax, beta_min, beta_max, 
                         T_bins, beta_bins,
                         II, TT, BB, a,
                         TCK, spline_smooth))
        else:
            p =  multiprocessing.Process(target = MBB_fit_LSQ,
                         args   = (UM, SS, dSS, fixed_T, fixed_beta, filters, cc_routine,
                         Tmin, Tmax, beta_min, beta_max, 
                         T_bins, beta_bins,
                         II, TT, BB, a,
                         TCK, spline_smooth))
            
        PROC.append(p)
        print('.... start %d ... %d-%d = %d pixels' % (i, a, b, b-a), SS[0].shape)
        p.start()
        
    # gather results
    # print 'gather'
    # T, B = zeros(npix, float32), zeros(npix, float32)
    for i in range(ncpus):
        PROC[i].join()
        print('gather ', i)
        
    II = asarray(II, float32)
    TT = asarray(TT, float32)
    BB = asarray(BB, float32)
    II.shape = shape
    TT.shape = shape
    BB.shape = shape
    return II, TT, BB



#########################################################################################
###########################  From MJ/Aux/ColourCorrect.py  ##############################
#########################################################################################


def CalculateColourCorrection(filterfile, f0, beta, T):
    """
    Calculate colour correction factor for modified black body curve BB(T)*f^beta
    based on given filter profile.
    Usage:
        k = CalculateColourCorrection(filterfile, f0, beta, T)
    Inputs:
        filterfile = name of an ascii file containing the filter profile,
                     the file should contain two columns, first giving 
                     wavelength [um], the second the response
                     (see, e.g., 'resp-p25.txt' in HOMEDIR+'/tt/PIA/beam/CCTABLES')
        f0         = nominal frequency of the filter [Hz]
        beta, T    = emissivity index and colour temperature for the modified 
                     black body
    Returns:
        k  =  the colour correction factor (divide by k to get true monochromatic values)
    Prerequisites:
        the filter file in directory HOMEDIR+'/tt/PIA/beam/CCTABLES' !!!
    """
    from MJ.mjDefs import HOMEDIR, C_LIGHT, PLANCK
    from numpy import isscalar
    # compare integrals of given spectrum and spectrum ~1/f and the
    # monochromatic value at the reference wavelength
    # what is the value at f0 if total energy is the same
    #  k =  ref(f0) -> int(ref) = int(this) => this(f0)
    #  k^-1 = this(f0)/int(this)   * int(ref)/ref(f0)
    #  k    = ref(f0)/int(ref)  * int(this)/this(f0)
    if (filterfile.find('/')<0):
        name = HOMEDIR+'/HERSCHEL/FILTER/filter_%s.txt' % filterfile
        if (not(os.path.exists(name))):
            name   =  HOMEDIR+'/tt/PIA/beam/CCTABLES/'+filterfile
        if (not(os.path.exists(name))): #  just in the current work directory
            name   =  filterfile
    else:
        name   =  filterfile
    # print name
    d      =  numpy.loadtxt(name)                        # filter response
    f, R   =  C_LIGHT/(1.0e-4*d[:,0]), d[:,1].copy()     # filter response
    if (isscalar(T)):
        Sthis  =  PlanckFunction(f, T)  * (f/1.0e12)**beta   # our spectrum
        val    =  PlanckFunction(f0, T) * (f0/1.0e12)**beta  # our value at nominal wavelength
        Ithis  =  simps(Sthis*R, f)                          # our integral over filter
        # reference spectrum
        Sref   =  1.0e12/f                                   # ref value at nominal wavelength
        Iref   =  simps(Sref*R, f)                           # integral of ref spectrum
        k      =  (Ithis / val)   /  ( Iref  / (1.0e12/f0) ) # colour correction factor
        return k
    else:
        if (isscalar(beta)):
            kk     =  zeros(len(T), float32)            
            vals   =  PlanckFunction(f0, T) * (f0/1.0e12)**beta  # our value at nominal wavelength
            # reference spectrum
            Sref    =  1.0e12/f                                  # ref value at nominal wavelength
            Iref    =  simps(Sref*R, f)                          # integral of ref spectrum
            for iii in range(len(T)):
                t       =  T[iii]
                Sthis   =  PlanckFunction(f, t)  * (f/1.0e12)**beta   # our spectrum
                # val     =  PlanckFunction(f0, t) * (f0/1.0e12)**beta  # our value at nominal wavelength
                val     =  vals[iii]
                Ithis   =  simps(Sthis*R, f)                          # our integral over filter
                kk[iii] =  (Ithis / val)   /  ( Iref  / (1.0e12/f0) ) # colour correction factor
            return kk
        else:
            # we need to loop over beta as well
            kk     =  zeros(len(T), float32)
            for iii in range(len(T)):
                t       =  T[iii]
                b       =  beta[iii]
                Sthis   =  PlanckFunction(f, t)  * (f/1.0e12)**b   # our spectrum
                val     =  PlanckFunction(f0, t) * (f0/1.0e12)**b  # our value at nominal wavelength
                Ithis   =  simps(Sthis*R, f)                       # our integral over filter
                # reference spectrum
                Sref    =  1.0e12/f                                # ref value at nominal wavelength
                Iref    =  simps(Sref*R, f)                        # integral of ref spectrum
                kk[iii] =  (Ithis/val) / (Iref / (1.0e12/f0) )     # colour correction factor
            return kk
            
            
            

def CalculateColourCorrectionGeneral(f0, filter_f, filter_value, spectrum_f, spectrum_value):
    """
    Calculate colour correction factor for a filter and a spectrum, both
    specified through input vectors (no files are used).
    Usage:
        k = CalculateColourCorrectionGeneral(f0, filter_f, filter_value, 
                                             spectrum_f, spectrum_value)
    Inputs:
        f0             =  nominal frequency of the filter [Hz].
        filter_f       =  vector of frequency values [Hz] covering the filter profile
        filter_value   =  relative response for frequencies in filter_f
        spectrum_f     =  vector of frequency values [Hz] covering the filter
        spectrum_value =  relative intensity of the spectrum at frequencies
                          listed in spectrum_f
    Output:
        k  =  the colour correction factor (divide by k to get true monochromatic values)
    Notes:
        Linear interpolation is used to get spectrum values onto filter_f frequency grid.
        If spectrum is smooth, spectrum_f grid needs not be dense but one must make sure
        it covers the filter profile.
        To get monochromatic values at f0, the instrumental values must be divided by
        the colour correction factor k.
    """
    # Compare integrals of given spectrum and spectrum ~1/f and the
    # monochromatic value at the reference wavelength.
    # What is the value at um0 if total energy is the same?
    # this spectrum
    f, R   =  filter_f, filter_value                   # filter response
    ip     =  interp1d(spectrum_f, spectrum_value, bounds_error=False, fill_value=0.0, kind='linear')    
    Sthis  =  ip(filter_f)                             # our spectrum at filter frequencies
    val    =  ip(f0)                                   # our spectrum at nominal wavelength
    Ithis  =  simps(Sthis*R, f)                        # our integral over filter profile
    # reference spectrum
    Sref   =  1.0e12/f
    Iref   =  simps(Sref*R, f)                             # reference spectrum integrated
    # print '>> ', Sthis[0:3], val, Ithis, Sref, Iref
    # print '      freq ', min(filter_f), max(filter_f), min(spectrum_f), max(spectrum_f), f0
    k      =  (Ithis / val)   /  ( Iref  / (1.0e12/f0) )
    return k
    


def CalculateColourCorrectionFilterfile(f0, filterfile, spectrum_f, spectrum_value, ref_expo=-1.0):
    """
    Calculate colour correction factor a given spectrum using a filter file.
    Usage:
        k = CalculateColourCorrectionFilterfile(f0, filterfile, spectrum_f, spectrum_value)
    Inputs:
        f0             =  nominal frequency of the filter [Hz].
        filterfile     =  name of the file containing [ um, response ]
        spectrum_f     =  vector of frequency values [Hz] covering the filter
        spectrum_value =  relative intensity of the spectrum at frequencies
                          listed in spectrum_f
        ref_expo       =  default -1.0,  the assumed reference spectrum powerlaw exponent
    Output:
        k  =  the colour correction factor (divide by k to get true monochromatic values)
    Notes:
        Linear interpolation is used to get spectrum values onto filter_f frequency grid.
        If spectrum is smooth, spectrum_f grid needs not be dense but one must make sure
        it covers the filter profile.
        To get monochromatic values at f0, the instrumental values must be divided by
        the colour correction factor k.
    """
    # read the filter profile
    d      =  loadtxt(filterfile)
    f      =  um2f(d[:,0])
    R      =  d[:,1].copy()
    # reverse the order so that frequency is increasing
    m      =  argsort(f)
    f      =  f[m].copy()
    R      =  R[m].copy()    
    # Compare integrals of given spectrum and spectrum ~1/f and the
    # monochromatic value at the reference wavelength.
    # What is the value at um0 if total energy is the same?
    # this spectrum
    # f, R   =  filter_f, filter_value                       # filter response
    ip     =  interp1d(spectrum_f, spectrum_value, bounds_error=False, fill_value=0.0)    
    Sthis  =  ip(f)                                        # our spectrum at filter frequencies
    val    =  ip(f0)                                       # our spectrum at nominal wavelength
    Ithis  =  simps(Sthis*R, f)                            # our integral over filter profile
    # reference spectrum
    Sref   =  1.0e12 * f**ref_expo
    Iref   =  simps(Sref*R, f)                             # reference spectrum integrated
    # print '>> ', Sthis[0:3], val, Ithis, Sref, Iref
    # print '      freq ', min(filter_f), max(filter_f), min(spectrum_f), max(spectrum_f), f0
    k      =  (Ithis / val)   /  ( Iref  / (1.0e12/f0) )
    return k
    

