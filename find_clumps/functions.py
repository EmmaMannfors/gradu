#!/usr/bin/python

##############################################
# Copied from subtract_average folder
##############################################
# Added dec/ra_to_deg
##############################################
from astropy.io import fits
#import aplpy
import numpy as np
from math import atan2, sin, cos, sqrt
from numpy import array, arctan2, sin, cos, sqrt
from Aux import DEGREE_TO_RADIAN, RADIAN_TO_DEGREE
 
def distance_on_sphere(lon1, lat1, lon2, lat2):
    """
    Calculates angular distance between two points on sky.
    Usage:
        d = distance_on_sphere(lon1, lat1, lon2, lat2)
    Inputs:
        lon1, lat1, lon2, lat2 = longitude and latitude of the two points [radians]
    Returns:
        distance on sphere [radians]
    """
    dlon = lon1-lon2
    nom  = (cos(lat1)*sin(dlon))**2.0 + (cos(lat2)*sin(lat1)-sin(lat2)*cos(lat1)*cos(dlon))**2.0
    den  = sin(lat2)*sin(lat1) + cos(lat2)*cos(lat1)*cos(dlon)
    return arctan2(sqrt(nom),den)
     
 
def distance_on_sphere_2(lon1, lat1, lon2, lat2):
    """
    Calculates approximate angular distance between two points on sky. Valid
    for small distances only.
    Usage:
        d = distance_on_sphere_2(lon1, lat1, lon2, lat2)
    Inputs:
        lon1, lat1, lon2, lat2 = longitude and latitude of the two points [radians]
    Returns:
        angular distance on sphere [radians]
    """
    return sqrt( ((lon1-lon2)*cos(lat1))**2.0 + (lat1-lat2)**2.0 )
 
 
def distance_on_sphere_deg(lon1, lat1, lon2, lat2):
    """
    As distance_on_sphere but with arguments in degrees.
    Returns distance in degrees.
    """
    return RADIAN_TO_DEGREE*distance_on_sphere(lon1*DEGREE_TO_RADIAN, lat1*DEGREE_TO_RADIAN, 
           lon2*DEGREE_TO_RADIAN, lat2*DEGREE_TO_RADIAN)
  

def ra_to_deg(ra):
	g = ra.split(':')
	h = float(g[0])
	m = float(g[1])
	s = float(g[2])
	h += (m/60.0)
	h += (s/3600.0)
	deg = h*(360.0/24.0)
	return round(deg,6)

def dec_to_deg(dec):
	y = dec.split(':')
	deg,amin,asec = float(y[0]),float(y[1]),float(y[2])
	if deg < 0:
		deg -= (amin/60.0)
		deg -= (asec/3600.0)
		return round(deg,6)
	deg += (amin/60.0)
	deg += (asec/3600.0)
	return round(deg,6) 




















