#!/usr/bin/python
 
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
     
