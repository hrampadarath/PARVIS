import re,sys
import numpy as np
import scipy.io
from scipy import interpolate


#some globals
deg2rad    = np.pi / 180
arcmin2rad = np.pi / (180 * 60)
arcsec2rad = np.pi / (180 * 60 * 60)
rad2arcmin = 1/(arcmin2rad)
rad2arcsec = 1/(arcsec2rad)


def parse_inp(filename):
	''' Parse the list of inputs given in the specified file. (Modified from evn_funcs.py, taken from eMERLIN_pipeline.py)'''
	INPUTFILE = open(filename, "r")
	control = dict()

	# a few useful regular expressions
	newline = re.compile(r'\n')
	space = re.compile(r'\s')
	char = re.compile(r'\w')
	comment = re.compile(r'#.*')

	# parse the input file assuming '=' is used to separate names from values
	for line in INPUTFILE:
		if char.match(line):
			line = comment.sub(r'', line)
			line = line.replace("'", '')
			(param, value) = line.split('=')

			param = newline.sub(r'', param)
			param = param.strip()
			param = space.sub(r'', param)

			value = newline.sub(r'', value)
			value = value.strip()
			valuelist = value.split(', ')
			control[param] = valuelist

	return control



def rotate2d(degrees,point,origin):
	"""
	A rotation function that rotates a point around a point
	to rotate around the origin use [0,0]
	"""
	degrees = degrees * 1
	x = point[0] - origin[0]
	yorz = point[1] - origin[1]
	newx = (x*np.cos(np.radians(degrees))) - (yorz*np.sin(np.radians(degrees)))
	#print newx
	newyorz = (x*np.sin(np.radians(degrees))) + (yorz*np.cos(np.radians(degrees)))
	#print newyorz
	
	newx += origin[0]
	newyorz += origin[1]
	
	return newx,newyorz





def float2str(inarray):
	strs = []
	for n in inarray:
		strs.append(str(n))
	return strs



def timetoHA(vis_time, vis_date, Array, RA):
    """
    Converting julian date to HA. 
    Copied from 
    http://uk.mathworks.com/matlabcentral/fileexchange/28233-convert-eci-to-ecef-coordinates/content/JD2GMST.m
    and http://www.csgnetwork.com/siderealjuliantimecalc.html
    and http://www.cv.nrao.edu/~rfisher/Ephemerides/times.html#LMST
    returns the HA 

    """

    #print 'converting normal time to HA'
    #print Array.lon
    Long = np.degrees(float((Array.lon))) # array longitude
    
    JD = vis_date + vis_time
    #1. determine the Julian date of the previous midninght, JD0
    JD0_max = np.floor(JD)+0.5
    if JD >= JD0_max:
        JD0 = JD0_max
    elif JD < JD0_max:
        JD0 = np.floor(JD)-0.5
    H = (JD-JD0)*24.      	#Time in hours past previous midnight
    D  = JD - 2451545.0    #Compute the number of days since J2000
    D0 = JD0 - 2451545.0   #Compute the number of days since J2000
    T  = D/36525.          #Compute the number of centuries since J2000
    #Calculate GMST in hours (0h to 24h)
    GMST = ((6.697374558 + 0.06570982441908*D0  + 1.00273790935*H + 0.000026*(T**2)) % 24. )
    LMST = GMST + (Long/15.) # determine LMST in hrs
    HA = LMST-RA
        
    return HA

def parangle(HA, dec_d, lat_d):
    
    	"""
    	from p. 91 of  treatise on spherical astronomy" By Sir Robert Stawell Ball
    	this is how CASA does the pb rotation (thanks to Preshant)
    	HA - hourangle in hours
    	lat_d - latitude of array in degrees
    	dec_d - declination in degrees
    	returns:
    	paralctic angle in degrees
    	"""
    	lat = np.radians(lat_d)
    	dec = np.radians(dec_d)   
    	HA = np.radians(HA*15.)
    	sin_eta_sin_z = np.cos(lat)*np.sin(HA)
    	cos_eta_sin_z = np.sin(lat)*np.cos(dec) - np.cos(lat)*np.sin(dec)*np.cos(HA)
    	
    	eta = np.arctan2(sin_eta_sin_z,cos_eta_sin_z)
    	return np.degrees(eta)







