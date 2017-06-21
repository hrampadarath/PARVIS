#!/usr/bin/env python


"""
SYNOPSIS

    TODO helloworld [-h,--help] [-v,--verbose] [--version]

DESCRIPTION

    TODO This describes how to use this script. This docstring
    will be printed by the script if there is an error or
    if the user requests help (-h or --help).

EXAMPLES

    TODO: Show some examples of how to use this script.

EXIT STATUS

    TODO: List exit codes

AUTHOR

    Hayden Rampadarath <hayden.rampadarath@manchester.ac.uk && haydenrampadarath@gmail.com>

LICENSE

    This script is in the public domain, free from copyrights or restrictions.

VERSION

    $Id$

NOTES
	This file was generated using the LoadPythonTemplate.py script using the pythonTemplate.py in ~/scripts
	Convert an array geocentric XYZ coordinates to local coordinates (i.e. relative to centre of array)
	

"""




#Some common imported modules - uncomment as required

import sys,os,traceback, optparse
import time
import re
import matplotlib.pyplot as plt
import numpy as np


#Define global constants


#Some useful Global args
#arcmin2rad = np.pi / (180 * 60)
#arcsec2rad = np.pi / (180 * 60 * 60)
#rad2arcmin = 1/(arcmin2rad)
#rad2arcsec = 1/(arcsec2rad)

# Figures defining the WGS84 ellipsoidal description of the Earth
R = 6378137.0 # equatorial radius (m)
f = 1/298.257223563 # flattening
# A useful derived value
#b = (1. - f)**2
e = 1-(1-f)**2

def main():
	outfile = 'VLA_C_hor_xyz_v1.txt'
	f = open(outfile,'ab')

	config = 'vla.c.cfg'
	ants_xyz = np.genfromtxt(config)
	x = []
	y = []
	z = []
	ref_XYZ = np.array([ants_xyz[:,0].mean(),ants_xyz[:,1].mean(),ants_xyz[:,2].mean()])
	for i in range (len(ants_xyz)):
		XYZ = ants_xyz[i]
		#1. geo long
		Long = np.arctan(XYZ[1]/XYZ[0])
		#print 'Long: ', np.degrees(Long)
		Lat = convergeU2(XYZ)
		#print 'Lat: ', np.degrees(Lat)
		Alt = altitude(Lat,XYZ)
		#print 'Alt: ', np.degrees(Alt)
		xi,yi,zi = relativeVLAAntPositions(ref_XYZ,Long,Lat,Alt)
		x.append(xi)
		y.append(yi)
		z.append(zi)
		f.write(str(xi)+' '+str(yi)+' '+str(zi))  
		f.write('\n')
	xyz = np.array([x,y,z])	

	just_xyz = np.genfromtxt('VLA_C_hor_xyz_v2')
	plt.scatter(just_xyz[:,0],just_xyz[:,1],marker='o',c='b')	
	plt.scatter(x,y,marker='o',c='r')
	plt.show()

	
	



def convergeU2(XYZ):

	s = np.sqrt(XYZ[0]**2 + XYZ[1]**2)
	B1 = np.arctan(XYZ[2]/(1-f)*s)
	U1 = lon(B1, XYZ)
	#print 'U1: ', U1
	#-----------------------
	B2 = np.arctan(((1-f)*np.sin(U1))/np.cos(U1))	
	U2 = lon(B2,XYZ)
	#print 'U2: ', U2
	#-----------------------
	B3 = np.arctan(((1-f)*np.sin(U2))/np.cos(U2))	
	U3 = lon(B3,XYZ)
	#print 'U3: ', U3
	#-----------------------
	B4 = np.arctan(((1-f)*np.sin(U3))/np.cos(U3))	
	U4 = lon(B4,XYZ)
	#print 'U4: ', U4
	#-----------------------
	return U4




def lon(B, XYZ):
	#2. initial guess of geo long, U1
	s = np.sqrt(XYZ[0]**2 + XYZ[1]**2)
	A1 = XYZ[2] + ((e**2*(1-f))/(1-e**2))*R*(np.sin(B))**3
	A2 = s - e**2 * R*(np.cos(B))**3
	U  = np.arctan(A1/A2)
	return U




def altitude(U,XYZ):
	"""
	evaluate the altitude above the planetary ellipsoid
	"""

	s = np.sqrt(XYZ[0]**2 + XYZ[1]**2)
	N = R/np.sqrt(1-e**2*(np.sin(U))**2)
	h = s * np.cos(U) + (XYZ[2] + e**2 * N * np.sin(U))*np.sin(U) - N
	return h




def relativeVLAAntPositions(ref_xyz,Long,Lat,Alt):
	#define centre of VLA array:
	#XYZ_cent = np.array([-1601.389e+03,-5042.634e+03,3553.876e+03])
	XYZ_cent = ref_xyz
	Long_cent = np.arctan(XYZ_cent[1]/XYZ_cent[0])
	Lat_cent = convergeU2(XYZ_cent)
	Alt_cent= altitude(Lat_cent,XYZ_cent)

	#x = deltaLong(Long,Long_cent,Lat_cent)
	#y = deltaLat(Lat,Lat_cent)
	x = np.cos(Lat_cent)*(Long - Long_cent) * R
	y = (Lat - Lat_cent) * R
	deltaAlt = Alt - Alt_cent
	return x, y, deltaAlt




def deltaLat(Lat,Lat_cent):
	
	deltaLat = Lat - Lat_cent
	a = np.sin(deltaLat/2)**2
	c = np.arctan2(np.sqrt(a),np.sqrt(1-a))
	d = R * c
	return d


def deltaLong(Long,Long_cent,Lat_cent):
	deltaLong = (Long - Long_cent)
	a = np.cos(Lat_cent)**2 * np.sin(deltaLong/2.)**2
	c = np.arctan2(np.sqrt(a),np.sqrt(1-a))
	d = R * c
	return d

	



if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'], version='$Id$')
        parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
        (options, args) = parser.parse_args()
        #if len(args) < 1:
        #    parser.error ('missing argument')
        if options.verbose: print time.asctime()
        main()
        if options.verbose: print time.asctime()
        if options.verbose: print 'TOTAL TIME IN MINUTES:',
        if options.verbose: print (time.time() - start_time) / 60.0
        sys.exit(0)
    except KeyboardInterrupt, e: # Ctrl-C
        raise e
    except SystemExit, e: # sys.exit()
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        os._exit(1)	


