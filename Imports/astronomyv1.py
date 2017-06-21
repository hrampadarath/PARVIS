#! /usr/bin/env python

######################################################################
# Module containing a set of astronomical conversion functions
######################################################################

import math
import numpy as np
from scipy.special import i0,i0e
from math import *



def newpos(dra,ddec,rah0,ram0,ras0,decd0,decm0,decs0):
	"""
	To determine a new position given"
	dra - offset in ra in arcseconds
	ddec - offset in dec arcseconds
	raho, ramo, raso - ra of old position
	decdo, decmo, decso - dec of old position  
	"""
	ra0 = hmstora(rah0,ram0,ras0)
	dec0 = dmstodec(decd0,decm0,decs0)
	#1 calculate the new dec
	decd1 = dec0 + (ddec/3600.)
	dec1 = dectodms(decd1)
	#2. calculate the new ra
	deltaD = math.cos(math.radians((dec0+decd1)/2.0))
	rad1 = (dra/3600.)/deltaD + ra0
	ra1 = ratohms(rad1)
	#print 'new RA =%0.10f degrees; Dec = %0.10f degrees' % (rad1,decd1)
	#print ra1,
	#print dec1
	#return rad1,decd1
	return ra1, dec1


def deltapos(rah1,ram1,ras1,decd1,decm1,decs1,rah2,ram2,ras2,decd2,decm2,decs2):
	"""To calculate the difference between 2 positions given in the tradiational
	Ra and dec format
	"""
	#position 1 in degs
	ra1 = hmstora(rah1,ram1,ras1)
	dec1 = dmstodec(decd1,decm1,decs1)
	#position 2 in degs
	ra2 = hmstora(rah2,ram2,ras2)
	dec2 = dmstodec(decd2,decm2,decs2)
	#1 calculate the separation in delta in arcsecs
	deltadec = deltasep(dec1,dec2)
	#2 calculate separation in ra arcsecs
	deltara = alphasep(ra1,ra2,dec1,dec2)
	#3 calculate the positional difference arcsecs
	sourcesep = angsep(ra1,dec1,ra2,dec2)
	#print "RA position difference (arcsecs) = %0.10f" %(deltara)
	#print "Dec position difference(arcsecs) = %0.10f" %(deltadec)
	#print "Source position difference (arcsecs)= ", sourcesep
	return deltara,deltadec,sourcesep


def ratohms(radegs):
	"""Convert RA in decimal degrees format to hours, minutes,
	seconds format.

	Keyword arguments:
	radegs -- RA in degrees format

	Return value:
	ra -- list of 3 values, [hours,minutes,seconds]

	"""

	rah = int(radegs/15)
	ram = int((radegs/15-rah)*60)
	ras = (((radegs/15-rah)*60)-ram)*60
	ra = [rah,ram,ras]
	return ra
# Converts a decimal RA to hms format	

def dectodms(decdegs):
	"""Convert Declination in decimal degrees format to hours, minutes,
	seconds format.

	Keyword arguments:
	decdegs -- Dec. in degrees format

	Return value:
	dec -- list of 4 values, [sign (+/-),degrees,minutes,seconds]

	"""
	
	sign = "+"
	if decdegs<0:
		sign = "-"
		decdegs *= -1
	decd = int(decdegs)
	decm = int((decdegs-decd)*60)
	decs = (((decdegs-decd)*60)-decm)*60
	dec = [sign,decd,decm,decs]
	return dec
# Converts a decimal dec to dms format

def hmstora(rah,ram,ras):
	"""Convert RA in hours, minutes, seconds format to decimal
	degrees format.

	Keyword arguments:
	rah,ram,ras -- RA values (h,m,s)

	Return value:
	radegs -- RA in decimal degrees

	"""
	radegs = 15*(float(rah)+(float(ram)/60)+(float(ras)/3600.0))
	
	return radegs
# Converts an hms format RA to decimal degrees

def dmstodec(decd,decm,decs):
	"""Convert Dec in degrees, minutes, seconds format to decimal
	degrees format.

	Keyword arguments:
	decd,decm,decs -- list of Dec values (d,m,s)

	Return value:
	decdegs -- Dec in decimal degrees

	"""
	sign = float(decd)
	if sign>0:
		decdegs = float(decd)+(float(decm)/60)+(float(decs)/3600.0)
	else:
		decdegs = float(decd)-(float(decm)/60)-(float(decs)/3600.0)
	
	return decdegs
# Converts a dms format Dec to decimal degrees	


def angsep(ra1,dec1,ra2,dec2):
        """Find the angular separation of two sources, in arcseconds,
	using the proper spherical trig formula

	Keyword arguments:
	ra1,dec1 - RA and Dec of the first source, in decimal degrees
	ra2,dec2 - RA and Dec of the second source, in decimal degrees

	Return value:
	angsep - Angular separation, in arcseconds

	"""
    
	b = (math.pi/2)-math.radians(dec1)
	c = (math.pi/2)-math.radians(dec2)

	return 3600*math.degrees(math.acos((math.cos(b)*math.cos(c))+(math.sin(b)*math.sin(c)*math.cos(math.radians(ra1-ra2)))))

# Find angular separation of 2 positions, in arcseconds

def alphasep(ra1,ra2,dec1,dec2):
    	"""Find the angular separation of two sources in RA, in arcseconds

	Keyword arguments:
	ra1,dec1 - RA and Dec of the first source, in decimal degrees
	ra2,dec2 - RA and Dec of the second source, in decimal degrees

	Return value:
	angsep - Angular separation, in arcseconds

	"""

	return 3600*(ra1-ra2)*math.cos(math.radians((dec1+dec2)/2.0))

# Find angular separation in RA of 2 positions, in arcseconds

def deltasep(dec1,dec2):
	"""Find the angular separation of two sources in Dec, in arcseconds

	Keyword arguments:
	dec1 - Dec of the first source, in decimal degrees
	dec2 - Dec of the second source, in decimal degrees

	Return value:
	angsep - Angular separation, in arcseconds

	"""

	return 3600*(dec1-dec2)

# Find angular separation in Dec of 2 positions, in arcseconds

def alpha(l,m,alpha0,delta0):
	"""Convert a coordinate in l,m into an coordinate in RA

	Keyword arguments:
	l,m -- direction cosines, given by (offset in cells) x cellsi (radians)
	alpha_0, delta_0 -- centre of the field

	Return value:
	alpha -- RA in decimal degrees

	"""
	return alpha0+(math.degrees(math.atan2(l,((math.sqrt(1-(l*l)-(m*m))*math.cos(math.radians(delta0)))-(m*math.sin(math.radians(delta0)))))))

# Find the RA of a point in a radio image, given l,m and field centre, for the SIN projection

def delta(l,m,alpha0,delta0):
	"""Convert a coordinate in l,m into an coordinate in Dec

	Keyword arguments:
	l,m -- direction cosines, given by (offset in cells) x cellsi (radians)
	alpha_0, delta_0 -- centre of the field

	Return value:
	delta -- Dec in decimal degrees

	"""
	return math.degrees(math.asin(m*math.cos(math.radians(delta0))+(math.sqrt(1-(l*l)-(m*m))*math.sin(math.radians(delta0)))))

# Find the declination of a point in a radio image, given l,m and field centre, for the SIN projection

def alpha_ncp(l,m,alpha0,delta0):
	"""Convert a coordinate in l,m into an coordinate in RA

	Keyword arguments:
	l,m -- direction cosines, given by (offset in cells) x cellsi (radians)
	alpha_0, delta_0 -- centre of the field

	Return value:
	alpha -- RA in decimal degrees

	"""
	return  alpha0+math.degrees(math.atan2(l,(math.cos(math.radians(delta0))-(m*math.sin(math.radians(delta0))))))

# Find the RA of a point in a radio image, given l,m and field centre,
# for the NCP projection

def delta_ncp(l,m,alpha0,delta0):
	"""Convert a coordinate in l,m into an coordinate in Dec

	Keyword arguments:
	l,m -- direction cosines, given by (offset in cells) x cellsi (radians)
	alpha_0, delta_0 -- centre of the field

	Return value:
	delta -- Dec in decimal degrees

	"""
	if delta0<0:
		sign = -1
	else:
		sign = +1

	alpha = alpha0+math.degrees(math.atan2(l,(math.cos(math.radians(delta0))-(m*math.sin(math.radians(delta0))))))
		
	return sign*math.degrees(math.acos((math.cos(math.radians(delta0))-(m*math.sin(math.radians(delta0)))).math.cos(math.radians(alpha-alpha0))))

# Find the declination of a point in a radio image, given l,m and
# field centre, for the NCP projection

def l(ra,dec,cra,cdec,incr):
	"""Convert a coordinate in RA,Dec into a direction cosine m

	Keyword arguments:
	ra,dec -- Source location
	cra,cdec -- centre of the field
	incr -- number of degrees per pixel (negative in the case of RA)

	Return value:
	l -- Direction cosine

	"""
	return (math.cos(math.radians(dec))*math.sin(math.radians(ra-cra)))/(math.radians(incr))

# Find the l direction cosine in a radio image, given an RA and Dec and the
# field centre

def m(ra,dec,cra,cdec,incr):
	"""Convert a coordinate in RA,Dec into a direction cosine m

	Keyword arguments:
	ra,dec -- Source location
	cra,cdec -- centre of the field
	incr -- number of degrees per pixel

	Return value:
	m -- direction cosine

	"""
	return ((math.sin(math.radians(dec))*math.cos(math.radians(cdec)))-(math.cos(math.radians(dec))*math.sin(math.radians(cdec))*math.cos(math.radians(ra-cra))))/math.radians(incr)

# Find the m direction cosine in a radio image, given an RA and Dec and the
# field centre, for the SIN projection

def m_ncp(ra,dec,cra,cdec,incr):
	"""Convert a coordinate in RA,Dec into a direction cosine m

	Keyword arguments:
	ra,dec -- Source location
	cra,cdec -- centre of the field
	incr -- number of degrees per pixel

	Return value:
	m -- direction cosine

	"""
	return((math.cos(math.radians(cdec))-(math.cos(math.radians(dec))*math.cos(math.radians(ra-cra))))/math.sin(math.radians(cdec)))

# Find the m direction cosine in a radio image, given an RA and Dec and the
# field centre, for the NCP projection

def eq_to_gal(ra,dec):
	"""Find the Galactic co-ordinates of a source given the equatorial
	co-ordinates

	Keyword arguments:
	(alpha,delta) -- RA, Dec in decimal degrees

	Return value:
	(l,b) -- Galactic longitude and latitude, in decimal degrees

	"""

	R = [[-0.054875539726,-0.873437108010,-0.483834985808],[0.494109453312,-0.444829589425,+0.746982251810],[-0.867666135858,-0.198076386122,+0.455983795705]]

	s = [math.cos(math.radians(ra))*math.cos(math.radians(dec)),math.sin(math.radians(ra))*math.cos(math.radians(dec)),math.sin(math.radians(dec))]
	
	sg = []
	
	sg.append(s[0]*R[0][0]+s[1]*R[0][1]+s[2]*R[0][2])
	sg.append(s[0]*R[1][0]+s[1]*R[1][1]+s[2]*R[1][2])
	sg.append(s[0]*R[2][0]+s[1]*R[2][1]+s[2]*R[2][2])

	b = math.degrees(math.asin(sg[2]))
	l = math.degrees(math.atan2(sg[1],sg[0]))

	if l<0:
		l = l+360

	return (l,b)

# Return the Galactic co-ordinates of a point given in equatorial co-ordinates

def gal_to_eq(l,b):
	"""Find the Galactic co-ordinates of a source given the equatorial
	co-ordinates

	Keyword arguments:
	(l,b) -- Galactic longitude and latitude, in decimal degrees

	Return value:
	(alpha,delta) -- RA, Dec in decimal degrees

	"""

	Rinv = [[-0.054875539692115144, 0.49410945328828509, -0.86766613584223429], [-0.87343710799750596, -0.44482958942460415, -0.19807638609701342], [-0.4838349858324969, 0.74698225182667777, 0.45598379574707293]]

	sg = [math.cos(math.radians(l))*math.cos(math.radians(b)),math.sin(math.radians(l))*math.cos(math.radians(b)),math.sin(math.radians(b))]
	
	s = []
	
	s.append(sg[0]*Rinv[0][0]+sg[1]*Rinv[0][1]+sg[2]*Rinv[0][2])
	s.append(sg[0]*Rinv[1][0]+sg[1]*Rinv[1][1]+sg[2]*Rinv[1][2])
	s.append(sg[0]*Rinv[2][0]+sg[1]*Rinv[2][1]+sg[2]*Rinv[2][2])

	dec = math.degrees(math.asin(s[2]))
	ra = math.degrees(math.atan2(s[1],s[0]))

	if ra<0:
		ra = ra+360

	return (ra,dec)

# Return the Galactic co-ordinates of a point given in equatorial co-ordinates

def ergstowatts(ergs):

	"""To convert ergs/second/Hertz to Watts/Hertz: 1 Watt = 1e7 ergs/s"""
	watts = ergs/(1.0e7)
	print 'Watts/second =' , '%10.4e ' % watts
	#return watts

def wattstoergs(watts):

	"""To convert ergs/second/Hertz to Watts/Hertz: 1 Watt = 1e7 ergs/s"""
	ergs = watts*(1.0e7)
	print 'ergs/second =' , '%10.4e ' % ergs
	return ergs


def distflux(ly,L):

	"""To calculate the Jy/Hz of a source given the distance in lightyears (ly) and luminosity (L) in W/Hz"""

	r = ly*9.4605284e15#m
	Jy = 1e-26#W/m**2/Hz
	S = ((L/(4*np.pi*r**2))/Jy)/1.0e-3
	print 'Flux in mJy at %02d lightyears = %02.10f' %(ly,S)
	#print S

def distflux2(pc,L):

	"""To calculate the mJy/Hz of a source given the distance in  parsecs and luminosity (L) in W/Hz"""

	r = pc*3.08568025e16
	Jy = 1e-26#W/m**2/Hz
	S = ((L/(4*np.pi*r**2))/Jy)/1.0e-3
	print 'Flux in mJy at %02d pc = %02.4f' %(pc,S)
	return S

def distlum(d,S):

	"""To calculate the luminosity of a source given the distance in parsecs and flux in Jy"""

	r = d*3.08568025e16 #m
	Jy = 1e-26#W/m**2/Hz
	s = S*Jy
	L = s*(4*np.pi*r**2)
	print 'Luminosity at %02d Parsecs = %02.4f W/Hz' %(d,L)
	return L
	#print S


def flux2lum(ly,s):

	""" to determine the luminosity of an isotropic emitter given it flux in mJy and distance in ly"""
	r = ly*9.4605284e15#m
	Jy = 1e-26#W/m**2/Hz
	S = s*1e-3#flux in Jy
	L = (4*np.pi*r**2)*(S*Jy)
	print 'Luminosity (W/Hz) of a %02.4f Jy source at  %02d lightyears= %10.4e ' %(s,ly,L)


def flux2lumpc(pc,s):

	""" to determine the luminosity of an isotropic emitter given it flux in mJy and distance in pc"""
	r = pc*3.08568025e16 #m
	Jy = 1e-26#W/m**2/Hz
	S = s*1e-3#flux in Jy
	L = (4*np.pi*r**2)*(S*Jy)
	#print 'Luminosity (W/Hz) of a %02.4f mJy source at  %02d pc = %10.4e ' %(s,pc,L)
	return L



def angdist(ly,d):
	"""To determine the angular distance of a source in arcseconds assuming no cosmology. The radius of the source should be given in AU and the distnce to source in light years"""
	r = ly*9.4605284e15
	S = d*1.49598e11#Au in metres
	theta = (S/r)*(180/np.pi)*3600.
	print 'Angular distance of source at %02d light years, and diameter of %02d AU = %0.24f arcseconds' %(ly,d,theta)

def angdia(D,r):
	"""Determine the angular diameter of a source in arcseconds (witout cosmological considerations)
	given its distance D in parsecs and its radius,r in AU
	"""
	Parsec = 3.08568025e16 #meters
	print Parsec
	AU = 1.49598e11 #meters
	Dm = 2.*D*Parsec
	dm = 2.*r*AU
	delta = 2.*atan2(dm,Dm) * (180./pi) * 3600.
	return delta

def angdia_pc(D,r):
	"""Determine the angular diameter of a source in arcseconds (witout cosmological considerations)
	given its distance D and its size/diameter in parsecs
	"""
	Parsec = 3.08568025e16 #meters
	print Parsec
	Dm = D*Parsec
	dm = r/2*Parsec
	delta = np.arctan2(dm,Dm) * (180./np.pi) * 3600.
	print "size of source/seperation of 2 sources at %0.24f parsecs, and physical diameter of %0.24f parcecs = %0.24f arcsecs" %(D,r,delta)


def phdia(angdia,D):
	"""Determine the physical diameter of a source or physical distance between 2 sources located on a 2-D plane in parsecs, given the angular diamter in arcseconds and the distance to the source in parsecs"""	
	Parsec = 3.08568025e16	#m
	Dm = D*Parsec
	angrad = (angdia/2)*(np.pi/(180*3600))
	print angrad
	r = 2 * (Dm*np.tan(angrad))
	rparsecs = r/Parsec
	print 'Physical dimeter of source at %0.24f parsecs, and angular diameter of %0.24f arcseconds = %0.24f parsecs' %(D,angdia,rparsecs)
	
def dmstodec2(sign,decd,decm,decs):
	"""Convert Dec in degrees, minutes, seconds format to decimal
	degrees format.

	Keyword arguments:
	decd,decm,decs -- list of Dec values (d,m,s)

	Return value:
	decdegs -- Dec in decimal degrees

	"""
	if sign>0:
		decdegs = float(decd)+(float(decm)/60)+(float(decs)/3600.0)
	else:
		decdegs = float(decd)-(float(decm)/60)-(float(decs)/3600.0)
	
	return decdegs
# Converts a dms format Dec to decimal degrees with sign
"""
def sline(x1,y1,x2,y2):
	#to determine the parameters of a straight line (y=mx+c) given 2 points

	p1 = (float(x1),float(y1))
	p2 = (float(x2),float(y2))

	grad = (p1[1]-p2[1])/(p1[0]-p2[0])
	yin = p2[1]-grad*p2[0]
		
	print grad
	#print "c=",c
#"""



def dphidt(phi,V,dudt):
	"""
	calculate the maximum dphidt for a source offset from the phase centre by phi milliarcseconds
	"""
	c = 3e8
	l = ((phi*1e-3)/3600)
	v = V*1e6
	dpdt = ((2*np.pi*v)/c)*l*dudt
	print 'maximum dpdt for a signal locates at %0.24f mas = %0.24f' %(phi, dpdt)

def SpecInd(v1,S1,v2,S2):
	"""
	Calculating the spectral index of source between 2 frequencies in GHz (v1, and v2) and the corresponding 
	flux (S1,S2) in Jy. Note the assumptions are S propo v**(alpha)
	"""

	alpha = (log(S1/S2))/(log(v1/v2))
	print 'Spectral Index of source with S1 = %0.2f Jy @ v1 = %0.4f GHz and  S2 = %0.2f Jy @ v2 = %0.4f GHz = %0.2f' %(S1,v1,S2,v2,alpha)
	return alpha

def res(nu,D):
	"""
	Calculate the resolution of a telescope/array im arcseconds, given the observing frequency im MHz and maximum baseline in km
	"""
	nu=nu*1e6
	c = 3e8
	la = c/nu
	print 'Lambda =',la
	d = D*1e3
	res = ((la/d)*(180/np.pi))*3600
	print 'Max. Resolution = ', res  

def angsep2(ra1,dec1,ra2,dec2):
	"""Find the angular separation of two sources in Dec, in arcseconds

	Keyword arguments:
	dec1 - Dec of the first source, in decimal degrees
	dec2 - Dec of the second source, in decimal degrees

	Return value:
	angsep - Angular separation, in arcseconds"""

	rad1 = np.pi/180. #radians to 1 degree
	ra1 = ra1*rad1
	dec1 = dec1*rad1
	ra2 = ra2*rad1
	dec2 = dec2*rad1

	sep1 = np.sqrt(((math.cos(2*dec2)+1)/2.)*(1-math.cos(2*(ra2-ra2)))/2. + (math.cos(dec1)*math.sin(dec2)-math.sin(dec1)*math.cos(dec2)*math.cos(ra2-ra2))**2)
	sep2 = math.sin(dec1)*math.sin(dec2) + math.cos(dec1)*math.cos(dec2)*math.cos(ra2-ra1)
	angsep = ((180/np.pi)*math.tan(sep1/sep2))*3600.
	return '%0.24f' %(angsep)



def brightnessTemp(Fnu,nu,bmaj):
	"""
	Estimateing the brightness temperature of a radio source
	nu   = frequency in MHz
	Fnu  = flux of source in Jy
	bmaj = size of major axis in arcseconds
	D    = distance to source in parsecs  
	"""
	c  = 299792458#light speed in m/s
	k  = 1.3806488e-23#boltzmann constant (m**2 kg /s**2/K)
	Jy = 1e-26#W/m**2/Hz
	Parsec = 3.08568025e16 #meters
	#D = D * Parsec
	nu  = nu*1e6
	Fnu = Fnu*Jy
	r = (bmaj/2)/3600. * (np.pi/180.)
	Tb = (c**2 * Fnu * 1**2)/(2 * np.pi * r**2 *k *nu**2)
	print "The Brightness Tempertaure of a %0.4f Jy source, at %0.4f MHz = %0.4f K" % (Fnu/Jy, nu/1e6, Tb)
	return Tb


def ThermalSen(SEFD,npol,N,tint,dnu):

	"""
	Estimate the thermal sensitivity of a homogenous array given:
	SEFD of the antenna
	nc - correlator efficiency (0.92 VLA; 0.7 EVN)
	npol - number of pols
	N - number of antennas
	tint - total integration on-source time, in hrs
	dnu - total bandwidth, MHz
	"""
	tint = tint * 60 * 60
	dnu  = dnu * 1e6
	nc = 0.92
	sigma = SEFD / (nc * np.sqrt(npol * N*(N-1)*tint*dnu)) 
	print '1-sigma thermal noise sensitivity for given configuration =', sigma


def ztod(z):
	"""
	Determine the distance in pc given the redshift
	"""
	h0=71#km/s/Mpc
	c=3e6#km/s
	return ((z*c)/h0)

def dtoz(d):
	"""
	Determine the redshift givebn the distanc in Mpc
	"""
	h0=71.#km/s/Mpc
	c=3.e6#km/s
	return (h0*d)/c

def smAngle(maxB,l):
	"""
	To calculate the small angle approximation for an array.
	i.e. the point at which w-projection needs to be considered.
	inputs -
	maxB = maximum baseline of the array in km (and most affected buy the w-term
	l = wavelength in m
	other terms -
	phi_e = maximum acceptable phase error. typically taken as 0.1 radians (i.e. 
	"""
	phi_e = 0.1#radians
	theta = l/(maxB*1e3)
	field_max  = ((phi_e * theta)/np.pi)**0.5 * (180/np.pi)*60
	print 'Maximum radial distance from phase centre that can be imaged with this array, assuming a maximum phase error of 0.1 radians is %0.2f arcmins' %(field_max)

def cellsize(maxB,l):
	"""
	Estimate the minimum cell size to use during imaging 
	given the maximum baseline (maxB) in km and the observing wavelength (l) in m
	"""
	cells = ((l/(maxB*1e3)) * (180/np.pi) * 3600 )/ 5 
	print 'the recommended cellsize is %0.2f arcseconds' %(cells)
	return cells


############################################
def parangle(HA, DEC, lat):
    """
    +
     NAME:
         parangle

     PURPOSE:
         To compute the parallactic angle at a given position on the sky.
	taken from http://www.mpia-hd.mpg.de/homes/ianc/python/_modules/spec.html#parangle

     CATEGORY:
         Spectroscopy    

     CALLING SEQUENCE:
         eta, za = parangle(HA, DEC, lat)

     INPUTS:
         HA  - Hour angle of the object, in decimal hours (0,24)
         DEC - Declination of the object, in degrees
         lat - The latitude of the observer, in degrees

     KEYWORD PARAMETERS:
         CANCEL - Set on return if there is a problem

     OUTPUTS:
         eta - The parallactic angle
         za  - The zenith angle

     PROCEDURE:
         Given an objects HA and DEC and the observers latitude, the
         zenith angle and azimuth are computed.  The law of cosines
         then gives the parallactic angle.  

     EXAMPLE:
         NA


     MODIFICATION HISTORY:
         2000-04-05 - written by M. Cushing, Institute for Astronomy,UH
         2002-08-15 - cleaned up a bit.
         2003-10-21 - changed to pro; outputs zenith angle as well - WDV
         2011-10-07 17:58 IJMC: Converted to Python
-"""

  #pro parangle, HA, DEC, lat, eta, za, CANCEL=cancel

    cancel = 0
    d2r = np.deg2rad(1.)
    r2d = np.rad2deg(1.)
    #"""
    #  If HA equals zero then it is easy.
    HA = HA % 24
    #  Check to see if HA is greater than 12.
    if hasattr(HA, '__iter__'):
        HA = np.array(HA, copy=False)
        HAind = HA > 12
        if HAind.any():
            HA[HAind] = 24. - HA[HAind]
    else:
        if HA>12.:
            HA = 24. - HA
    #"""
    HA = HA*15.

    #  Determine Zenith angle and Azimuth
    cos_za = np.sin(lat*d2r) * np.sin(DEC*d2r) + \
             np.cos(lat*d2r) * np.cos(DEC*d2r) * np.cos(HA*d2r)
    za     = np.arccos(cos_za) * r2d
    cos_az = (np.sin(DEC*d2r) - np.sin(lat*d2r)*np.cos(za*d2r)) / \
             (np.cos(lat*d2r) * np.sin(za*d2r))
    az     = np.arccos(cos_az)*r2d

    if hasattr(az, '__iter__'):
        azind = az==0
        if azind.any() and DEC<lat:
            az[azind] = 180.
    else:
        if az==0. and DEC<lat:
            az = 180.

    tan_eta = np.sin(HA*d2r)*np.cos(lat*d2r) / \
              (np.cos(DEC*d2r)*np.sin(lat*d2r) - \
               np.sin(DEC*d2r)*np.cos(lat*d2r)*np.cos(HA*d2r))
    eta     = np.arctan(tan_eta)*r2d

    if hasattr(eta, '__iter__'):
        etaind = eta < 0
        ezind = (eta==0) * (az==0)
        zaind = za > 90
        if etaind.any():
            eta[etaind] += 180.
        elif ezind.any():
            eta[ezind] = 180.
        if zaind.any():
            eta[zaind] = np.nan
    else:
        if eta < 0:
            eta = eta#+= 180.
        elif eta==0 and az==0:
            eta = 180.
        if za>90:
            eta = np.nan

    HA = HA/15.0

    return eta, za

