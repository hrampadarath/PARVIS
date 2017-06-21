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





#############Primary beam codes####################


def makeJonesMatrix(fov,nu,cell_size,primary_beam):

	"""
	SKA DVA 1 primary beam
	generate and store the Mueller terms
	fov - in radians
	pol - required polarisation, XX or YY only
	nu - obsering frequency (MHz)
	ff_size - sampling 
	"""
	

	#TODO - need to incorporate beam selection for frequencies between the above
	# specified freqs. Most important for wide bandwidths.

	#load the primary beams and convert from matlab format to numpy
	

	
	#freq1 = np.array([660,670,680,690,700,710,720,730,740,750])
	finc = 10
	freq1 = np.arange(660,760,finc)
	beam1 = primary_beam

	
	if (nu >= freq1.min()-finc/2) and (nu < freq1.max()+finc/2):
		fx = (np.abs(freq1-nu))
    		ind = np.where(fx==fx.min())
    		if len(ind[0]) >= 2:
        		f = ind[0][1:][0]
		else:
			f = ind[0][0]
        		
	
	mat = scipy.io.loadmat(beam1)
	
	beamfreq = mat['freq'][0]
	#print 'beam frequency (MHz) = ', beamfreq[f]
	
	#load the Jones matrix.
	#Jones Matrix = [Jpv Jph; Jqv Jqh], where each Jxy is a 3 dimensional matrix of size length(th) x length(ph) x length(freq)
	#Jpv is the vertical co-polarized fields and Jph the vertical cross-polarized fields.
	# Jqh and Jqv are the horizontal co- and cross-polarized components
	Jxy = mat['Jqv'] 
	Jx = mat['Jqh']
	Jy = mat['Jpv']
	Jyx = mat['Jph']

	
	#create a 2D grid over the required fov to interpolate the beams onto
	#these are simply the x and y positions of the Jones and Mueller terms in the imnage plane. 
	#This is required to obtain the value of the primary beam at the location of the source

	ff_size=int(np.degrees(fov)*3600/cell_size)
	x1 = np.linspace(-fov/2.,fov/2.,ff_size) #np.linspace(-0.13,0.13)
	y1 = x1
	x,y = np.meshgrid(x1,y1)

	phands = 1.
	chands = 1.

	Jx  = generateJones(Jx[0::1,1:,f]*phands,fov,mat,cell_size)
	Jy  = generateJones(Jy[0::1,1:,f]*phands,fov,mat,cell_size)

	Jxy = generateJones(Jxy[0::1,1:,f]*chands,fov,mat,cell_size)
	Jyx = generateJones(Jyx[0::1,1:,f]*chands,fov,mat,cell_size)
		
	#print 'Jx/Jxy ', 10*np.log10(np.amax(np.abs(Jx))/np.amax(np.abs(Jxy)))
	#print 'Jy/Jyx', 10*np.log10(np.amax(np.abs(Jy))/np.amax(np.abs(Jyx)))

	return x, y, Jx, Jxy, Jyx, Jy


def generateJones(J1,fov,mat,cell_size):


	ff_size=int(np.degrees(fov)*3600/cell_size)
	x1 = np.linspace(-fov/2.,fov/2.,ff_size) #np.linspace(-0.13,0.13)
	y1 = x1
	x,y = np.meshgrid(x1,y1)

	J1s = np.fliplr(J1)
	Jp = np.concatenate((J1s,J1),axis=1)
	
	th = mat['th'][:,0]*deg2rad #theta angle 0 -> 8 degrees convert to radians
	ph1 = mat['ph'][:,0][1:]*deg2rad #the phi angle goes from 0 to 180. Assuming bilateral symmetry
	
	#The beams are assumed to be bilateral symmetrical (i.e. left-right)
	#To obtain the full beam extend the range of phi to -180 to 180
	ph2 = ph1[::-1]*-1#reverse and multiply by -1
	ph = np.concatenate((ph2,ph1),axis=0)

	
	#determine the phi and theta values of the grid
	
	pph = np.arctan2(x,y) #+ rotate
	pth = (np.arccos( np.cos(x) * np.cos(y) ))
	
	###Since the beams are in spherical coordinates, to plot as a 2D image in python
	###I had to do an interpolation. 
	
	#interpolate the beams over a rectangular mesh.
	bcplx = np.zeros((len(pth),len(pth)),dtype=complex)
	
	interpolatorReal = interpolate.RectBivariateSpline(th,ph,Jp.real)
	pJreal = interpolatorReal.ev(pth,pph)

	bcplx.real = pJreal
	
	interpolatorImag = interpolate.RectBivariateSpline(th,ph,Jp.imag)
	pJimag = interpolatorImag.ev(pth,pph)
	bcplx.imag = pJimag

	return bcplx
	


def MuellerMatrix(x,y,Jones,pol):


	Jx, Jxy, Jyx, Jy = normaliseJones(Jones[0], Jones[1], Jones[2], Jones[3])

	#The full Mueller matrix for each of the 4 polarisation	
	#M_xx = (Mxx, Mxxy, Mxyx, Mxyxy) = [Jx*np.conj(Jx), Jx*np.conj(Jxy), Jxy*np.conj(Jx), Jxy*np.conj(Jxy)]
	#M_yy = (Myxyx, Myxy, Myyx, Myy) = [Jy*np.conj(Jy), Jyx*np.conj(Jy), Jy*np.conj(Jyx), Jy*np.conj(Jy)]
	#M_xy = (Mxyx, Mxy, Mxyyx, Mxyy) = [Jx*np.conj(Jyx), Jx*np.conj(Jy), Jxy*np.conj(Jyx), Jxy*np.conj(Jy)]
	#M_yx = (Myxx, Myxxy, Myx, Myxy) = [Jyx*np.conj(Jx), Jyx*np.conj(Jxy), Jy*np.conj(Jx), Jy*np.conj(xy)]  
	

	if pol == 'XX':		
		#calculate the Mueller terms related to the Xpol
		M0_cplx = np.abs(Jx*np.conj(Jx))
		M0 = np.abs(M0_cplx)

		M1_cplx = np.abs(Jx*np.conj(Jxy))
		M1 = np.abs(M1_cplx)

		M2_cplx = np.abs(Jxy*np.conj(Jx))
		M2 = np.abs(M2_cplx)

		M3_cplx = np.abs(Jxy*np.conj(Jxy))
		M3 = np.abs(M3_cplx)


	elif pol == 'XY':
		#calculate the Mueller terms related to the XYpol

		M0_cplx = np.abs(Jx*np.conj(Jyx))
		M0 = np.abs(M0_cplx)

		M1_cplx = np.abs(Jx*np.conj(Jy))
		M1 = np.abs(M1_cplx)

		M2_cplx =  np.abs(Jxy*np.conj(Jyx)) 
		M2 = np.abs(M2_cplx)

		M3_cplx =  np.abs(Jxy*np.conj(Jy))
		M3 = np.abs(M3_cplx)

	elif pol == 'YX':

		#calculate the Mueller terms related to the YXpol

		M0_cplx =  np.abs(Jyx*np.conj(Jx))
		M0 = np.abs(M0_cplx)

		M1_cplx =  np.abs(Jyx*np.conj(Jxy))
		M1 = np.abs(M1_cplx)

		M2_cplx =  np.abs(Jy*np.conj(Jx))
		M2 = np.abs(M2_cplx)

		M3_cplx =  np.abs(Jy*np.conj(Jxy))
		M3 = np.abs(M3_cplx)


	elif pol == 'YY':
		
		#calculate the Mueller terms related to the Ypol

		M0_cplx =  np.abs(Jyx*np.conj(Jyx))
		M0 = np.abs(M0_cplx)

		M1_cplx =  np.abs(Jyx*np.conj(Jy))
		M1 = np.abs(M1_cplx)

		M2_cplx =  np.abs(Jy*np.conj(Jyx))
		M2 = np.abs(M2_cplx)

		M3_cplx =  np.abs(Jy*np.conj(Jy))
		M3 = np.abs(M3_cplx)


	return x, y, M0, M1, M2, M3

def normaliseJones(Jx, Jxy, Jyx, Jy):

	#find the determinant of the central pixel of the 2 x 2 Jones matrix 
	xc = np.shape(Jx)[0]/2
	yc = np.shape(Jx)[1]/2

	det = Jx[xc,yc] * Jy[xc,yc] - Jxy[xc,yc] * Jyx[xc,yc]
	
	#determine the inverse of teh Jones matrix at the central pixel
	invXX = Jy[xc,yc]/det
	invXY = -Jxy[xc,yc]/det
	invYX = -Jyx[xc,yc]/det
	invYY = Jx[xc,yc]/det
	
	#normalise the Jones matrix by the inverse of the JOnes matrix at the central pixel
	#i.e. invJones*Jones

	normJx  = invXX * Jx  + invXY * Jyx
	normJxy = invXX * Jxy + invXY * Jy
	normJyx = invYX * Jx  + invYY * Jyx
	normJy  = invYX * Jxy + invYY * Jy 

	#print 'maxVal normJx ', np.amax(np.abs(normJx))
	#print 'maxVal normJxy ', np.amax(np.abs(normJxy))	
	#print 'maxVal normJyx ', np.amax(np.abs(normJyx))
	#print 'maxVal normJy ', np.amax(np.abs(normJy))

	return normJx, normJxy, normJyx, normJy


def findMuellerValue(xsrc,ysrc,Mterms):

	"""
	To obtain the beam attenuation of a source located at l = xsrc and m = ysrc
	for the sph beam
	xsrc - xposition of source in radians
	ysrc - ypositon of source in radians
	Mterms - four Mueller terms (Mterms[2:6]) as related to a single polarisation
		including the x (Mterms[0]) and y (Mterms[1]) positions
	"""

	l = find_nearest(Mterms[0],xsrc)
	m = find_nearest(Mterms[1],ysrc)
	loc = np.ravel(np.where((Mterms[0]==l) & (Mterms[1]==m)))
	M0 = Mterms[2][loc[0]][loc[1]]
	M1 = Mterms[3][loc[0]][loc[1]]
	M2 = Mterms[4][loc[0]][loc[1]]
	M3 = Mterms[5][loc[0]][loc[1]]
	return M0,M1,M2,M3



########unused but important codes


def findsphbeamvalue(xsrc,ysrc,sbeam):

	"""
	To obtain the beam attenuation of a source located at l = xsrc and m = ysrc
	for the sph beam
	xsrc - radians
	ysrc - radians
	rotate - paralactic rotation in radians
	pol - polarisation XX or YY
	nu - frequency of obs in MHz
	"""
	#if np.abs(xsrc) >= np.abs(ysrc):
	#	fov=(np.abs(xsrc)+0.1*np.abs(xsrc))*2
	#elif np.abs(ysrc) > np.abs(xsrc):
	#	fov=(np.abs(ysrc)+0.1*np.abs(ysrc))*2

	#sbeam = sphbeam(fov,pol,nu)

	l = find_nearest(sbeam[0],xsrc)
	m = find_nearest(sbeam[1],ysrc)
	loc = np.ravel(np.where((sbeam[0]==l) & (sbeam[1]==m)))
	return sbeam[3][loc[0]][loc[1]]


def find_nearest(array,value):
	if np.abs(value) <= np.max(np.abs(array)):
		idx = (np.abs(array-value)).argmin()
	else:
		print 'Source is beyond the current FOV of the beam'
		print np.abs(value),
		print np.max(np.abs(array))
		sys.exit()
    	return np.ravel(array)[idx]




def sphbeam(fov,pol,nu):

	"""
	SKA DVA 1 primary beam
	fov - in radians
	rotate - in radians
	pol - required polarisation, XX or YY only
	nu - obsering frequency (MHz)
	ff_size - sampling 
	"""
	

	#TODO - need to incorporate beam selection for frequencies between the above
	# specified freqs. Most important for wide bandwidths.

	#load the primary beams and convert from matlab format to numpy
	

	
	#freq1 = np.array([660,670,680,690,700,710,720,730,740,750])
	finc = 10
	freq1 = np.arange(660,760,finc)
	beam1 = 'Data1.mat'

	
	if (nu >= freq1.min()-finc/2) and (nu < freq1.max()+finc/2):
		fx = (np.abs(freq1-nu))
    		ind = np.where(fx==fx.min())
    		if len(ind[0]) >= 2:
        		f = ind[0][1:][0]
		else:
			f = ind[0][0]
        		
	
	mat = scipy.io.loadmat(beam1)
	
	beamfreq = mat['freq'][0]
	#print 'beam frequency (MHz) = ', beamfreq[f]

	th = mat['th'][:,0]*deg2rad #theta angle 0 -> 8 degrees convert to radians
	ph1 = mat['ph'][:,0][1:]*deg2rad #the phi angle goes from 0 to 180. Assuming bilateral symmetry
	
	#The beams are assumed to be bilateral symmetrical (i.e. left-right)
	#To obtain the full beam extend the range of phi to -180 to 180
	ph2 = ph1[::-1]*-1#reverse and multiply by -1
	ph = np.concatenate((ph2,ph1),axis=0)
	
	
	#load the Jones matrix.
	#Jones Matrix = [Jpv Jph; Jqv Jqh], where each Jxy is a 3 dimensional matrix of size length(th) x length(ph) x length(freq)
	#Jpv is the vertical co-polarized fields and Jph the vertical cross-polarized fields.
	# Jqh and Jqv are the horizontal co- and cross-polarized components
	Jqv = mat['Jqv'] 
	Jqh = mat['Jqh']
	Jpv = mat['Jpv']
	Jph = mat['Jph']


	if pol == 'XX':
		J1 = Jqh[0::1,1:,f]
	elif pol == 'YY':
		J1 = Jpv[0::1,1:,f]


	J1s = np.fliplr(J1)
	Jp = np.concatenate((J1s,J1),axis=1)
	
	
	###Since the beams are in spherical coordinates, to plot as a 2D image in python
	###I had to do an interpolation. 
	
	#create a 2D grid over the required fov to interpolate the beams onto
	ff_size=int(np.degrees(fov)*3600/100.)
	x1 = np.linspace(-fov/2.,fov/2.,ff_size) #np.linspace(-0.13,0.13)
	y1 = x1
	x,y = np.meshgrid(x1,y1)
	
	#determine the phi and theta values of the grid
	
	pph = np.arctan2(x,y) #+ rotate
	pth = (np.arccos( np.cos(x) * np.cos(y) ))
	

	
	#interpolate the beams over a rectangular mesh.
	bcplx = np.zeros((len(pth),len(pth)),dtype=complex)
	
	interpolatorReal = interpolate.RectBivariateSpline(th,ph,Jp.real)
	pJreal = interpolatorReal.ev(pth,pph)

	bcplx.real = pJreal
	
	interpolatorImag = interpolate.RectBivariateSpline(th,ph,Jp.imag)
	pJimag = interpolatorImag.ev(pth,pph)
	bcplx.imag = pJimag

	Acplx = bcplx * np.conj(bcplx)
	
	pJ = np.abs(Acplx)/np.max(np.abs(Acplx))

	return x,y,Acplx,pJ







def Gaussbeam(pol,nu):

	"""
	A perfect Gaussian beam
	fov - field of view in radians
	pol - polarisation. THis have no effect on the beam. 
	nu - frequency in MHz
	"""

	d = 15 # 15m dish
	#freq = 700e6
	c = 299792458
	lam = c/(nu*1e6)
	fov = 450 * arcmin2rad
	#Obtain 1-D gaussian beam
	ff_size=1000#int(np.degrees(fov)*3600/1000.)
	#ff_size=int(np.degrees(fov)*3600/100.)
	#print ff_size
	#ff_size=int(fov/100.)
	x2= np.linspace(-fov/2,fov/2,ff_size)#distance from pb centre in arcmins on the x axi
	y2= x2[:,np.newaxis]#np.arange(-60,60,0.12)#distance from pb centre in arcmins on the y axis
	x,y = np.meshgrid(x2,y2)
	

	theta = 1.02*(lam/d)#*rad2arcmin
	#print 'FWHM = ',theta
	thetal = theta/(2*np.sqrt(2*np.log(2)))
	#Amax = l**2/FWHM
	G = np.exp(- (x**2 + y**2) /(2*thetal**2))
	return x,y,G




def findGaussbeamvalue(xsrc,ysrc,sbeam):

	"""
	To obtain the beam attenuation of a source located at l = xsrc and m = ysrc
	for the sph beam
	xsrc - radians
	ysrc - radians
	rotate - paralactic rotation in radians
	pol - polarisation XX or YY
	nu - frequency of obs in MHz
	"""
	#if np.abs(xsrc) >= np.abs(ysrc):
	#	fov=(np.abs(xsrc)+0.1*np.abs(xsrc))*2
	#elif np.abs(ysrc) > np.abs(xsrc):
	#	fov=(np.abs(ysrc)+0.1*np.abs(ysrc))*2

	#sbeam = Gaussbeam(pol,nu)
	l = find_nearest(sbeam[0],xsrc)
	m = find_nearest(sbeam[1],ysrc)
	loc = np.ravel(np.where((sbeam[0]==l) & (sbeam[1]==m)))
	return sbeam[2][loc[0]][loc[1]]


