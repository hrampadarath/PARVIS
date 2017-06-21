#!/usr/bin/env python


"""
SYNOPSIS

    TODO helloworld [-h,--help] [-v,--verbose] [--version]

DESCRIPTION
    PARVIS -  Python skA Radio Visibility fIts Simulator 

    To simulate a full radio interferometric observation and export to FITSidi format to be 
    imaged within CASA.
    Can include a primary beam.  Currently ony the SKA1MID primary beam can be used, 
    although it is not difficult to include others
   
    The SKA-MID PB is only generated for XX, YY, YX and XY
    and at freqs - 	#in MHz for the SKA sim the options are 
    	      	#freq1 = np.array([660,670,680,690,700,710,720,730,740,750])
    		#freq2 = np.array([690,691,692,693,694,695,696,697,698,699])
    
    The simulator is executed as 
    >python2.7 PARVIS.py  doFITSsim.inputs
    
    where doFITSsim.inputs is an input file that specifies the different 
    observational parameters
    Note before running this script, a new .xml file in generated via 
    >python2.7 createxml.py doFITSsim.inputs
    
    this is to create the FITS header file
    the full simulator can be executed using the SimulatorMueller.sh script
    
    Update: Now can include the full Mueller matrix into the simulation and simulate full polarisation

    This version is a more general version of the previous simulator createDVA1fitsidi-fullMueller
    And can be used to simulate observations using the VLA/JVLA and SKA1MID

EXAMPLES
    >python2 PARVIS.py -f PARVIS.inputs
 
    
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

"""




#Some common imported modules - uncomment as required

import traceback, optparse
import time
import re
import numpy as np
import os, sys, datetime, time
import math
from numpy.random import normal
import matplotlib.pyplot as plt

sys.path.append('./Imports')
from pyFitsidi import *
from astroCoords import *
from funcs import *
from astronomyv1 import deltasep, alphasep,hmstora,dmstodec

start = time.time()

deg2rad    = np.pi / 180
rad2deg    = 1/deg2rad
arcmin2rad = np.pi / (180 * 60)
arcsec2rad = np.pi / (180 * 60 * 60)
rad2arcmin = 1/(arcmin2rad)
rad2arcsec = 1/(arcsec2rad)

global earth_radius, light_speed, pi
earth_radius = 6371 * 10**3     # Earth's radius
light_speed = 299792458         # Speed of ligh

#define directories

models = './models'
antconfig = './antconfig'
pb = './primary_beams'
tools = './tools'

def main(inFile):

	print '------Using inputs specified in', inFile
	if os.path.isfile(inFile) :
		control = parse_inp(inFile)
	else :
		print "Error:" + inFile + "does not exist, quitting."
		sys.exit()	


	#first create the .xml file for the header
	try:
		os.system('python2.7 '+tools+'createxml.py '+inFile)
	except (RuntimeError, TypeError, NameError):
		print "Error creating xml file"


	####specify params from the .inputs file
	#notes 		= control[notes][0]
	array 		= control['array'][0]
	diameter	= control['diameter'][0]


	RAcent     = control['RAcent'][0].replace(',',' ').split()
	DECcent    = control['DECcent'][0].replace(',',' ').split()
	
	lat 		= float(control['lat'][0])
	Tobs   		= float(control['Tobs'][0])
	tstep  		= float(control['tstep'][0])
	trange 		= np.arange(0,Tobs,tstep/60.)
	t_len		= int(len(trange))


	######FREQUENCY & POLARISATION
	pol_len 	= int(control['pol_len'][0]) #number of pols
	ri_len  	= int(control['ri_len'][0])#number of data points
	chan_width 	= float(control['chan_width'][0]) #MHz
	ref_freq  	= float(control['start_freq'][0])# in MHz assume this is the start freq in the 1st channel
	chan_len 	= int(control['chan_len'][0]) #number of channels

	npol  		= int(control['pol_len'][0])
	nband 		= int(control['nband'][0])
	nchan 		= int(control['chan_len'][0])
	
	######ADVANCED OPTIONS
	doPB 		= control['doPB'][0] #simulate antenna voltage pattern??
	wterm 		= control['wterm'][0] # turn wterm on?
	derotation 	= control['derotation'][0] # derotate the primary beam?
	derotate	= float(control['derotate'][0])


	######FILES	

	outfolder = control['outfolder'][0]
	fitsfile = outfolder+control['fitsfile'][0]
	configxml = outfolder+control['configxml'][0]
	src_list = models+control['srclist'][0]
	array_config = antconfig+control['array_config'][0]

	array_geometry = np.genfromtxt(array_config)
	N_ants = len(array_geometry) # number of antennas
	bl_len = int(N_ants*(N_ants-1)/2.)
	
	DIR = control['DIR'][0]
	amper_array = DIR + control['amper_array'][0]
	
	
	pointing = {'RAcent': RAcent,
		    'DECcent': DECcent,
		    'lat':lat}
	params = np.array([t_len,chan_len,pol_len,bl_len])
	freq_params = np.array([ref_freq,chan_width,nband,nchan])
	advanced = {'doPB':doPB,
                    'wterm':wterm,
                    'derotation':derotation,
                    'derotate':derotate}
	
	
	print '################################################'
	print '		------ Observation specs -------- '
	print '         Array = ', array
	print '		RA = ', RAcent 
	print '         Dec = ', DECcent
	print '		Total observing time = %s hrs' % str(Tobs)
	print '		scan length = %s seconds' % str(tstep*60)
	print '		number of polarisations = ', pol_len
	print '		Starting frequency = %s MHz' % str(ref_freq)
	print '		channel width = %s MHz' % str(chan_width)
	print '		number of channels = ', chan_len
	print '		numer of IFs = ', nband
	print '		Simulating primary beam?: ', doPB
	print '		Simulating the wterm?: ', wterm
	print '         Derotate the pb?: ', derotation
	if derotation == 'T':
		print '         derotation every %s minutes' % str(derotate)
	print '		output file name: ', fitsfile
	print '		output directory: ', outfolder
	print '		array config: ', array_config
	print '		source list: ', src_list
	print '###################################################'



	"""
	Main function call. This is the conductor.
	"""
	
	print('\nConfiguring the output fits-idi')
	
	
	print('\nConfiguring Array geography')
	print('--------------------------')	
		
	# Let's make ourselves an Array (pyEphem observer)
	print('Location: '+ array)		
	OBSERVATORY = getObservatory(array,array_geometry)

	
	print('\nConfiguring phase source')
	print('--------------------------')
	# The source is our phase centre for UVW coordinates
	ImS = makeSource(
	    name="ImS",
	    ra = RAcent[0]+':'+RAcent[1]+':'+RAcent[2],
            dec = DECcent[0]+':'+DECcent[1]+':'+DECcent[2],flux='0')

	source = ImS
	source.compute(OBSERVATORY)
	print "Name: %s \nRA: %s \nDEC: %s"%(source.name,source.ra,source.dec)
	
	# Make a new blank FITS HDU
	print('\nCreating PRIMARY HDU')
	print('------------------------------------')
	hdu = make_primary(config=configxml)
	#print hdu.headerascardlist()
	
	# Go through and generate required tables
	print('\nCreating ARRAY_GEOMETRY')
	print('------------------------------------')
	# TODO: edit the config file
	tbl_array_geometry = make_array_geometry(config=configxml, num_rows=N_ants)# pyfitidi
	tbl_array_geometry = config_array_geometry(tbl_array_geometry,array_geometry,diameter)
	#print tbl_array_geometry['STABXYZ']
	
	print('\nCreating FREQUENCY')
	print('------------------------------------')
	tbl_frequency = make_frequency(config=configxml, num_rows=1)
	tbl_frequency = config_frequency(tbl_frequency,freq_params)
	#print tbl_frequency.header.ascardlist()
	print('\n')
	
	print('\nCreating SOURCE')
	print('------------------------------------')
	tbl_source = make_source(config=configxml, num_rows=1)
	tbl_source = config_source(tbl_source, source)
	#print tbl_source.header.ascardlist()
	print('\n')
	
	print('\nCreating ANTENNA')
	print('------------------------------------')
	tbl_antenna = make_antenna(config=configxml, num_rows=N_ants)
	tbl_antenna = config_antenna(tbl_antenna,array)
	#print tbl_antenna.header.ascardlist()
	print('\n')
	
	
	#################################
	print('\nCreating UV_DATA')
	print('------------------------------------')
	
	print('Data dimensions: %i dumps, %i chans, %i baselines, %i pols, %i data (real/imag)'\
	%(t_len, chan_len, bl_len, pol_len, ri_len))
	
	print('Generating blank UV_DATA rows...')
	tbl_uv_data = make_uv_data(config=configxml, num_rows=t_len*bl_len)
	
	print('Now filling FITS file with numpy generated visibilities ...')
	# The config function is in a seperate file, so import it
	#tbl_uv_data = config_uv_data(h5,tbl_uv_data, medicina, source)
	tbl_uv_data = config_uv_data(params, tbl_uv_data, OBSERVATORY, source, pointing, trange, freq_params, advanced, array_config, src_list, t_len, amper_array)
	
	#print tbl_uv_data.header.ascardlist()
	print('\n')
	
	##############################
	
	
	hdulist = pf.HDUList(
	            [hdu, 
	            tbl_array_geometry,
	            tbl_frequency,
	            tbl_antenna,
	            tbl_source, 
	            tbl_uv_data
	            ])
	
	print('Verifying integrity...')            
	hdulist.verify()
	
	if(os.path.isfile(fitsfile)):
	  print('Removing existing file...')
	  os.remove(fitsfile)
	print('Writing to file... %s ') % fitsfile
	hdulist.writeto(fitsfile)
	
	
	end = time.time()
	print 'elapsed time = %0.4f mins' % ((end - start)/60.)
	
	print('Done.')







####################Define funcs and classes



def getObservatory(array,array_geometry):

	arrayList = {'VLA':{'Latitude':'34:04:43.497', 'Longitude':'-107:37:05.819', 'Elevation':2124},
		      'KAT-7':{'Latitude':'-30:43:15.6','Longitude':'21:24:39.6', 'Elevation':1100}}
	
	now = '2016/3/23 20:0:0' #datetime.datetime.now()

	if array in arrayList.keys():
		longitude = arrayList[array]['Longitude']
		latitude  = arrayList[array]['Latitude']
		elevation = arrayList[array]['Elevation']

	obs = Array(lat=latitude, long=longitude,elev=elevation,date=now, antennas=array_geometry)
	return obs

class Array(ephem.Observer):
	""" An antenna array class.
	
	Based on pyEphem's Observer class.
	Probably very similar to the one in AIPY.
	
	Parameters
	----------
	lat: dd:mm:ss
	  latitude of array centre, e.g. 44:31:24.88
	long: dd:mm:ss
	  longitude of array centre, e.g. 11:38:45.56
	elev: float
	  elevation in metres of array centrem e.g. 28.0
	date: datetime object
	  date and time of observation, e.g. datetime.now()
	antennas: np.array([x,y,z])
	  numpy array of antenna positions, in xyz coordinates in meters,
	  relative to the array centre.
	    
	"""
	def __init__(self, lat, long, elev, date, antennas):
		super(Array, self).__init__()
		self.lat = lat
		self.long = long
		self.elev = elev
		self.date = date
		self.antennas = antennas
	def update(self, date):
		"""Update antenna with a new datetime"""
		self.date = date






def makeSource(name,ra,dec,flux=0,epoch=2000):
	""" Create a pyEphem FixedBody
	
	Parameters
	----------
	name: string
	  Name of source, e.g. CasA
	ra: hh:mm:ss
	  right ascension, e.g. 23:23:26
	dec: dd:mm:ss
	  declination e.g. 58:48:22.21
	flux: float
	  flux brightness in Jy (not actually used here)
	epoch: J2000
	  Defaults to J2000, i.e. 2000"
	"""
	line = "%s,f,%s,%s,%s,%d"%(name,ra,dec,flux,epoch)
	body = ephem.readdb(line)
	return body


def config_array_geometry(tbl, antenna_array,diameter):
	"""  Configures the array_geometry table with Medicina values
	
	Parameters
	----------
	tbl: pyfits.hdu
	  table to be configured
	antenna_array: np.array
	  an array of xyz coordinates of the antenna locations (offsets) in METERS
	  from the array centre (this is a keyword in the header unit)
	  use ECF2localConfig.py to convert from geocentric to local reference positions. 
	  e.g. 
	"""
	
	geometry = tbl.data
	
	# X-Y-Z in metres
	xyz_m = antenna_array

	for i in range(0,tbl.data.size):
		geometry[i]['ANNAME']  = 'VLA_%i'%i
		geometry[i]['STABXYZ'] = xyz_m[i]
		#print xyz_m[i],
		#print i
		geometry[i]['DERXYZ']  =  0
		#geometry[i]['ORBPARM'] = 0
		geometry[i]['NOSTA']   = i
		geometry[i]['MNTSTA']  = 1 
		# NOTE: Aperture arrays are given code 6, but not supported by CASA
		geometry[i]['STAXOF']  = np.array([0,0,0])
		geometry[i]['DIAMETER'] = diameter #change for different arrays
		
		tbl.data = geometry
	
	return tbl


def config_frequency(tbl,freq_params):
	"""
	Configures the frequency table.
	#TODO: needs reconfiguring for multiple frequency setups ie. multiple bands 
	
	Parameters
	----------
	tbl: pyfits.hdu
	  table to be configured
	"""
	
	ref_freq = freq_params[0]
	chan_width = freq_params[1]	
	nband = freq_params[2]
	nchan = freq_params[3]


	frequency = tbl.data[0]
	
	frequency['FREQID']         = 1
	frequency['BANDFREQ']       = 0         # This is offset from REF_FREQ, so zero!
	frequency['CH_WIDTH']       = chan_width * 10**6
	frequency['TOTAL_BANDWIDTH']= chan_width * 10**6 * nchan * nband

	frequency['SIDEBAND']       = 1
	
	tbl.data[0] = frequency
	
	return tbl  

def config_source(tbl, source):
	"""  Configures the source table.
	
	Parameters
	----------
	tbl: pyfits.hdu
	  table to be configured
	source: ephem.fixedBody
	  source to be phased to (use makeSource())
	"""
	
	# Stupidly using source as a variable name twice
	source_ra   = np.rad2deg(source._ra)
	source_dec  = np.rad2deg(source._dec)
	source_name = source.name
	
	print('Source is: %s'%source.name)
	
	source = tbl.data[0]
	
	source['SOURCE_ID'] = 1
	source['SOURCE']    = source_name
	source['VELDEF']    = 'RADIO'
	source['VELTYP']    = 'GEOCENTR'
	source['FREQID']    = 1
	source['RAEPO']     = source_ra
	source['DECEPO']    = source_dec
	source['EQUINOX']   = 'J2000'
	
	# Things I'm just making up
	source['IFLUX']    = 0
	source['QFLUX']    = 0
	source['UFLUX']    = 0
	source['VFLUX']    = 0
	source['ALPHA']    = 0
	source['FREQOFF']  = 0
	
	tbl.data[0] = source
	
	return tbl

def config_antenna(tbl,array):
	""" Configures the antenna table.
	
	Parameters
	----------
	tbl: pyfits.hdu
	  table to be configured
	"""
	
	antenna = tbl.data
	for i in range(0,tbl.data.size):
	  
		antenna[i]['ANNAME']      = array+'_%i'%i
		antenna[i]['ANTENNA_NO']  = i
		antenna[i]['ARRAY']       = 1
		antenna[i]['FREQID']      = 1
		antenna[i]['NO_LEVELS']   = 12
		antenna[i]['POLTYA']      = 'X'#change for dual
		antenna[i]['POLTYB']      = 'Y'
		antenna[i]['POLAA']       = 0
		antenna[i]['POLAB']       = 0
		  
		  
		tbl.data = antenna
	
	return tbl




def generate_visibitites(uvws, elapsed, julian_midnight, array, params, pointing, freq_params, advanced, src_list, t_len, amper_array):

	"""
	generate visibilities in multidimensional format:
	time, channels, baselines, polarisation, then data=(real, imaginary) 
	
	"""
	print '--------------generating visibilities--------------'
	#create multi-dimensional empty array to hold data

	RAcent = pointing['RAcent']
	DECcent = pointing['DECcent']
	RAc = np.array([int(RAcent[0]),int(RAcent[1]),float(RAcent[2])])
	DECc = np.array([int(DECcent[0]),int(DECcent[1]),float(DECcent[2])])
	
	vis_time = elapsed
	vis_date = np.ones(t_len)*julian_midnight

	t_len=params[0]
	chan_len=params[1]
	pol_len=params[2]
	bl_len=params[3]
	uv = np.zeros((t_len,chan_len,pol_len,bl_len,2))
	
	
	ref_freq = freq_params[0]
	chan_width = freq_params[1]	
	nband = freq_params[2]
	nchan = freq_params[3]

	doPB = advanced['doPB']
	derotation = advanced['derotation']
	derotate = advanced['derotate']

	cell_size = 5.# arcseconds
	if doPB == 'T':
		#generate the Mueller Terms and keep in memory
		#for multi-freq sims ... this needs rethinking
		#The Mueller Terms are dependednt on freq
		fov = 450 * arcmin2rad
		freq = (ref_freq + chan_width )
		xyJones = makeJonesMatrix(fov,freq,cell_size)
		x = xyJones[0]
		y = xyJones[1]
		Jones = xyJones[2:]
		
		#startMueller = time.time()
		#endMueller = time.time()
		#print 'Time taken to create the 4 XX Mueller Terms = %0.4f mins' % ((endMueller - startMueller)/60.)

		MtermsXX = MuellerMatrix(x,y,Jones,'XX')
		MtermsYY = MuellerMatrix(x,y,Jones,'YY')
		MtermsXY = MuellerMatrix(x,y,Jones,'XY')
		MtermsYX = MuellerMatrix(x,y,Jones,'YX')
		Mterms = (MtermsXX,MtermsXY,MtermsYX,MtermsYY)
	elif doPB == 'F':
		Mterms = np.zeros(4)



	srclist = np.loadtxt(src_list,comments='#')
	#vp_all = np.zeros((t_len,pol_len,len(srclist)+1))
	vp_all = np.zeros((t_len,pol_len))

	ha_all = np.zeros(t_len)
	pa_all = np.zeros(t_len)
	for ti in range (int(t_len)):
		ha = timetoHA(vis_time[ti], vis_date[ti], array, RAc[0])
		ha_all[ti] = ha	

		basel_uvw = uvws[ti] * light_speed

		#for each freq channel
		for c in range (int(chan_len)):
			freq = (ref_freq + chan_width * c)*1e6
			freq_scale = freq/light_speed


			for p in range(int(pol_len)): 
				if p == 0:
					pol = 'XX'
				elif p == 1:
					pol = 'YY'
				elif p == 2:
					pol = 'XY'
				elif p == 3:
					pol = 'YX'
				else:
					print 'Polarisation number not recognised.'
					#sys.exit()


				if np.shape(srclist)[0] == 9:
					srclist = np.reshape(srclist,(1,9))

				vis=[]
				i = 0
				for s in srclist:

					visa,vp,pa = sourceVisibility(s,RAc,DECc,pol,Mterms,basel_uvw,ha,pointing['lat'],advanced,freq_scale)
					vis.append(visa)
					i+=1
					vp_all[ti,p] = vp
					pa_all[ti] = pa				
				
				vis=np.sum(vis,axis=0)
				#simulate Gaussian noise
				mu = 0.
				sigma = 1.e-3*np.sqrt(bl_len)*np.sqrt(t_len) 
				#sigma = 0
				points = len(vis)
				noise = np.zeros((points), dtype=np.complex64)
				noise.real = normal(mu, sigma, points)
	        		noise.imag = normal(mu, sigma, points)
				

				
				#add noise to data and finalise the visibilities

				vis_re = vis.real + noise.real 
				vis_im = vis.imag + noise.imag
				for b in range(bl_len):
					uv[ti,c,p,b,0] = vis_re[b]
					uv[ti,c,p,b,1] = vis_im[b]
	#savePBData(vp_all,amper_array,srclist)
	#DR = '/raid/scratch/haydenr/Dropbox/Research/A_Proj/github/DR_Limit_Analytical/'
	#np.savetxt(DR+'PB_FullMueller_HA_HWHM_1mint.txt',ha_all,fmt='%1.4e')
	#np.savetxt(DR+'PB_FullMueller_PA_HWHM_1mint.txt',pa_all,fmt='%1.4e')
	return uv





def savePBData(data, outfile, srclist):

	if len(srclist)==1:
		np.savetxt(outfile,data,fmt='%1.4e')
	else:
		# Write the array to disk
		with file(f, 'w') as outfile:
			# I'm writing a header here just for the sake of readability
			# Any line starting with "#" will be ignored by np.loadtxt
			outfile.write('# Array shape: {0}\n'.format(data.shape))
			
			# Iterating through a ndimensional array produces slices along
			# the last axis. This is equivalent to data[i,:,:] in this case
			for data_slice in data:
		
				# The formatting string indicates that I'm writing out
				# the values in left-justified columns 7 characters in width
				# with 2 decimal places.  
				np.savetxt(outfile, data_slice, fmt='%-7.6f')
				
				# Writing out a break to indicate different slices...
				outfile.write('# New slice\n')

	




def sourceVisibility(s,RAc,DECc,pol,Mterms,basel_uvw,ha,lat,advanced,freq_scale):


	MtermsXX = Mterms[0]
	MtermsXY = Mterms[1]
	MtermsYX = Mterms[2]
	MtermsYY = Mterms[3]

	doPB = advanced['doPB']
	wterm = advanced['wterm']


	x = (np.cos((hmstora(s[0],s[1],s[2])*deg2rad) 
		- (hmstora(RAc[0],RAc[1],RAc[2])*deg2rad))
		* np.cos((dmstodec(s[3],s[4],s[5])*deg2rad)))

	y = (np.sin((hmstora(s[0],s[1],s[2])*deg2rad) 
		- (hmstora(RAc[0],RAc[1],RAc[2])*deg2rad))
		* np.cos((dmstodec(s[3],s[4],s[5])*deg2rad)))

	z = np.sin(dmstodec(s[3],s[4],s[5])*deg2rad)
	xyz = np.array([x,y,z])

	lmn = computeUVW(xyz,0,dmstodec(DECc[0],DECc[1],DECc[2])*deg2rad)
	l = lmn[0]
	m = lmn[1]
	n = lmn[2]
	#"""
	src_pos = np.array([l, m , n - 1]) 
	I = s[6]
	p0 = s[7] * I
	Q = p0*np.cos(s[8]*np.pi/180)
	U = p0*np.sin(s[8]*np.pi/180)
	V = 0
	xxFlux = (I + Q)
	xyFlux = U #U + jV; V = 0
	yxFlux = U #U - jV; V = 0
	yyFlux = (I - Q)
	#print 'src_pos = ', src_pos[0:2]
	#sys.exit()
	if doPB == 'T':
		#rotating pb is introducing artificial errors
		#rotate the source and find the corresponding pb attenuation
		dech = DECc[0]
		PA = parangle(ha, dech, lat)

		orig = (0,0)
		newx,newy =  rotate2d(PA,src_pos[0:2],orig)
		if pol == 'XX':
			Mxx,Mxxy,Mxyx,Mxyxy = findMuellerValue(newx,newy,MtermsXX)
			vp = Mxx*xxFlux - Mxxy*xyFlux - Mxyx*yxFlux + Mxyxy*yyFlux
		elif pol == 'YY':
			Myxyx,Myxy,Myyx,Myy = findMuellerValue(newx,newy,MtermsYY)
			vp = Myxyx*xxFlux - Myxy*xyFlux - Myyx*yxFlux + Myy*yyFlux
		elif pol == 'XY':
			Mxyx, Mxy, Mxyyx, Mxyy = findMuellerValue(newx,newy,MtermsXY) 
			vp = Mxyx*xxFlux + Mxy*xyFlux - Mxyyx*yxFlux - Mxyy* yyFlux
		elif pol == 'YX':
			Myxx, Myxxy, Myx, Myxy = findMuellerValue(newx,newy,MtermsYX)
			vp = Myxx*xxFlux - Myxxy*xyFlux + Myx*yxFlux - Myxy*yyFlux
	
	elif doPB == 'F':
		vp = 1 * xxFlux
		PA = 0
		#print 'vp = ', vp

	if wterm == 'F':
		#print 'not simulating the wterm'
		visa=(vp * np.exp(-2j*np.pi *  np.dot(basel_uvw[:,0:2] 
		* freq_scale, src_pos[0:2]))) 
			
	elif wterm == 'T':
		#print 'simulating the wterm'
		visa=(vp * np.exp(-2j*np.pi *  np.dot(basel_uvw[:,0:3] 
		* freq_scale, src_pos[0:3]))) 						
		
	return visa, vp, PA



def float2str(inarray):
	strs = []
	for n in inarray:
		strs.append(str(n))
	return strs



def basels(configfile):

	ants_xyz=np.genfromtxt(configfile)
	#print ants_xyz
	#sys.exit()
	N = ants_xyz.shape[0]

	bl_pairs = []
	bl_vector = []
	for i in range(N):
    		for j in range(i+1, N):
   			bl_pairs.append((i,j))
			bl_vector.append(ants_xyz[i]-ants_xyz[j])
	#print np.shape(bl_vector)
	#print np.array(bl_vector)[:,0]
	#plt.scatter(np.array(bl_vector)[:,0],np.array(bl_vector)[:,1])
	#plt.show()
	return bl_pairs,bl_vector




def computeUVW(xyz,H,d):
	""" Converts X-Y-Z coordinates into U-V-W
	
	Uses the transform from Thompson Moran Swenson (4.1, pg86)
	
	Parameters
	----------
	xyz: should be a np array [x,y,z]
	H: float (degrees)
	  is the hour angle of the phase reference position
	d: float (degrees)
	  is the declination
	"""
	sin = np.sin
	cos = np.cos
	xyz = np.matrix(xyz) # Cast into a matrix

	
	trans= np.matrix([
	  [sin(H),         cos(H),        0],
	  [-sin(d)*cos(H), sin(d)*sin(H), cos(d)],
	  [cos(d)*cos(H), -cos(d)*sin(H), sin(d)]
	])
	
	uvw = trans * xyz.T
	
	uvw = np.array(uvw)
	
	return uvw[:,0]


########################################################

def config_uv_data(params, tbl_uv_data, antenna_array, source, pointing, trange, freq_params, advanced, array_config, src_list, t_len, amper_array):

	t_len=params[0]
	chan_len=params[1]
	pol_len=params[2]
	bl_len=params[3]
	ri_len=2
	#weights = [1  in range(0,2048)]
	#---------------------1 determine time stamps-----------#
	
	
	#start time of observation in seconds since J2000 epoch
	#obtain from time.time() or  specify a date and time using http://www.epochconverter.com/
	
	timestart = 1458763200 # Wed, 23 Mar 2016 20:00:00 (UTC)
	timestamps = []
	
	for t in range(int(0),int(t_len)):
		timestamp = trange[t]*3600 + timestart # in seconds
	    	timestamps.append(timestamp)
	
	
	firststamp = timestamps[0]
	julian = ephem.julian_date(time.gmtime(firststamp)[:6])
	
	# Ephem returns julian date at NOON, we need at MIDNIGHT
	julian_midnight = int(julian)+1
	
	elapsed = []
	for timestamp in timestamps:
		elapsed.append((ephem.julian_date(time.gmtime(timestamp)[:6]) - julian_midnight))
	
	
	#--------------2. create baseline ids-------------------#
	
	
	
	
	bl_pairs, bl_vectors = basels(array_config)
	baselines = []
	for bl in range(int(0),int(bl_len)):
		# Baseline is in stupid 256*baseline1 + baseline2 format
		ant1, ant2 = bl_pairs[bl][0], bl_pairs[bl][1] 
		bl_id = 256*ant1 + ant2
		bl_vector = bl_vectors[bl]
		baselines.append((bl_id,bl_vector)) 
		#baselines.append((bl_id,-bl_vector)) 
			
	#-------------------3. compute the UVW coordinates-----------------
	
	
	
	print('Computing UVW coordinates...\n')
	# Extract the timestamps and use these to make source our phase centre
	uvws = []
	for timestamp in timestamps:
		t = datetime.datetime.utcfromtimestamp(timestamp)
		antenna_array.update(str(t))
		source.compute(antenna_array)
		for baseline in baselines:
			vector = baseline[1]
			#print vector
			H1,d1 = (antenna_array.sidereal_time() - source.ra,source.dec)
	    		uvws.append(computeUVW(vector,H1,d1)) 
			
	#----------------------4. Compute the visibilities-------------------------------------------	

	# This array has shape t_len, num_ants, 3
	# and units of SECONDS
	uvws = np.array(uvws)
	uvws = uvws.reshape(uvws.size/bl_len/3,bl_len,3) / light_speed
	#generate the visibilities
	uvdata = generate_visibitites(uvws, elapsed, julian_midnight, antenna_array, params, pointing, freq_params, advanced, src_list, t_len, amper_array)

	
	print('\nReformatting numpy vis format -> FITS IDI UV_DATA')
	print('--------------------------------------------')
	 
	# The actual data matrix is stored per row as a multidimensional matrix
	# with the following mandatory axes:
	# COMPLEX     Real, imaginary, (weight)
	# STOKES      Stokes parameter
	# FREQ        Frequency (spectral channel)
	# RA          Right ascension of the phase center
	# DEC         Declination of the phase center 
	flux = np.ndarray(shape=(chan_len,pol_len,ri_len))
	print 'shape(flux)', np.shape(flux)
	print 'shape(uvdata)', np.shape(uvdata)
	 
	# This step takes ages.
	# I imagine there is some way to massage the hdf5 array
	# to do this a lot quicker than iterating over the indexes 
	print('\nCreating multidimensional UV matrix...')
	for t in range(0,int(t_len)):
		print('processing time sample set %i/%i'%(t+1,t_len))
		for bl in range(0,bl_len):
	    
			# Create a 1D index for the uv_data table
			i = t*bl_len + bl
			
			# Swap real and imaginary
			flux[:,:,0] = uvdata[t,:,:,bl,0]
			flux[:,:,1] = uvdata[t,:,:,bl,1]
			
			tbl_uv_data.data[i]['FLUX']     = flux.ravel()
			#tbl_uv_data.data[i]['WEIGHT']   = weights
			
			tbl_uv_data.data[i]['UU']       = uvws[t][bl][0]
			tbl_uv_data.data[i]['VV']       = uvws[t][bl][1]
			tbl_uv_data.data[i]['WW']       = uvws[t][bl][2]
			
			# baselines is a list: [ (id, vec), (id,vec) ... ]
			tbl_uv_data.data[i]['BASELINE'] = baselines[bl][0]
			
			# Date and time
			# Date is julian date at midnight that day
			# The time is seconds since midnight
			tbl_uv_data.data[i]['DATE']     = julian_midnight
			tbl_uv_data.data[i]['TIME']     = elapsed[t]
			
			tbl_uv_data.data[i]['SOURCE']   = 1
			tbl_uv_data.data[i]['FREQID']   = 1
			tbl_uv_data.data[i]['INTTIM']   = 3
	  
	print('\nData reformatting complete')
	return tbl_uv_data  
	print('DONE.')








if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'], version='$Id$')
        parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
	parser.add_option("-f", "--file", dest="filename",help="read data from FILENAME")
        (options, args) = parser.parse_args()
        if options.verbose: print time.asctime()
        #main()

	#Comment out if need an input file
	#change def main() to def main(infile)
	if options.filename:
		print "reading %s..." % options.filename
		main(options.filename)

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


