import xml.etree.ElementTree as ET
import os,sys
import datetime
import re
import time

"""
pyfitsidi documentation - http://fits.gsfc.nasa.gov/registry/fitsidi.html
"""

tools = './tools/'

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






if len(sys.argv)==2:

	inFile = sys.argv[1]
	print 'Using inputs specified in', inFile
	if os.path.isfile(inFile) :
		control = parse_inp(inFile)
	else :
		print "Error:" + inFile + "does not exist, quitting."
		sys.exit()
else:
	print "Error: no input parameter file included ... quitting"
	print 'Usage: "python code.py inputfile.inputs"'
	sys.exit()

inFile = sys.argv[1]
print ' Creating XML file to be used in pyfitsidi.py'

array 		= control['array'][0]

outfolder	= control['outfolder'][0]
templatexml 	= tools+control['templatexml'][0]
outxml 		= outfolder+control['configxml'][0]

npol  		= control['pol_len'][0]
nband 		= control['nband'][0]
nchan 		= control['chan_len'][0]


chan_width 	= control['chan_width'][0] #MHz
ref_freq  	= control['start_freq'][0]# for now assume this is the start freq in the 1st channel

notes 		= control['notes'][0]

#today = datetime.date.today()
#t = today.timetuple()
timeobs = 1458756000 + (2*3600)# Wed, 23 Mar 2016 20:00:00 (UTC)
t = time.gmtime(timeobs)



os.system('cp '+templatexml+' '+outxml)

tree = ET.parse(outxml)

root = tree.getroot()

root.set('name','DVA1/JVLA-A')

#edit the notes 

#root[0].text = " 'Example configuration file made with pyxml' "
root[0].text = " '"+ notes + "' "

# edit the PARAMETERS data
"""
    A number of arrays have dimensions that depend on the parameters of the data set
    or of the table to which they belong. The notations used for these parameters are
    listed below (and in the FITS=IDI documentation)

    The values here are mainly used to help set up the table columns in the HDU.
    For example, if you had 8 bands, with 32 channels, then the UV_DATA table needs
    to know about this, and it should assign an array of size 8*32*2(real & imag)*32bits.

    Data Set Parameters
    nstokes The number of Stokes parameters in the data set
    nband   The number of bands in the data set
    nchan   The number of channels in the data set
    ntone   The maximum number of pulse-cal tones in a PHASE-CAL table
    norb    The number of orbital parameters in an ANTENNA_GEOMETRY table
    npoly   The number of terms in a delay polynomial in an INTERFEROMETER_MODEL table
    ntab    The maximum number of tabulated vals/terms for in a GAIN_CURVE table

    <PARAMETERS>
      <NSTOKES>2     </NSTOKES> 
      <NBAND>  4     </NBAND>  
      <NCHAN>  32  </NCHAN> 
      <NTONE>  1     </NTONE>
      <NORB>   1     </NORB>
      <NPOLY>  1     </NPOLY>
      <NTAB>   1     </NTAB>
      <NPCAL>  1     </NPCAL>
    </PARAMETERS>

Mostly would only have to change the first 3 params

"""

root[1][0].text = ' '+str(npol)+' ' # NSTOKES
root[1][1].text = ' '+str(nband)+' ' #NBANDS
root[1][2].text = ' '+str(nchan)+' ' 


#edit the common header values

"""
    common header values
    ====================
    
    These values are shared across all tables and are mandatory.
    THESE SHOULD BE THE SAME FOR EACH TABLE AND SHOULD ONLY BE SET HERE!

    TABREV      Revision number of the table definition (normally 1)
    NO_STKD     The number of Stokes parameters
    STK_1       The first Stokes parameter
    NO_BAND     The number of bands
    NO_CHAN     The number of spectral channels
    REF_FREQ    The file reference frequency in Hz
    CHAN_BW     Channel BW in Hz for the 1st band in freq. group with freq. ID 1
    REF_PIXL    The reference pixel for the frequency axis
    RDATE       Reference date: the date for which the time params in ARRAY_GEOMETRY apply

    Each table also has a EXTNAME, but this is set up automatically so you don't
    have to worry about it.
    
    Numeric Codes for Stokes Parameters:
    1:I, 2:Q, 3:U, 4:V, -1:RR, -2:LL, -3:RL, -4:LR, -5:XX, -6:YY, -7:XY, -8:YX

    <COMMON>
       <OBSCODE>   'MEDI0001'      </OBSCODE>           
       <RDATE>     '2011-04-25'    </RDATE>         
       <NO_STKD>   params['NSTOKES']</NO_STKD>
       <STK_1>     -1              </STK_1>
       <NO_BAND>   params['NBAND'] </NO_BAND>
       <NO_CHAN>   params['NCHAN'] </NO_CHAN>
       <REF_FREQ>  398.0E+06       </REF_FREQ>
       <CHAN_BW>   0.01953125E+06  </CHAN_BW>
       <REF_PIXL>  1.00E+00        </REF_PIXL>
       <TABREV>    1               </TABREV>
       <EXTVER>    1               </EXTVER>
    </COMMON



"""

root[3][0].text = " 'pySimulator' "
root[3][1].text = " '"+str(t[0])+'-'+str(t[1])+'-'+str(t[2])+"' "
root[3][2].text = ' '+str(npol)+' ' # NSTOKES
root[3][3].text = ' -5 '
root[3][6].text = ' '+str(ref_freq)+'E+06'+' '
root[3][7].text = ' '+str(chan_width)+'E+06'+ ' '





"""
 ARRAY_GEOMETRY table header
   ===========================
   
   The ARRAY_GEOMETRY tables define the arrays used in the file. Each ARRAY_GEOMETRY table 
   lists the antennas that are part of that array together with their coordinates. It also
   provides information about the time system used for that array.
   
   EXTVER  Array number
   ARRNAM  Array name
   FRAME   Coordinate frame
   ARRAYX  x coordinate of array center in m (important, check the convention!!)
   ARRAYY  y coordinate of array center in m
   ARRAYZ  z coordinate of array center in m
   NUMORB  Number of sattelites. Normally safe to assume this is zero.
   FREQ    Reference frequency
   TIMSYS  Time system
   RDATE   Reference date
   GSTIA0  Grenwich Sidereal time at 0 hrs
   DEGPDY  The Earth's rotation rate in deg per day

   UT1UTC  UT1 - UTC
   IATUTC  IAT - UTC
   POLARX  x coordinate of north pole
   POLARY  y coordinate of north pole
   GSTIA0, DEGPDY, UT1UTC notes: The default values for the time system 
   are taken from http://fits.gsfc.nasa.gov/registry/fitsidi.html
   I haven't checked these are actually correct (on my todo list)!
   
   ARRAYX, ARRAYY, ARRAYZ notes: The values below are for Medicina.
   These are VERY important and you'll have to change them. This might help:
   http://www.oc.nps.edu/oc2902w/coord/llhxyz.htm   
   
   -->

   <ARRAY_GEOMETRY>
     <EXTNAME>  'ARRAY_GEOMETRY' </EXTNAME>
     <ARRAYX>   4461.122E+03   </ARRAYX>
     <ARRAYY>   919.46E+03     </ARRAYY>
     <ARRAYZ>   4449.776E+03   </ARRAYZ>
     <ARRNAM>   'MEDICINA'     </ARRNAM>
     <NUMORB>   0              </NUMORB>
     <FREQ>     408.0E+06      </FREQ> 
     <FRAME>    'GEOCENTRIC'   </FRAME>
     <TIMSYS>   'UTC'          </TIMSYS>
     <TIMESYS>  'UTC'          </TIMESYS>
     <GSTIA0>   3.30909596261338038E+02 </GSTIA0>
     <DEGPDY>   3.60985644973299998E+02 </DEGPDY>
     <POLARX>   2.08099999999999996E-01 </POLARX>
     <POLARY>   2.80019999999999989E-01 </POLARY>
     <UT1UTC>  -1.63126999999999995E-01 </UT1UTC>
     <IATUTC>   3.30000000000000000E+01 </IATUTC>
   </ARRAY_GEOMETRY>


"""



root[5][1].text = ' -1601.389E+03 '
root[5][2].text = ' -5042.634E+03 '
root[5][3].text = ' 3553.876E+03 '
root[5][4].text = " '"+array+ "' "
root[5][6].text = ' '+str(ref_freq)+'E+06'+' '
           

"""
UV_DATA table header
===================

A UV_DATA table contains a set of visibility data matrices. If there is more than
one UV_DATA table in the file then no two tables shall contain data for overlapping
times and the tables shall appear in time order in the file1.

This one is the biggest, most complicated and most important table. Spend a little time
making sure you've got this right.

TABREV  2: This should OVERRIDE the common value that we set above as 1
NMATRIX 1: Don't think we've got a choice here

The UV_DATA is a multidimensional array (6 levels in general)
MAXIS   M = number axes in regular matrix
MAXISm I Number pixels on axis m = 1 to M
CTYPEm A Name of regular axis m = 1 to M
CDELTm E Coordinate increment on axis m = 1 to M
CRPIXm E Reference pixel on axis m = 1 to M
CRVALm E Coordinate value at reference pixel on axis m = 1 to M

Notes about the axes, with Medicina for an example:
* first axis is complex, real imag
* second axis is stokes, we only have 1
* third axis is number of frequency chans (1024 for us)
* this axis is the number of discrete bands (1 for medicina)
* RA  - not really sure why they decided to make this an axis
* DEC - this is part of the charm of FITS IDI

Which column is the visibility matrix in? I've hard coded this to 11
TMATXn L T - column n contains the visibility matrix

Finally, some values that you might to change: 
EQUINOX   Mean equinox (probably J2000)
WEIGHTYP  Type of data weights
DATE-OBS  Observing date
TELESCOP  Telescope name
#OBSERVER  Observer's name
VIS SCAL  Visibility scale factor
SORT      Sort order, * does no sorting (thus the quickest?)

-->

<UV_DATA>
  <EXTNAME>   'UV_DATA'   </EXTNAME>
  <TABREV>    2           </TABREV>
  <DATE-OBS>  '2011-04-25'</DATE-OBS>
  <TELESCOP>  'BEST-2'    </TELESCOP>
  <OBSERVER>  'DPRICE'    </OBSERVER>
  <EQUINOX>   'J2000'     </EQUINOX>
  <WEIGHTY>   'NORMAL'    </WEIGHTY>
  <SORT>      '*'         </SORT>
  <NMATRIX>   1           </NMATRIX>
  <MAXIS>     6           </MAXIS>
  <MAXIS1>    2           </MAXIS1>
  <CTYPE1>    'COMPLEX'   </CTYPE1>
  <CDELT1>    1.000E+00   </CDELT1>
  <CRPIX1>    1.000E+00   </CRPIX1>
  <CRVAL1>    1.000E+00   </CRVAL1>
  <MAXIS2>    params['NSTOKES'] </MAXIS2>
  <CTYPE2>    'STOKES'    </CTYPE2> 
  <CDELT2>    -1.000E+00  </CDELT2>
  <CRPIX2>    1.0000E+00  </CRPIX2>
  <CRVAL2>    -1.000E+00  </CRVAL2>
  <MAXIS3>    params['NCHAN'] </MAXIS3>
  <CTYPE3>    'FREQ'      </CTYPE3>
  <CDELT3>    0.01953125E+06</CDELT3>
  <CRPIX3>    1.00000E+00 </CRPIX3>
  <CRVAL3>    398.000E+06 </CRVAL3>
  <MAXIS4>    params['NBAND'] </MAXIS4>
  <CTYPE4>    'BAND'      </CTYPE4>
  <CDELT4>    1.000E+00   </CDELT4>
  <CRPIX4>    1.000E+00   </CRPIX4>
  <CRVAL4>    1.000E+00   </CRVAL4>
  <MAXIS5>    1           </MAXIS5>
  <CTYPE5>    'RA'        </CTYPE5>
  <CDELT5>    0.000E+00   </CDELT5>
  <CRPIX5>    1.000E+00   </CRPIX5>
  <CRVAL5>    0.000E+00   </CRVAL5>
  <MAXIS6>    1           </MAXIS6>
  <CTYPE6>    'DEC'       </CTYPE6>
  <CDELT6>    0.000E+00   </CDELT6>
  <CRPIX6>    1.000E+00   </CRPIX6>
  <CRVAL6>    0.000E+00   </CRVAL6>
  <TMATX11>   T           </TMATX11>
</UV_DATA>

"""

root[8][2].text = " '"+str(t[0])+'-'+str(t[1])+'-'+str(t[2])+"' "
root[8][3].text = " '"+array+ "' "
root[8][4].text = " 'HRAMPADARATH' "
#root[8][17].text = " -5.00E+00 " #firsr pol XX
root[8][19].text = " -5.00E+00 " #firsr pol XX
root[8][22].text = ' '+str(chan_width)+'E+06'+ ' '
root[8][24].text = ' '+str(ref_freq)+'E+06'+' '



tree.write(outxml)

print 'XML file configured and saved to:', outxml

 
