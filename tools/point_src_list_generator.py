import numpy as np
from astronomyv1 import deltasep, alphasep, hmstora, dmstodec, newpos
import sys,os



srcfile = 'single_nopol_src'

if os.path.isfile(srcfile) == True:
	print 'file exists ... removing'
	os.system('rm -rf '+ srcfile)


print 'Generating a list of random sources'

f = open(srcfile, 'ab')
a = "#list of sources to be used in the A-projection simulation"
b = "#centre:     ra='05:00:00', dec='45:00:00'"
c = "#RAh RAm RAs DECd DECm DECs Stokes_I (Jy) frac_pol pol_angle " 
f.write(a)
f.write('\n')
f.write(b)
f.write('\n')
f.write(c)
f.write('\n')


dech = 45


RAc = np.array([5.,0.,0.])
DECc = np.array([dech,0.,0.])

#3c286 pol properties
#frac_pol = 0.086 
#pol_angle = 33 
frac_pol = 0
pol_angle = 0

c = 299792458.   # speed of light
nu = 700.e6      # observing frequency
l = (np.degrees(1.02 * ((c/nu) / 15))*60.)/4. #radial position of source from pointing centre. Here we define the source to be at the FWHM

dDEC = np.array([l])
#dDEC = [-1*l,l] #np.arange(-90,90,5) # random Declination offset in arcmins
#dDEC = np.arange(-l,l+5,5)
dRA = np.zeros(len(dDEC)) #np.arange(-90,90,5) # random RA offset in arcmins

random_flux = 1#np.random.random_integers(1,21)   # random source flux between 2 to 195 mJy


for i in range(len(dDEC)):
	
	sRA,sDEC = newpos(dRA[i]*60,dDEC[i]*60,RAc[0],RAc[1],RAc[2],DECc[0],DECc[1],DECc[2])
	src = str(sRA[0])+' '+str(sRA[1])+' '+str(sRA[2])+' '+str(sDEC[1])+' '+str(sDEC[2])+' '+str(sDEC[3])+' '+str(random_flux)+' '+str(frac_pol)+' '+str(pol_angle)
	f.write(src)
	f.write('\n')


print 'Completed'
