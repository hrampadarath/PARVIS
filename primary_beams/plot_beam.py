import scipy.io
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import sys
 
 
 
deg2rad = np.pi / 180
 
 
#load the primary beams and convert from matlab format to numpy
beam1 = 'Data1.mat'
mat = scipy.io.loadmat(beam1)
 
 
freq = mat['freq'][0]
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
Jpv = mat['Jpv']
Jph = mat['Jph']
Jqv = mat['Jqv']
Jqh = mat['Jqh']

J1 = np.abs(Jpv[0::1,1:,0])
J1s = np.fliplr(J1)
Jp = np.concatenate((J1s,J1),axis=1)
 
 
###Since the beams are in spherical coordinates, to plot as a 2D image in python
###I had to do an interpolation.
 
#create a 2D grid to interpolate the beams to
x1 = np.linspace(-0.13,0.13)
y1 = np.linspace(-0.13,0.13)
x,y = np.meshgrid(x1,y1)
 
#determine the phi and theta values of the grid
pph = np.arctan2(x,y)
pth = (np.arccos( np.cos(x) * np.cos(y) ))
 
#interpolate the beams over a rectangular mesh.
interpolator = interpolate.RectBivariateSpline(th,ph,Jp)
pJ = interpolator.ev(pth,pph)
 
#log and plot
pJ_log = np.log10(pJ)
plt.contourf(x,y,pJ_log,50)
plt.axis('scaled')
 
plt.ylabel('Y (radians)')
plt.xlabel('X (radians)')
plt.colorbar()
 
#plt.savefig('Jpv_660MHz.jpg')
 
plt.show()
