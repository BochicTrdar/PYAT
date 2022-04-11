#======================================================================
# 
# Bellhop3D: tl test
# Faro, Qua Mai 26 19:08:53 WEST 2021
# Written by Orlando Camargo Rodriguez 
#
#======================================================================

from os import *
from numpy import *
from scipy.io import *
from matplotlib.pyplot import *
from wbellhop3denv import *
from readshd import *

case_title = '''Wedge problem''' 

freq = 122.0
ray_step = 1.0
slope = 4.55; k = tan( slope*pi/180 )
z0 = 44.4
x1 = -z0/k 
x2 = -x1 
z1 = k*x1 + z0 
z2 = k*x2 + z0
xmax = max([x1,x2])
zmax = max([z1,z2])

#==================================================================
#  
#  Source:
#  
#==================================================================

xs  = array([0.0]) 
ys  = array([0.0]) 
zs  = array([8.3])
box = array([ray_step, xmax/1000.0, 8.0, zmax])  # [step(m) xbox(km) ybox(km) zbox(m)]
thetas = array([401,6,-24.0,24.0]) 
phi    = array([501,46,90.0,110.0])
bearings = array([1,90.0,0.0]) 
p    = [] 
comp = [] 

source_data = {"xs":xs, "ys":ys, "zs":zs, "box":box, "f":freq, "thetas":thetas, "phi":phi, "bearings":bearings, "p":p, "comp":comp}

#==================================================================
#  
#  Surface:
#  
#==================================================================

aunits = ''' ''' 
itype = ''' '''
itype = [] 
xati  = zeros(1)
yati  = zeros(1)
zati  = zeros(1)
p     = zeros(1)

surface_data = {"itype":itype,"xati":xati,"yati":yati,"zati":zati,"p":p,"units":aunits}
 
#==================================================================
#  
#  Sound speed profile:
#  
#==================================================================

c0 = 1488.2

z = array([0.0,zmax])
c = array([c0,c0])
r = zeros(1)

ssp_data = {"r":r,"z":z,"c":c}

#==================================================================
#  
#  Bottom:
#  
#==================================================================
 
bunits = '''W''' # Bottom absorption units: (dB/m)*kHZ 
itype = '''R'''  # Not used, but required...
xbty = array([x1,x2])/1000.0 
ybty = array([-8.0,8.0]) 
zbty = array([[z1,z1],[z2,z2]])
p    = array([1700.0,0.0,1.99,0.5])

bottom_data = {"itype":itype,"xbty":xbty,"ybty":ybty,"zbty":zbty,"p":p,"units":bunits}

#==================================================================
#  
#  Array:
#  
#==================================================================

yarraykm = arange(0.035,5.000,0.005); nya = yarraykm.size
zarray = array([9.3])

options1 = '''SVWT''' # No ati file expected
options2 = '''A*''' # bty file expected
options3 = '''C    3''' # Coherent pressure + 3D run
options4 = ' '

nrd      = 1 
nr       = nya
rd       = zarray
r        = array([min(yarraykm), max(yarraykm)])

options = {"options1":options1,"options2":options2,"options3":options3,"options4":options4,\
           "nrd":nrd,"nr":nr,"rd":rd,"r":r}

print("Writing environmental file...")

wbellhop3denv('wedge',case_title,source_data,surface_data,ssp_data,bottom_data,options)

print( "Running Bellhop3D..." )

system("bellhop3d.exe wedge")

print( "Reading output data..." )

filename = 'wedge.shd'
pressure,geometry = readshd(filename,xs,ys)

p = squeeze( pressure, axis=(0,1) )
p = squeeze( p )

tl = 20*log10( abs( p ) )

thetitle = 'Bellhop3D @ ' + str(freq) + ' Hz'

figure(1)
grid(True)
plot( yarraykm, tl, linewidth=2 )
xlabel('ACROSS SLOPE RANGE [in km]')
ylabel('TLOSS [in dB]')
title( thetitle )

show()

print('done.')
