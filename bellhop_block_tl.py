#=======================================================================
# 
# Bellhop: Block problem
# Faro, Seg 11 Abr 2022 12:45:18 WEST 
# Written by Orlando Camargo Rodriguez 
# 
#=======================================================================

from os import system
from numpy import *
from scipy.io import *
from matplotlib.pyplot import *
from wbellhopenvfil import *
from readshd import *

case_title = "Block problem"

freq   = 100.0
Rmaxkm =   5.0; Rmax = Rmaxkm*1000.0 
Dmax   = 1000.0 
cw   = 1500.0 # sound speed in water
cb   = 1700.0 # sound speed in lower halfspace
rhow =  1.0   # density in water
rhob =  2.0   # density in lower halfspace

source_nrays    =  301 # number of propagation rays considered #
source_aperture = 20.0 # maximum launching angle (degrees) #
source_ray_step =  5.0 # ray calculation step (meters) #

#==================================================================
#  
#  Source properties 
#  
#==================================================================

nzs  = 1
zs   = array([500.0])
rs   = array([  0.0])
zbox = 1001.0 
rbox = 5.1  # km!!!!!
box  = array([source_ray_step,zbox,rbox])
thetas = array([source_nrays,-source_aperture,source_aperture]) 
p      = zeros(1)
comp   = ''

source_data = {"zs":zs, "box":box, "f":freq, "thetas":thetas, "p":p, "comp":comp}

#==================================================================
#  
#  Surface definition:
#  
#==================================================================

itype = ''
xati  = []  # The *.ati file won't be written  
p     = []  # Surface properties
aunits= ''

surface_data = {"itype":itype,"x":xati,"p":p,"units":aunits}

#==================================================================
#  
#  Sound speed:
#  
#==================================================================

z = array([0.0,Dmax])
c = array([cw,cw])
r = []

ssp_data = {"r":r,"z":z,"c":c}

#==================================================================
#  
#  Bottom:
#  
#==================================================================

rbty   = array([rs[0]-2,2000,2010,2990,3000,Rmax+2]); rbtykm = rbty/1000.0
zbty   = array([Dmax,Dmax,500,500,Dmax,Dmax])
itype  = '''L''' # RID properties
bunits = '''W'''
xbty   = array([rbtykm,zbty]) # Bottom coordinates
p      = array([2000.0,0.0,2.0,0.5,0.0]) # Bottom properties

bottom_data = {"itype":itype,"x":xbty,"p":p,"units":bunits}

#==================================================================
#  
#  Array: 
#  
#==================================================================

options1    = '''CVW''' # No ati file expected  
options2    = '''A*'''
options3    = '''C''' # Coherent TL (rewrite final block accordingly) 
options4    = []

rarray = linspace(0,Rmax,501); rarraykm = rarray/1000.0 
zarray = linspace(0,Dmax,101)

options = {"options1":options1,"options2":options2,"options3":options3,"options4":options4,"rarray":rarraykm,"zarray":zarray}

print("Writing environmental file...")

wbellhopenvfil('block',case_title,source_data,surface_data,ssp_data,bottom_data,options)

print( "Running Bellhop..." )

system("bellhop.exe block")

print( "Reading output data..." )

filename = 'block.shd'
xs = nan
ys = nan
pressure,geometry = readshd(filename,xs,ys)

p = squeeze( pressure, axis=(0,1) )
p = where( p == 0, nan, p )
tl = -20*log10( abs( p ) )

figure(1)
pcolormesh(rarray,zarray,tl,vmin=60,vmax=90,cmap='jet_r',shading='auto')
fill_between(rbty,zbty,Dmax)
colorbar()
plot(rs[0],zs[0],marker="<",markersize=16,color="k")
xlim(0,Rmax)
ylim(Dmax,0)
xlabel('Range (m)')
ylabel('Depth (m)')
title('Bellhop - Block problem')

show()

print("done.")
