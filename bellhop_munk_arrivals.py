#=======================================================================
# 
# Bellhop: Munk profile + arrivals
# Faro, Qua Mai 26 19:08:18 WEST 2021
# Written by Orlando Camargo Rodriguez 
# 
#=======================================================================

from os import system
from numpy import *
from scipy.io import *
from matplotlib.pyplot import *
from wbellhopenvfil import *
from read_arrivals_asc import *
from mpl_toolkits.mplot3d import Axes3D

case_title = "Munk profile"

freq   = 50.0
Rmaxkm = 15.0; Rmax = Rmaxkm*1000.0 
Dmax   = 5000.0 
cb   = 1600.0 # sound speed in lower halfspace
rhob =  1.8   # density in lower halfspace

source_nrays    = 1000 # number of propagation rays considered #
source_aperture = -80.0 # maximum launching angle (degrees) #
source_ray_step =  50.0 # ray calculation step (meters) #

#==================================================================
#  
#  Source properties 
#  
#==================================================================

nzs  = 1
zs   = array([4700.0])
rs   = array([   0.0])
zbox = Dmax + 1.0 
rbox = Rmaxkm + 1.0  # km!!!!!
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

nza = 1001
z = linspace(0,Dmax,nza)

c1 = 1500.0
z1 = 1300.0
B = 1.3e3; BxB = B*B
epsilon = 7.37e-3
eta = 2.0*( z - z1 )/B
c   = 1.0 + epsilon*( eta + exp( -eta ) - 1.0 )
c   = c1*c

r = []

ssp_data = {"r":r,"z":z,"c":c}

#==================================================================
#  
#  Bottom:
#  
#==================================================================

rbty   = array([rs[0]-2,Rmax+2]); rbtykm = rbty/1000.0
zbty   = array([Dmax,Dmax])
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
options2    = '''A '''  # No bty file expected  
options3    = '''A'''
options4    = []

nra = 21
rarray = linspace(0,Rmax,nra); rarraykm = rarray/1000.0
zarray = zs

options = {"options1":options1,"options2":options2,"options3":options3,"options4":options4,"rarray":rarraykm,"zarray":zarray}

print("Writing environmental file...")

wbellhopenvfil('munk',case_title,source_data,surface_data,ssp_data,bottom_data,options)

print( "Running Bellhop..." )

system("bellhop.exe munk")

print( "Reading output data..." )

Arr, Pos = read_arrivals_asc( 'munk.arr' )

Narr   = int_( squeeze( Arr['Narr'] ) )
delay  = real( squeeze( Arr['delay'] ) )
A      = abs( squeeze( Arr['A'] ) )
rarray = squeeze( Pos['receiver_ranges'] )

Nrr = rarray.size

fig = figure()
ax = fig.gca(projection='3d')

for i in range(1,Nrr):
    tau = delay[i,0:Narr[i]]
    amp = abs( A[i,0:Narr[i]] )
    rangei = rarray[i]*ones(Narr[i])
    for j in range(Narr[i]):
        ax.plot([tau[j],tau[j]], [rangei[j],rangei[j]], [0.0,amp[j]],'b')    
        ax.plot([tau[j],tau[j]], [rangei[j],rangei[j]], [0.0,amp[j]],'bo')
xlabel('Travel time (s)')
ylabel('Range (in m)')
title('Bellhop - Munk profile')
show()

print('done.')
