#==================================================================
#  
#  SCOOTER: Elba waveguide
#  Gambelas, qui 20 jun 2024 13:18:21 
#  Written by Orlando Camargo Rodriguez 
#  
#==================================================================

import os
from numpy import *
from scipy.io import *
from matplotlib.pyplot import *
from fieldsco import *
from readshd import *
from scalecol import *
from taper import *
from wkrakenenvfil import *

print('Elba waveguide:') 

case_title = 'Elba waveguide'

freq =  100.0 # frequency in Hz
Dmax =  150.0 # depth in meters
cw   = 1500.0 # sound speed in water
cb   = 1700.0 # sound speed in lower halfspace
rhow =    1.0 # density in water
rhob =    2.0 # density in lower halfspace

zs = array([10.0]) # Source depth
rs = array([ 0.0]) # Auxiliary parameter 

rmax = 1000.0; rmaxkm = rmax/1000
nra = 201; rarray = linspace(0,rmax,nra); rarraykm = rarray/1000
zarray = arange(0,Dmax+1,1); nza = zarray.size

layer_info = array([[98.0,1600.0,130.0,1.49,0.1,1.70],
                   [103.0,1800.0,500.0,1.90,0.90,2.50],
		   [128.0,2500.0,900.0,2.40,0.01,0.01]])

#==================================================================
#  
#  Source data
#  
#==================================================================

source_data = {"zs":zs, "f":freq}

#==================================================================
#  
#  Surface data
#  
#==================================================================

bc = 'V'
properties = [] # Not required (vacuum over surface)
reflection = [] # Not required for this case

surface_data = {"bc":bc,"properties":properties,"reflection":reflection}

#==================================================================
#  
#  Scatter data
#  
#==================================================================

# Scatter is not required for this case: 
bumden = [] # Bump density in ridges/km 
eta    = [] # Principal radius 1 of bump 
xi     = [] # Principal radius 2 of bump 

scatter_data= {"bumden":bumden,"eta":eta,"xi":xi}

#==================================================================
#  
#  Sound speed data
#  
#==================================================================

sspdata = loadtxt("elba.ssp")

z = sspdata[:,0]; zmax = max( z )
c = sspdata[:,1]
nz = z.size
cs = zeros(nz)
rho = ones(nz)
apt = cs
ast = cs
type  = 'H'
itype = 'N'
# Number of mesh points to use initially, should be about 10 per vertical wavelenght:
nmesh = 0
sigma = 0.0 # RMS roughness at the surface 
clow  = 1500.0
chigh = 1800.0

cdata = array([z,c,cs,rho,apt,ast])

ssp_data = {"cdata":cdata,"type":type,"itype":itype,"nmesh":nmesh,"sigma":sigma,"clow":clow,"chigh":chigh,"zbottom":zmax}

#==================================================================
#  
#  Bottom data
#  
#==================================================================

layerp     = array([[nmesh, 0.0, layer_info[1,0]],
                    [nmesh, 0.0, layer_info[2,0]]])
layert     = 'HH'
properties = array( layer_info[2,] )

m1 = array( [layer_info[0,],layer_info[0,]] )
m2 = array( [layer_info[1,],layer_info[1,]] )

bdata       = array( [m1,m2] )
bdata[0,1,0] = layer_info[1,0]
bdata[1,1,0] = layer_info[2,0]

units      = 'W';
bc	   = 'A';
sigma      = 0.0 # Interfacial roughness

bottom_data = {"n":3,"layerp":layerp,"layert":layert,"properties":properties,"bdata":bdata,"units":units,"bc":bc,"sigma":sigma}

#==================================================================
#  
#  Field data
#  
#==================================================================

rp    =   0 
np    =   1
m     = 999
rmodes = 'A'
stype  = 'R'
thorpe = 'T'
finder = ' '
dr     = zeros( nza )

field_data = {"rmax":rmaxkm,"nrr":nra,"rr":rarraykm,"rp":rp,"np":np,"m":m,"rmodes":rmodes,"stype":stype,"thorpe":thorpe,"finder":finder,"rd":zarray,"dr":dr,"nrd":nza}

print("Writing environmental file...")

wkrakenenvfil('elbaS',case_title,source_data,surface_data,scatter_data,ssp_data, bottom_data,field_data)

print( "Writing the flp file..." )

flpfil = 'elbaS.flp'
arrayinfo = str( rarraykm[0] ) + ' ' + str( rarraykm[-1] ) + ' /'

fid = open(flpfil, 'w')
fid.write("'RP'");fid.write('\n')
fid.write(str(nra));fid.write('\n')
fid.write(arrayinfo);fid.write('\n')
fid.close()

print( "Running SCOOTER..." )

os.system("scooter.exe elbaS")

print( "Reading output data..." )

fieldsco("elbaS.grn", zs, zarray, array([freq]))

data = loadmat("elbaS.shd.mat")

pressure = data["pressure"]

p = squeeze( pressure, axis=(0,1) )
p = where( abs( p ) < 1e-6, 1e-6, p )

tl = -20*log10( abs( p ) )

figure(1)
imshow(tl,extent=[0,rmax,Dmax,0], aspect='auto',cmap='jet_r',vmin=30,vmax=80)
colorbar()
plot([0,rmax],[98  , 98],'k',linewidth=2)
plot([0,rmax],[103, 103],'k',linewidth=2)
plot([0,rmax],[128, 128],'k',linewidth=2)
plot(rs,zs,marker="<",markersize=16,color="k")
xlabel('Range (m)')
ylabel('Depth (m)')
title('SCOOTER - Elba waveguide')
ylim(Dmax,0)

show()

print("done.")
