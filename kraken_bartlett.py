#==================================================================
#  
#  KRAKEN: Elba waveguide
#  Source localization with the Bartlett estimator
#  Faro, Ter 28 Jun 2022 19:37:45 WEST 
#  Written by Orlando Camargo Rodriguez 
#  
#==================================================================

from os import *
import sys
from numpy import *
from scipy.io import *
from matplotlib.pyplot import *
from wkrakenenvfil import *
from readshd import *

print('Elba waveguide:') 

case_title = 'Elba waveguide'

freq =  170.0 # frequency in Hz

zs = array([79.0]) # Source depth
rs = array([ 0.0]) # Auxiliary parameter 

rmax   = 5437.0
rmaxkm = rmax/1000
nra    = 1
rarray = array([rmax]); rarraykm = rarray/1000
first_hyd = 17.7
dhyd      = 2 
last_hyd  = 89.7
zarray = arange(first_hyd,last_hyd+dhyd,dhyd)
nza = zarray.size

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

sspdata = loadtxt("ssp.dat")

z = sspdata[:,0]; zmax = max( z )
c = sspdata[:,1]
nz = z.size
cs = zeros(nz)
rho = ones(nz)
apt = cs
ast = cs
type  = 'H'
itype = 'N'
nmesh = 0 # Number of mesh points (let KRAKEN decide)
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

print( "Generating the observation..." )

wkrakenenvfil('elba',case_title,source_data,surface_data,scatter_data,ssp_data, bottom_data,field_data)

system("krakenc.exe elba")
system("cp field.flp elba.flp")
system("field.exe elba < elba.flp")

filename = 'elba.shd'
pressure,geometry = readshd(filename,nan,nan,nan)

p = squeeze( pressure )
d = p/linalg.norm( p )
R = outer( d,conj(d) )

print( "Generating the replicas..." )

rmin = 4000.0; rminkm = rmin/1000
rmax = 7000.0; rmaxkm = rmax/1000
nra = 201; rarray = linspace(rmin,rmax,nra); rarraykm = rarray/1000

nzs = 124
zs = linspace(1,124,nzs)

source_data = {"zs":zs, "f":freq}
field_data = {"rmax":rmaxkm,"nrr":nra,"rr":rarraykm,"rp":rp,"np":np,"m":m,"rmodes":rmodes,"stype":stype,"thorpe":thorpe,"finder":finder,"rd":zarray,"dr":dr,"nrd":nza}

wkrakenenvfil('elba',case_title,source_data,surface_data,scatter_data,ssp_data, bottom_data,field_data)

system("krakenc.exe elba")
system("cp field.flp elba.flp")
system("field.exe elba < elba.flp")

filename = 'elba.shd'
pressure,geometry = readshd(filename,nan,nan,nan)

p = squeeze( pressure )

bartlett = zeros((nzs,nra))

print('Estimating source range and depth...')

for i in range(nzs):
    pii = squeeze( p[i,:,:] )
    for j in range(nra):
        pj = squeeze( pii[:,j] )
        e = pj/linalg.norm( pj )
        prod = inner( conj(e), d )
        eR = dot(conj(e),R)
        bartlett[i,j] = abs( dot(eR,e) )

maxb = amax( bartlett )
ji = where( bartlett == maxb )
j = ji[0]
i = ji[1]
rsb = rarray[i][0]
zsb = zs[j][0]
thetitle = 'Source located at (' + str(rsb) + ',' + str(zsb) + ') m'

figure()
imshow(bartlett,cmap='jet',extent=[rarraykm[0],rarraykm[-1],zs[-1],zs[0]],aspect='auto')
colorbar()
xlabel('Range (km)',fontsize=18)
ylabel('Depth (m)',fontsize=18)
title(thetitle)

show()

print("done.")
