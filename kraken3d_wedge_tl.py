#======================================================================
# 
# KRAKEN3D: wedge problem
# Faro, Seg 11 Abr 2022 13:29:31 WEST 
# Written by Orlando Camargo Rodriguez 
#
#======================================================================

import os
from numpy import *
from scipy.io import *
from matplotlib.pyplot import *
from wkrakenenvfil import *
from readshd import *

print('Analytic Solution:') 

case_title = '''Analytic Solution''';

refdata = loadtxt("3D_PE_of_F_Sturm_with_cross_term.dat")

yarray = refdata[:,0]; yarraykm = yarray/1000.0
tlref  = refdata[:,1]

rmax = max( yarray ); rmaxkm = rmax/1000.0

#==================================================================
#  
#  Define source data:
#  
#==================================================================

slope = 0.5
k = tan( slope*pi/180 ); z0 = 90

cw = 1500
cp = 2000

rhob = 2; ab = 0.5

#==================================================================
#  
#  Define wedge bathymetry:
#  
#==================================================================

nx = 35; ny = 35

xmin = -rmax; xmax = rmax;
ymin = -rmax; ymax = rmax;

x = linspace(xmin,xmax,nx)

zslope = k*x + z0; zmax = max( zslope )

zbty = zeros((ny,nx))

for i in range(ny):
    zbty[i,:] = zslope

y = linspace(ymin,ymax,ny)

#==================================================================
#  
#  Source data
#  
#==================================================================

freq = 50
xs  = array([ 0.0]) 
ys  = array([ 0.0]) 
zs  = array([10.0])

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

z = array([0,zmax])
c = array([cw,cw])
cs = array([0,0])
rho = array([1.0,1.0])
apt = cs
ast = cs
mtype  = 'H'
itype = 'N'
# Number of mesh points to use initially, should be about 10 per vertical wavelenght:
nmesh = 1001
sigma = 0.0 # RMS roughness at the surface 
clow  =    0.0
chigh = 2000.0

cdata = array([z,c,cs,rho,apt,ast])

ssp_data = {"cdata":cdata,"type":mtype,"itype":itype,"nmesh":nmesh,"sigma":sigma,"clow":clow,"chigh":chigh,"zbottom":zmax}

#==================================================================
#  
#  Bottom data
#  
#==================================================================

layerp     = array([0, 0, zmax])
layert     = 'R'
properties = array( [zmax, cp, 0.0, rhob, ab, 0.0] )
bdata      = []
units      = 'W';
bc	   = 'A';
sigma      = 0.0 # Interfacial roughness

bottom_data = {"n":1,"layerp":layerp,"layert":layert,"properties":properties,"bdata":bdata,"units":units,"bc":bc,"sigma":sigma}

#==================================================================
#  
#  Field data
#  
#==================================================================

zarray = zs
nra   = yarray.size
nza   =   1
rp    =   0 
np    =   1
m     = 999
rmodes = 'A'
stype  = 'R'
thorpe = 'T'
finder = ' '
dr     = zeros( nza )

field_data = {"rmax":rmaxkm,"nrr":nra,"rr":yarraykm,"rp":rp,"np":np,"m":m,"rmodes":rmodes,"stype":stype,"thorpe":thorpe,"finder":finder,"rd":zarray,"dr":dr,"nrd":nza}

#==================================================================
#  
#  Write the files and get the modes:
#  
#==================================================================

print('Writing env files and getting the modes...')

nza = 101
dr = zeros(nza)

for i in range(nx):
   print(nx,i+1)
   Di = zbty[0,i]
   zi = array([0,Di])
   zarrayi = linspace(0.0,Di,nza)
   cdata = array([zi,c,cs,rho,apt,ast])
   layerp  = array([0, 0, Di])
   properties = array([Di, cp, 0, rhob, ab, 0])
   rd         = zarrayi
   nrd        = nza
   ssp_data = {"cdata":cdata,"type":mtype,"itype":itype,"nmesh":nmesh,"sigma":sigma,"clow":clow,"chigh":chigh,"zbottom":Di}
   bottom_data = {"n":1,"layerp":layerp,"layert":layert,"properties":properties,"bdata":bdata,"units":units,"bc":bc,"sigma":sigma}
   field_data = {"rmax":rmaxkm,"nrr":nra,"rr":yarraykm,"rp":rp,"np":np,"m":m,"rmodes":rmodes,"stype":stype,"thorpe":thorpe,"finder":finder,"rd":zarrayi,"dr":dr,"nrd":nza}
   wkrakenenvfil('dummy',case_title,source_data,surface_data,scatter_data,ssp_data, bottom_data,field_data)
   os.system("kraken.exe dummy")
   thecommand = "mv dummy.env range1_" + str(i+1) + ".env"
   os.system( thecommand )
   thecommand = "mv dummy.mod range1_" + str(i+1) + ".mod"
   os.system( thecommand )
   
print('Organizing files...');

for i in range(nx):
    print (nx,i+1)
    for j in range(2,ny+1):
        thecommand = "cp range1_" + str(i+1) + ".env range" + str(j) + "_" + str(i+1) + ".env"
        os.system( thecommand )
        thecommand = "cp range1_" + str(i+1) + ".mod range" + str(j) + "_" + str(i+1) + ".mod"
        os.system( thecommand )

#*******************************************************************************
#
# Write the 3D field file
#
#*******************************************************************************

nzr = 1
 
filename = 'range'

#*******************************************************************************
#
# options(1:3) 'STD' -> Standard (no horizontal refraction)
#              'PDQ' -> Quick (?)
#              'GBT' -> Gaussian Beam (yes horizontal refraction)
# options(4:4) 'T' -> Tesselation check
#              'I' -> Incoherent
#              'F' -> ?
# options(5:5) 'F' -> Space filling in far field
#              'M' -> Minimum width at Rmax
#
#*******************************************************************************

options    = 'GBTFM' # horizontal refraction
#options    = 'GBTFF' # horizontal refraction
nmodes     = 9999
nsources   = 1
nranges    = nra
ranges     = yarraykm
dtheta     = 1.0
thetamin   = 0
thetamax   = 360.0-dtheta 
thetas = arange(thetamin,thetamax+dtheta,dtheta); nthetas = thetas.size; nalphas = nthetas
alphamin   = thetamin 
alphamax   = thetamax
nsteps     = 2001
step       = rmax/nsteps
epmult     =  1.0; # ?
ntriangles =  2*( nx - 1 )*( ny - 1 ) # write_field3dflp.m line 118

#*******************************************************************************  
print( "Writing the flp file..." )

flpfil = '3dpe.flp'

fid = open(flpfil, 'w')

fid.write(case_title);fid.write('\n')
fid.write('\'');fid.write(options);fid.write('\'\n')
fid.write(str(nmodes));fid.write('\n')
fid.write('1\n')
fid.write(str(xs[0]));fid.write('\n')
fid.write('1\n');
fid.write(str(ys[0]));fid.write('\n')
fid.write('1\n')
fid.write(str(zs[0]));fid.write('\n')
fid.write('1\n');
fid.write(str(zarray[0]));fid.write('\n')
fid.write(str(nranges));fid.write('\n')
fid.write(str(min(yarraykm)));fid.write(' ')
fid.write(str(max(yarraykm)));fid.write(' /\n')
fid.write(str(nthetas));fid.write('\n')
fid.write(str(thetamin));fid.write(' ')
fid.write(str(thetamax));fid.write(' /\n')
fid.write(str(ny*nx));fid.write('\n')

for i in range(ny):

    for j in range(nx):
        
        filename_at_node = 'range' + str(i+1) + '_' + str(j+1)
        fid.write( str( x[j]/1000.0 ) ); fid.write(' ')
        fid.write( str( y[i]/1000.0 ) ); fid.write(' ')
        fid.write('\'');fid.write(filename_at_node);fid.write('\'\n')
	
fid.write(str(ntriangles));fid.write('\n')

# write_field3dflp.m lines 125 - 132

inode = 1

for iy in range( 1, ny ):
    for ix in range( 1, nx ):
        fid.write( str( inode          ) ); fid.write(' ')
        fid.write( str( inode + 1      ) ); fid.write(' ')
        fid.write( str( inode + nx     ) ); fid.write('\n')
        fid.write( str( inode + 1      ) ); fid.write(' ')
        fid.write( str( inode + nx     ) ); fid.write(' ')
        fid.write( str( inode + nx + 1 ) ); fid.write('\n')
        inode = inode + 1
    inode = inode + 1

fid.write(str(alphamin));fid.write(' ')
fid.write(str(alphamax));fid.write(' ')
fid.write(str(nalphas));fid.write('\n')
fid.write(str(step));fid.write(' ')
fid.write(str(nsteps));fid.write('\n')
fid.write(str(epmult));fid.write('\n')

fid.close()

print( "Running FIELD3D..." )

os.system("field3d.exe 3dpe > field3d.prt")

print( "Reading output data..." )

filename = '3dpe.shd'
pressure,geometry = readshd(filename,xs,ys)

thetas = geometry["thetas"]
rarray = geometry["rarray"]
zarray = geometry["zarray"]

R,Thetas= meshgrid(rarray,thetas)

X = R*cos( Thetas*pi/180.0 )
Y = R*sin( Thetas*pi/180.0 )

f = abs( thetas - 90.0 )

useless,i90deg = f.min(0),f.argmin(0)

pressure = squeeze( pressure )
tldisk = 20*log10( abs( pressure ) )

tldisk[:,0] = tldisk[:,1]

tl = squeeze( tldisk[i90deg, : ] )

mismatch = std( tlref[1:-1] - tl[1:-1] )

thetitle = 'KRAKEN3D - wedge problem @ ' + str(freq) + ' Hz'

figure(1)
pcolormesh(X,Y,tldisk,shading='auto'), colorbar()
xlabel('X (m)')
ylabel('Y (m)')
title( thetitle )

thetitle = 'KRAKEN3D @ ' + str(freq) + ' Hz'

figure(2)
grid(True)
plot( yarraykm, tlref, linewidth=2, label="3DPE" )
plot( yarraykm, tl, 'k--', linewidth=2, label="KRAKEN3D" )
legend()
ylim(-100,-20)
xlabel('ACROSS SLOPE RANGE [in km]')
ylabel('TLOSS [in dB]')
title( thetitle )

show()

print('done.')
