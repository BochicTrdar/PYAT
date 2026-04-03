# ==================================================================
#  
#  KRAKEN: range dependent calculations
#  Mexilhoeira Grande, sex 03 abr 2026 15:12:21 
#  Written by Tordar 
#  
# ==================================================================

import numpy as np
import matplotlib.pyplot as plt
import warnings
from os import system
from scipy.io import loadmat
from read_modes    import * 
from read_shd      import *
from wkrakenenvfil import *

# It works like this: 
# (1) load the transect;
# (2) write a *.env for every pair (rbottom,zbottom) along the transect;
# (3) merge all the *.env into a single "BIG" env; 
# (4) write the *.flp for RD calculations;
# (5) get the field.

warnings.filterwarnings('ignore') # Perhaps not a good idea, but I hate warnings, anyway... 

# Geometry parameters
thickness = 3.3
cp  = np.array( [1505, 1556, 1576] )
ap  = np.array( [0.11, 0.18] )
rho = np.array( [2.0 , 1.6 ] )

first_hyd = 17.7
dhyd = 2
last_hyd = 111.7

zarray = [i for i in range(int(first_hyd), int(last_hyd) + 1, int(dhyd))]
# More precise way to create zarray:
# zarray = list(range(int(first_hyd), int(last_hyd) + int(dhyd), int(dhyd)))
# Or using numpy for floating point precision:
# import numpy as np
# zarray = np.arange(first_hyd, last_hyd + dhyd, dhyd)

nza = len(zarray)

freq = 169
w = 2*np.pi*freq

case_title = 'TRANSECT'

# =====================================================================
# General stuff goes here:

zs = np.array([10.0]) # Single source
nzs = len(zs)

source_data = {
    'n': nzs,
    'zs': zs,
    'f': freq
}

surface_data = {
    'bc': 'V',
    'properties': [],  # Not required due to vacuum over surface
    'reflection': []   # Not required for this case
}

scatter_data = {
    'bumden': [],  # Bump density in ridges/km
    'eta'   : [],  # Principal radius 1 of bump
    'xi'    : []   # Principal radius 2 of bump
}

nmesh = 0

layert     = 'HH'

units      = 'W';
bc	   = 'A';
sigma      = 0.0 # Interfacial roughness

nra = 1

field_data = {
    'rmax': 1000, # defined, but not used...
    'rr': np.array([1]), # defined, but not used...
    'nrr': nra,
    'rp': 0,  # Range of first profile
    'np': 1,  # Number of profiles
    'm': 999,
    'rmodes': 'A',
    'stype': 'R',
    'thorpe': 'T',
    'finder': ' ',
    'rd': zarray,
    'dr': np.zeros_like(zarray),
    'nrd': nza
}
# =====================================================================

data = loadmat('atransect.mat') 

ri = np.squeeze( data['ri'] ); rikm = ri/1000
zi = np.squeeze( data['zi'] )
nenvs = np.size(ri) # Aaaaattention: this is less than 100... 

ssp = np.loadtxt('transect.ssp')

zssp = ssp[:,0]
cssp = ssp[:,1]

npoints = 101 # Don't like it? Use your own one... 

for ienv in range(nenvs):
    print(nenvs,ienv+1)
    Dmax = zi[ienv]
    zenv = np.linspace(0,Dmax,npoints)
    cenv = np.interp(zenv, zssp, cssp)
    csw  = np.zeros_like(zenv)
    rhow = np.ones_like( zenv)
    apw  = np.zeros_like(zenv)
    asw  = np.zeros_like(zenv)
    ssp_data = {
        'type': 'H',
        'itype': 'N',
        'nmesh': nmesh,  # Number of mesh points (about 10 per vertical wavelength)
        'sigma': 0,      # RMS roughness at the surface
        'clow': 0.0,
        'chigh': 5000.0,
        'cdata': np.vstack([zenv, cenv, csw, rhow, apw, asw]),
        'zbottom': Dmax
    }
    layer_info = np.array([[Dmax       ,cp[0],0,rho[0],ap[0],0], # [z,cp,cs,RHO,ap,as]
                        [Dmax+thickness,cp[2],0,rho[1],ap[1],0]])
    layerp     = np.array([[nmesh, 0.0, Dmax+thickness]])
    properties = np.array( layer_info[1,] )
    m1 = np.array( [layer_info[0,],layer_info[0,]] )
    bdata        = np.array( [m1] )
    bdata[0,1,0] = layer_info[1,0]
    bottom_data = {"n":2,"layerp":layerp,"layert":layert,"properties":properties,"bdata":bdata,"units":units,"bc":bc,"sigma":sigma}
    wkrakenenvfil('dummy', case_title, source_data, surface_data, scatter_data, ssp_data, bottom_data, field_data)
    if ienv < 10:
       system("mv dummy.env e0" + str(ienv) + ".env")
    else:
       system("mv dummy.env e"  + str(ienv) + ".env")
       
system('cat e??.env > atransect.env')
system('rm e??.env')

#==================================================================
#  
#  Write the flp file (clumsy, but functional):
#  
#==================================================================

max_number_of_modes = 201
nzs = 1 # One source 
nra = 201
Rmax = max( ri )
rarray = np.linspace(0,Rmax,nra); rarraykm = rarray/1000

fid = open('atransect.flp','w')
fid.write('TRANSECT Calculations\n')
fid.write('RA\n')
fid.write(str(max_number_of_modes)); fid.write('\n') 
fid.write(str(nenvs));fid.write('\n') 
fid.write(str(rikm[0]));fid.write(" "); fid.write(str(rikm[-1])); fid.write(" /\n")
fid.write(str(nra)); fid.write('\n') 
fid.write(str(rarraykm[0]));fid.write(" "); fid.write(str(rarraykm[-1])); fid.write(" /\n")
fid.write(str(nzs)); fid.write('\n') 
fid.write(str(zs[0])); fid.write(' /\n') 
fid.write(str(nza)); fid.write('\n') 
fid.write(str(zarray[0]));fid.write(" "); fid.write(str(zarray[-1])); fid.write(" /\n")
fid.write(str(nza)); fid.write('\n') # Yup, nza again... useful for non-vertical arrays  
fid.write('0.0 0.0 /\n')
fid.close( )

# Calculations: 

system('krakenc.exe atransect')
system('field.exe atransect')

# Plotting: 

PlotTitle, PlotType, freqVec, freq0, atten, Pos, p = read_shd( 'atransect' )

p = np.squeeze( p, axis=(0,1) )
p = np.where( p == 0, np.nan, p )
tl = -20*np.log10( abs( p ) )


plt.figure(1,dpi=300)
plt.plot(ri,zi)
plt.xlabel('Range (m)')
plt.ylabel('Depth (m)')
plt.title('TRANSECT')
plt.ylim([max(zi),0])
plt.grid(True)

plt.figure(2,dpi=300)
plt.pcolormesh(rarraykm,zarray,tl,vmin=60,vmax=90,cmap='jet_r',shading='auto')
cbar = plt.colorbar()
cbar.ax.invert_yaxis()
plt.xlim(rarraykm[0],rarraykm[-1])
plt.ylim(zarray[-1],zarray[0])
plt.xlabel('Range (km)')
plt.ylabel('Depth (m)')
plt.title('KRAKEN Transect calculations')

plt.show()

print("done.")
