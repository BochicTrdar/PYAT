from numpy import *
from matplotlib.pyplot import *

def plotray(filename=None):    
    #*******************************************************************************
    # Faro, Seg 11 Abr 2022 12:40:41 WEST 
    # Written by Orlando Camargo Rodriguez
    # Based on plotray.m by Michael Porter
    #*******************************************************************************

    Nsxyz       = zeros(3) 
    NBeamAngles = zeros(2)
    
    fid = open(filename,'r')
    title =        fid.readline()
    freq  = float( fid.readline() )
    theline = str( fid.readline() )
    datai = theline.split()
    Nsxyz[0] = int( datai[0] )
    Nsxyz[1] = int( datai[1] )
    Nsxyz[2] = int( datai[2] )
    theline = str( fid.readline() )
    datai = theline.split()
    NBeamAngles[0] = int( datai[0] )
    NBeamAngles[1] = int( datai[1] )
    DEPTHT = float( fid.readline() )
    DEPTHB = float( fid.readline() )
    Type   = fid.readline()
    Nsx = int( Nsxyz[0] )
    Nsy = int( Nsxyz[1] )
    Nsz = int( Nsxyz[2] )
    Nalpha = int( NBeamAngles[0] )
    Nbeta  = int( NBeamAngles[1] )
    # axis limits
    rmin =  1.0e9
    rmax = -1.0e9
    zmin =  1.0e9
    zmax = -1.0e9
    for isz in range(Nsz):
        for ibeam in range(Nalpha):
	   #alpha0    = float( fid.readline() )
            theline = str( fid.readline() )
            l = len( theline )
            if l > 0:
               alpha0 = float( theline )
               theline = str( fid.readline() )
               datai = theline.split()
               nsteps    = int( datai[0] )
               NumTopBnc = int( datai[1] )
               NumBotBnc = int( datai[2] )
               r = zeros(nsteps)
               z = zeros(nsteps)
               for j in range(nsteps):
                   theline = str(fid.readline())
                   rz = theline.split()
                   r[j] = float( rz[0] )
                   z[j] = float( rz[1] )        
               rmin = min( [ min(r), rmin ] )
               rmax = max( [ max(r), rmax ] )
               zmin = min( [ min(z), zmin ] )
               zmax = max( [ max(z), zmax ] )
               plot( r, -z )
               axis([rmin,rmax,-zmax,-zmin])
    fid.close()
