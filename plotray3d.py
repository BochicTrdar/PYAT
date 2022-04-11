from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

def plotray3d(filename=None):    
    #*******************************************************************************
    # Faro, Qua Mai 26 19:07:39 WEST 2021
    # Written by Orlando Camargo Rodriguez
    # Based on plotray3d.m by Michael Porter
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
    xmin = +1e9;
    xmax = -1e9;
    ymin = +1e9;
    ymax = -1e9;
    zmin = +1e9;
    zmax = -1e9;
    
    fig = figure(1)
    ax = fig.gca(projection='3d')
    
    for isx in range(Nsx):
        for isy in range(Nsy):        
	    for ibeam2 in range(Nbeta):
                for ibeam in range(Nalpha):
		    theline   = str( fid.readline() )
		    l = len( theline )
		    if l > 0:
	               alpha0    = float( theline )
 	               theline = str( fid.readline() )
 	               datai = theline.split()
 	               nsteps    = int( datai[0] )
 	               NumTopBnc = int( datai[1] )
 	               NumBotBnc = int( datai[2] )
	               x = zeros(nsteps)
		       y = zeros(nsteps)
	               z = zeros(nsteps)
	               for j in range(nsteps):
	                   theline = str(fid.readline())
                           xyz = theline.split()
                           x[j] = float( xyz[0] )
                           y[j] = float( xyz[1] )
		           z[j] = float( xyz[2] )        
                       xmin = min( [ min(x), xmin ] )
                       xmax = max( [ max(x), xmax ] )
                       ymin = min( [ min(y), ymin ] )
                       ymax = max( [ max(y), ymax ] )		       
                       zmin = min( [ min(z), zmin ] )
                       zmax = max( [ max(z), zmax ] )
		       ax.plot(x,y,-z)
    
    fid.close()
