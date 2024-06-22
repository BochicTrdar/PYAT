from numpy import *
from pylab import *

def readbrc(filename=None):    
    #*******************************************************************************
    # Boston, sex jun 30 14:09:20 WEST 2017
    # Written by Tordar
    #*******************************************************************************
    
    fid = open(filename,'r')
    nthetas = int( fid.readline() )
#   print nthetas
    thetas = zeros( nthetas ) 
    absR   = zeros( nthetas )
    phaseR = zeros( nthetas )
    for i in range(nthetas):
        theline = str( fid.readline() )
        datai = theline.split()
        thetas[i] = float( datai[0] )
        absR[i]   = float( datai[1] )
        phaseR[i] = float( datai[2] )     
    fid.close()
    R = absR*exp( 1j*phaseR )
    return thetas,R
