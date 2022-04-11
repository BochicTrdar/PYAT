from numpy import *

def readmod( filename=None, freq=None ):
    #*******************************************************************************
    # Faro, Qua Mai 26 19:06:46 WEST 2021
    # Written by Orlando Camargo Rodriguez
    # Based on read_modes.m & read_modes_bin by Michael Porter
    #*******************************************************************************
    Modes = []
    bc    = []
    iRecProfile = 1
    fid = open(filename,'rb')
#======================================================================
    lrecl  = 4*int( fromfile( fid, int32, 1 ) )
    rec = iRecProfile - 1
    fid.seek( rec * lrecl + 4 )
    title = fid.read(80)
    Nfreq  =   int( fromfile( fid,   int32, 1 ) )
    Nmedia =   int( fromfile( fid,   int32, 1 ) )
    Ntot   =   int( fromfile( fid,   int32, 1 ) )
    NMat   =   int( fromfile( fid,   int32, 1 ) )
    if Ntot < 0: 
       return Modes,bc
    rec = iRecProfile  
    fid.seek( rec * lrecl )
    N = zeros( Nmedia )
    Mater = []
    for Medium in range(Nmedia):
        N[ Medium ] = int( fromfile( fid, int32, 1 ) )
        Mi = fid.read(8)
        Mater.append( Mi )    
    # depth and rho
    rec = iRecProfile + 1
    fid.seek( rec * lrecl )
    bulk = zeros((2,Nmedia))
    for i in range(Nmedia):
        bulk[0,i]      = float( fromfile( fid, float32, 1 ) )
        bulk[1,i]      = float( fromfile( fid, float32, 1 ) )
    depth = bulk[0,]
    rho   = bulk[1,]
    # frequencies
    rec = iRecProfile + 2
    fid.seek( rec * lrecl )
    freqVec = zeros(Nfreq)
    for i in range(Nfreq):
        freqVec[i] = float( fromfile( fid, float32, 1 ) )
    # z
    rec = iRecProfile + 3
    fid.seek( rec * lrecl )
    z = zeros(Ntot)
    for i in range(Ntot):
        z[i] = float( fromfile( fid, float32, 1 ) )
    # read in the modes
    # (sorry, 1 frequency only!)
    # number of modes, m
    iRecProfile = iRecProfile + 4
    rec = iRecProfile
    fid.seek( rec * lrecl )
    M = int( fromfile( fid, int32, 1 ) )
    rec = iRecProfile + 1
    # Top
    Top_bc    = fid.read(1)
    cp        = zeros( 2 )
    cs        = zeros( 2 )
    cp[0]      = float( fromfile( fid, float32, 1 ) )
    cp[1]      = float( fromfile( fid, float32, 1 ) )
    Top_cp    = complex( cp[ 0 ], cp[ 1 ] )
    cs[0]      = float( fromfile( fid, float32, 1 ) )
    cs[1]      = float( fromfile( fid, float32, 1 ) )
    Top_cs    = complex( cs[ 0 ], cs[ 1 ] )
    Top_rho   = float( fromfile( fid, float32, 1 ) )
    Top_depth = float( fromfile( fid, float32, 1 ) )
    # Bottom
    Bot_bc    = fid.read(1)
    cp[0]      = float( fromfile( fid, float32, 1 ) )
    cp[1]      = float( fromfile( fid, float32, 1 ) )
    Bot_cp    = complex( cp[ 0 ], cp[ 1 ] )
    cs[0]      = float( fromfile( fid, float32, 1 ) )
    cs[1]      = float( fromfile( fid, float32, 1 ) )
    Bot_cs    = complex( cs[ 0 ], cs[ 1 ] )
    Bot_rho   = float( fromfile( fid, float32, 1 ) )
    Bot_depth = float( fromfile( fid, float32, 1 ) )
    rec = iRecProfile
    fid.seek( rec * lrecl )
    if ( M == 0 ):
       phi = []
       k   = []
    else:
       phi = zeros( (NMat,M) ) + 1j*zeros( (NMat,M) )
       for ii in range(M):
           rec = iRecProfile + 1 + ii
           fid.seek( rec * lrecl )
           for jj in range(NMat):
               phi1 = float( fromfile( fid, float32, 1 ) )
               phi2 = float( fromfile( fid, float32, 1 ) )
               phi[ jj, ii-1 ] = complex( phi1, phi2 )
       rec = iRecProfile + 2 + M
       fid.seek( rec * lrecl )
       k = zeros( M ) + 1j*zeros( M )
       for kk in range(M):
           k1 = float( fromfile( fid, float32, 1 ) )
           k2 = float( fromfile( fid, float32, 1 ) )
           k[ kk ] = complex( k1, k2 )
#======================================================================
    fid.close()
    
    bc    = {"Top_bc":Top_bc,"Bot_bc":Bot_bc,"Top_cp":Top_cp,"Bot_cp":Bot_cp,"Top_cs":Top_cs,"Bot_cs":Bot_cs,
             "Top_depth":Top_depth,"Bot_depth":Bot_depth,"Top_rho":Top_rho,"Bot_rho":Bot_rho}
    Modes = {"freq":freqVec, "Nmedia":Nmedia, "Ntot":Ntot, "NMat":NMat,"N":N,"Mater":Mater,"depth":depth,"rho":rho,"z":z,
             "M":M,"phi":phi,"k":k}
    
    return Modes,bc
