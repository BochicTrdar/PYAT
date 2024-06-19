from scipy.special import hankel1,hankel2
from scipy.io import *
from numpy    import *
from readshd  import *
from scalecol import *
from taper    import *

def fieldsco( filename=None, zs=None, zarray=None, freqVec=None ):
# Calculates the field using the Green's function produced by SCOOTER
#
# This version uses the trapezoidal rule directly to do a DFT, rather than an FFT
#
# You could save storage by doing the transform for one row at a time out
# of the Green's function file, i.e. a single source/receiver depth.
# This version does the entire matrix at once. If that fits comfortably in
# memory then it is potentially faster in terms of using multiple cores
# efficiently
#
# usage:
#    fieldsco( filename )
#  or
#    fieldsco( filename, PlotTitle, freq, atten, Pos, Gtemp, option, Rminkm, Rmaxkm, Nrr )
#
# You must include the file extension, if it exists
#
# Optionally reads control info from a file fields.flp
# mbp, Dec. 2012 based on fieldsco

# phase speed limits for optional tapering
# this is for user play (at your own risk)
# Based on fieldsco.m by Michael Porter

    xs = nan
    ys = nan
    cmin = 1e-10
    cmax = 1e30
    fileroot = filename[ 0 : -4 ] # *.grn

#   Read the *.flp:
    fileflp =  fileroot + '.flp'
    fid = open(fileflp, 'r')
    Opt     = fid.readline()
    theline = fid.readline(); Nrr = int( theline )
    theline = fid.readline(); i = theline.find('/'); rr = fromstring( theline[0:i], dtype=float, sep=' ' )
    fid.close()
    rarraykm = linspace(rr[0],rr[1],Nrr); Rr = 1000*rarraykm; 
    if ( Rr[0] < 1e-6 ):
         Rr[0] = 1.0
    i = Opt.find('*')
    SBP = Opt[i]
    if ( SBP == '*' ):
       print( '-----------------------------------' )
       print( 'Using source beam pattern file' )
       sbpfil = fileroot + '.sbp'
       fid = open( sbpfil, 'r' )    
       theline = fid.readline; NSBPPts = int( theline )
       print( 'Number of source beam pattern points =', NSBPPts )
       SrcBmPat = zeros( (NSBPPts, 2) )
       print( ' ' )
       print( ' Angle (degrees)  Power (dB)' )
       for I in range( NSBPPts):
           theline = fid.readline()
           print( theline )
           SrcBmPat[ I, : ] = fromstring( theline, dtype=float, sep=' ' )
       fid.close()
    else:   # no pattern given, use omni source pattern
       SrcBmPat = zeros( (2, 2) )
       SrcBmPat[ 0, : ] = array([ -180.0, 0.0 ])
       SrcBmPat[ 1, : ] = array([  180.0, 0.0 ])

    SrcBmPat[ :, 1 ] = 10**( SrcBmPat[ :, 1 ] / 20 )  # convert dB to linear scale
    
    Nsz   = size( zs     )  # Number of source   depths
    Nrz   = size( zarray )  # Number of receiver depths
    Nfreq = size( freqVec ) # Number of frequencies
    
    pressure = zeros( (Nfreq, Nsz, Nrz, Nrr) ) + 1j*zeros( (Nfreq, Nsz, Nrz, Nrr) )
    
    for ifreq in range(Nfreq):
        freq   = freqVec[ ifreq ]
        omega  = 2.0 * pi * freq
        kleft  = omega / cmax # left  limit for tapering
        kright = omega / cmin # right limit for tapering
   
   # read in the Green's function for this frequency
        Gtemp,geometry = readshd( filename, xs, ys, freq );
        cVec = geometry['rarray']
        k = 2 * pi * freq / cVec

   # for a SPARC run the wavenumber vector is independent of frequency
#   if ( contains( PlotTitle( 1 : 6 ), 'SPARC' ) )
#      k    = 2 * pi * freq0 ./ cVec;
#   end
   
        Nk = size( k ) # Number of wavenumbers
   
        deltak = ( k[-1] - k[0] ) / ( Nk - 1 )

        if ( Rr[-1] * deltak > 10 ):
           print( 'The wavenumber sampling is too coarse to accurately calculate the field at the largest range' )
   
#   if ( contains( PlotTitle( 1 : 7 ), 'SCOOTER' ) )
#      atten  = deltak;   % this is calculated here because SCOOTER doesn't write it for each frequency
#   end
        atten = deltak
        ck = k + 1j * atten
#   x  = ck.' * abs( Rr' );   % this is a matrix of size nk * nrr
        RR,CK = meshgrid(Rr,ck); x = RR*CK
        if Opt[1] == 'X': # case 'X' => line source
            X  = exp( -1j * x )  # e^( -i k r ) matrix
            X2 = exp(  1j * x )  # e^( +i k r ) matrix
            factor1 = 1;
            factor2 = deltak/sqrt( 2 * pi )
        elif Opt[1] == 'R': # point source
            X  = exp( -1j * ( x - pi / 4 ) ) # e^( -i k r ) matrix
            X2 = exp(  1j * ( x - pi / 4 ) ) # e^( +i k r ) matrix
            factor1 = sqrt( ck )
            factor2 = deltak / sqrt( 2 * pi * abs( Rr ) )
        elif Opt[1] == 'S': # point source with cylindrical spreading removed
            X  = exp( -1j * ( x - pi / 4 ) ) # e^( -i k r ) matrix
            X2 = exp(  1j * ( x - pi / 4 ) ) # e^( +i k r ) matrix
            factor1 = sqrt( ck )
            factor2 = deltak / sqrt( 2 * pi )
        elif Opt[1] == 'H': # exact Bessel transform
            X  = hankel1( 0, x )  # Hankel function
            X2 = hankel2( 0, x )  # Hankel function
            factor1 = ck
            factor2 = deltak / 2
        else:
            print( 'fieldsco: Unknown first letter option in second line of .flp file' )
   
        for isz in range(Nsz):
             G = squeeze( Gtemp[ isz, 0, :, : ] )
             s = shape( G )
#             if s[1] == 1 # if G is a vector, expand it into a matrix with one row
#                G = reshape( G, 1, length( G ) )
#      end
      
      # apply the source beam pattern
             if ( SBP == '*' ):
                c = 1500  # reference sound speed, should be speed at the source depth
                kz2 = omega**2 / c**2 - k**2 # vertical wavenumber squared
                kz2 = where( kz2 < 0, 0, kz2 ) # remove negative values
                theta = arctan( sqrt( kz2*180/pi )/ k ) # calculate the angle in degrees
                S = interp( theta, SrcBmPat[ :, 0 ], SrcBmPat[ : , 1 ] )
                G = scalecol( G, S )
      
             G = taper( G, k, Nk, kleft, kright )
             G = where( abs( G ) < 1e-40, 0, G ) # avoid underflows--- they slow the following down a huge amount
        
             if Opt[2] == 'P':
                 G = scalecol( G, factor1 )            
                 Y = -matmul( G ,  X )# here's where the DFT is done
             elif Opt[2] == 'N': # check this !!!
                 G = scalecol( G, factor1 )                        
                 Y = -matmul( G , X2 )
             elif Opt[2] == 'N': # check this !!!
                 G = scalecol( G, factor1 )                       
                 Y = -matmul( G , ( X + X2 ) )
             else:
                 print( 'fieldsco: Unknown second letter option in second line of .flp file' )
        
        Y = scalecol( Y, factor2 )

        pressure[ ifreq, isz, :, : ] = Y # cylindrical spreading, etc.

        print( 'Transform completed for source depth',zs[isz],'frequency', freq )

    PlotType = "rectilin  "
#    data = {"PlotTitle": "SCOOTER", "PlotType":PlotType,"freqVec":freq,"atten":0,'Pos':Rr,'pressure':pressure}
    PlotTitle = 'SCOOTER'
    atten = 0.0
    data = {'PlotTitle':PlotTitle,'PlotType':PlotType,'freqVec':freq,"atten":atten,'pressure':pressure}
    filemat = fileroot + '.shd.mat'
    savemat( filemat, data )
    
    return
