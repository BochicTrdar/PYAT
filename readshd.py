from numpy import *

def readshd( filename=None, xs=None, ys=None, freq=None ):
    #*******************************************************************************
    # Faro, dom 16 jun 2024 19:37:41 
    # Written by Orlando Camargo Rodriguez
    # Based on read_shd_bin.m by Michael Porter
    #*******************************************************************************
   fid = open(filename,'rb')
   recl  = int( fromfile( fid, int32, 1 ) )
   title = fid.read(80)
   fid.seek( 4*recl )
   PlotType = fid.read(10)
   fid.seek( 2*4*recl ); # reposition to end of second record
   Nfreq  = int(   fromfile( fid, int32  , 1 ) )
   Ntheta = int(   fromfile( fid, int32  , 1 ) )
   Nsx    = int(   fromfile( fid, int32  , 1 ) )
   Nsy    = int(   fromfile( fid, int32  , 1 ) )
   Nsz    = int(   fromfile( fid, int32  , 1 ) )
   Nrz    = int(   fromfile( fid, int32  , 1 ) )
   Nrr    = int(   fromfile( fid, int32  , 1 ) )
   atten  = float( fromfile( fid, float32, 1 ) )
   fid.seek( 3 * 4 * recl ); # reposition to end of record 3
   freqVec = fromfile( fid, float32, Nfreq  )
   fid.seek( 4 * 4 * recl ); # reposition to end of record 4
   theta  = fromfile( fid, float32, Ntheta )
   if  ( PlotType[ 0 : 1 ] != 'TL' ):
       fid.seek( 5 * 4 * recl ); # reposition to end of record 4
       Xs     = fromfile( fid, float32, Nsx )
       fid.seek( 6 * 4 * recl );  # reposition to end of record 5
       Ys     = fromfile( fid, float32, Nsy )
   else:   # compressed format for TL from FIELD3D
       fid.seek( 5 * 4 * recl ); # reposition to end of record 4
       Pos_S_x     = fromfile( fid, float32, 2 )
       Xs          = linspace( Pos_S_x[0], Pos_S_x[1], Nsx )
       fid.seek( 6 * 4 * recl ); # reposition to end of record 5
       Pos_S_y     = fromfile( fid, float32, 2 )
       Ys          = linspace( Pos_S_y[0], Pos_S_y[1], Nsy )
   fid.seek( 7 * 4 * recl ) # reposition to end of record 6
   zs = fromfile( fid, float32, Nsz )
   fid.seek( 8 * 4 * recl ) # reposition to end of record 7
   zarray =  fromfile( fid, float32, Nrz )
   fid.seek( 9 * 4 * recl ) # reposition to end of record 8
   rarray =  fromfile( fid, float64, Nrr )
   if PlotType == 'rectilin  ':
       pressure = zeros( (Ntheta, Nsz, Nrz, Nrr) ) + 1j*zeros( (Ntheta, Nsz, Nrz, Nrr) )
       Nrcvrs_per_range = Nrz
   elif PlotType == 'irregular ':
       pressure = zeros( (Ntheta, Nsz,   1, Nrr) ) + 1j*zeros( (Ntheta, Nsz,   1, Nrr) )
       Nrcvrs_per_range = 1
   else:
       pressure = zeros( (Ntheta, Nsz, Nrz, Nrr) ) + 1j*zeros( (Ntheta, Nsz, Nrz, Nrr) )
       Nrcvrs_per_range = Nrz
   if isnan( xs ):
      ifreq = 0
      if isnan(freq) == False: 
         freqdiff = abs( freqVec - freq )
         ifreq = freqdiff.argmin( )
      for itheta in range(Ntheta):
          for isz in range( Nsz ):
              for irz in range( Nrcvrs_per_range ):
                  recnum = 10 + ifreq*Ntheta*Nsz*Nrcvrs_per_range + itheta*Nsz*Nrcvrs_per_range + isz* Nrcvrs_per_range + irz
                  status = fid.seek( recnum * 4 * recl ) # Move to end of previous record
                  if ( status == -1 ):
                     print('Seek to specified record failed in readshd...')
                  temp = fromfile( fid, float32, 2 * Nrr ) # Read complex data
                  indexes = arange(0,2*Nrr,2)
                  pressure[ itheta, isz, irz, : ] = temp[indexes] + 1j*temp[indexes+1]

   else:    
       xdiff = abs( Xs - xs * 1000.0 )
       idxX  = xdiff.argmin(0)
       ydiff = abs( Ys - ys * 1000.0 )
       idxY  = ydiff.argmin(0)
       for itheta in range(Ntheta):
           for isz in range(Nsz):
               for irz in range( Nrcvrs_per_range ):
                   recnum = 10 + idxX * Nsy * Ntheta * Nsz * Nrcvrs_per_range + idxY * Ntheta * Nsz * Nrcvrs_per_range + itheta * Nsz * Nrcvrs_per_range + isz * Nrcvrs_per_range + irz
                   status = fid.seek( recnum * 4 * recl ) # Move to end of previous record
                   if ( status == -1 ):
                      print('Seek to specified record failed in read_shd_bin')
                   temp = fromfile( fid, float32, 2 * Nrr ) # Read complex data
                   indexes = arange(0,2*Nrr,2)
                   pressure[ itheta, isz, irz, : ] = temp[indexes] + 1j*temp[indexes+1]
               
       fid.close()
   geometry = {"zs":zs, "f":freqVec,"thetas":theta,"rarray":rarray,"zarray":zarray}
   return pressure,geometry
