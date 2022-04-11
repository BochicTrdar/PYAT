from numpy import *

def read_arrivals_asc(filename=None):    
    #*******************************************************************************
    # Faro, Qua 26 Mai 2021 20:54:06 WEST 
    # Written by Orlando Camargo Rodriguez
    # Based on read_arrivals_asc.m by Michael Porter
    #*******************************************************************************
    Arr = []
    Pos = []
    Narrmx = 100
    maxnarr = 0
    fid = open(filename,'r')
    flag = str( fid.readline() )
    if flag[2] == '2':
       freq  = float( fid.readline() )
       theline = str( fid.readline() )
       datai = theline.split()
       Nsz   = int( datai[0] )
       source_depths = zeros(Nsz)
       for i in range(Nsz):
           source_depths[i] = float( datai[i+1] )

       theline = str( fid.readline() )
       datai = theline.split()
       Nrz   = int( datai[0] )
       receiver_depths = zeros(Nrz)
       for i in range(Nrz):
           receiver_depths[i] = float( datai[i+1] )

       theline = str( fid.readline() )
       datai = theline.split()
       Nrr   = int( datai[0] )
       receiver_ranges = zeros(Nrr)
       for i in range(Nrr):
           receiver_ranges[i] = float( datai[i+1] )

       Narr      = zeros( (Nrr,         Nrz, Nsz) )	
       A         = zeros( (Nrr, Narrmx, Nrz, Nsz) ) + 1j*zeros( (Nrr, Narrmx, Nrz, Nsz) )
       delay     = zeros( (Nrr, Narrmx, Nrz, Nsz) ) + 1j*zeros( (Nrr, Narrmx, Nrz, Nsz) )
       SrcAngle  = zeros( (Nrr, Narrmx, Nrz, Nsz) )
       RcvrAngle = zeros( (Nrr, Narrmx, Nrz, Nsz) )
       NumTopBnc = zeros( (Nrr, Narrmx, Nrz, Nsz) )
       NumBotBnc = zeros( (Nrr, Narrmx, Nrz, Nsz) )

       for isd in range(Nsz):
           Narrmx2 = int( fid.readline() )
           for ird in range(Nrz):
               for irr in range(Nrr):
	           narr = int( fid.readline() )
                   Narr[ irr, ird, isd ] = narr
                   maxnarr = max( narr, maxnarr )
		   if narr > 0:
		      narr = min( narr, Narrmx )
		      for k in range(narr):
                          theline = str( fid.readline() )
                          datai = theline.split()
		          amp   = float( datai[0] )
		          phase = float( datai[1] )
                          A[ irr, k, ird, isd ] = amp*exp( 1j*phase*pi/180.0 )
		          rtau = float( datai[2] )
		          itau = float( datai[3] )
 		          delay[ irr, k, ird, isd ] = rtau + 1j*itau
		          source_angle = float( datai[4] ) 
                          SrcAngle[ irr, k, ird, isd ] = source_angle
		          receiver_angle = float( datai[5] ) 
		          RcvrAngle[ irr, k, ird, isd ] = receiver_angle
		          bounces = int( datai[6] )
                          NumTopBnc[ irr, k, ird, isd ] = bounces
		          bounces = int( datai[7] )		       
                          NumBotBnc[ irr, k, ird, isd ] = bounces
       A         = A[        :,0:maxnarr,:,:]
       delay     = delay[    :,0:maxnarr,:,:]
       SrcAngle  = SrcAngle[ :,0:maxnarr,:,:]
       RcvrAngle = RcvrAngle[:,0:maxnarr,:,:]
       NumTopBnc = NumTopBnc[:,0:maxnarr,:,:]
       NumBotBnc = NumBotBnc[:,0:maxnarr,:,:]       
       
       Pos = {'freq':freq,'source_depths':source_depths,'receiver_depths':receiver_depths,'receiver_ranges':receiver_ranges}
       Arr = {'Narr':Narr,'A':A,'delay':delay,'SrcAngle':SrcAngle,'RcvrAngle':RcvrAngle,'NumTopBnc':NumTopBnc,'NumBotBnc':NumBotBnc}
    else:

       freq  = float( fid.readline() ) 

       theline = str( fid.readline() )
       datai = theline.split()
       Nsx   = int( datai[0] )
       sourcex = zeros(Nsx)
       for i in range(Nsx):
           sourcex[i] = float( datai[i+1] )

       theline = str( fid.readline() )
       datai = theline.split()
       Nsy   = int( datai[0] )
       sourcey = zeros(Nsy)
       for i in range(Nsy):
           sourcey[i] = float( datai[i+1] )

       theline = str( fid.readline() )
       datai = theline.split()
       Nsz   = int( datai[0] )
       sourcez = zeros(Nsz)
       for i in range(Nsz):
           sourcez[i] = float( datai[i+1] )

       theline = str( fid.readline() )
       datai = theline.split()
       Nrz   = int( datai[0] )
       receiver_depths = zeros(Nrz)
       for i in range(Nrr):
           receiver_depths[i] = float( datai[i+1] )

       theline = str( fid.readline() )
       datai = theline.split()
       Nrr = int( datai[0] )
       receiver_ranges = zeros(Nrr)
       for i in range(Nrr):
           receiver_ranges[i] = float( datai[i+1] )

       theline = str( fid.readline() )
       datai = theline.split()
       Nrtheta = int( datai[0] )
       receiver_thetas = zeros(Nrtheta)
       for i in range(Nrtheta):
           receiver_thetas[i] = float( datai[i+1] )

       Narr      = zeros( (Nrr,         Nrz, Nrtheta, Nsz) )	
       A         = zeros( (Nrr, Narrmx, Nrz, Nrtheta, Nsz) ) + 1j*zeros( (Nrr, Narrmx, Nrz, Nrtheta, Nsz) )
       delay     = zeros( (Nrr, Narrmx, Nrz, Nrtheta, Nsz) ) + 1j*zeros( (Nrr, Narrmx, Nrz, Nrtheta, Nsz) )
       SrcDeclAngle = zeros( (Nrr, Narrmx, Nrz, Nrtheta, Nsz) )
       SrcAzimAngle = zeros( (Nrr, Narrmx, Nrz, Nrtheta, Nsz) )
       RcvrDeclAngle= zeros( (Nrr, Narrmx, Nrz, Nrtheta, Nsz) )
       RcvrAzimAngle= zeros( (Nrr, Narrmx, Nrz, Nrtheta, Nsz) )
       NumTopBnc    = zeros( (Nrr, Narrmx, Nrz, Nrtheta, Nsz) )
       NumBotBnc    = zeros( (Nrr, Narrmx, Nrz, Nrtheta, Nsz) )

       for isd in range(Nsz):
           Narrmx2 = int( fid.readline() )
           for irtheta in range(Nrtheta):
               for irz in range(Nrz):
                   for irr in range(Nrr):
                       Narr[ irr, ird, irtheta, isd ] = narr
                       maxnarr = max( narr, maxnarr )
		       if narr > 0:
		          narr = min( narr, Narrmx )
		          for k in range(narr):
                              theline = str( fid.readline() )
                              datai = theline.split()
		              amp   = float( datai[0] )
		              phase = float( datai[1] )
                              A[ irr, k, ird, irtheta, isd ] = amp*exp( 1j*phase*pi/180.0 )
		              rtau = float( datai[2] )
		              itau = float( datai[3] )
 		              delay[ irr, k, ird, irtheta, isd ] = rtau + 1j*itau
		              theangle = float( datai[4] ) 
                              SrcDeclAngle[ irr, k, ird, irtheta, isd ] = theangle
		              theangle = float( datai[5] ) 
		              SrcAzimAngle[ irr, k, ird, irtheta, isd ] = theangle
		              theangle = float( datai[6] ) 
                              RcvrDeclAngle[ irr, k, ird, irtheta, isd ] = theangle
		              theangle = float( datai[7] )
                              RcvrAzimAngle[ irr, k, ird, irtheta, isd ] = theangle
                              bounces = int( datai[8] )
                              NumTopBnc[ irr, k, ird, irtheta, isd ] = bounces
		              bounces = int( datai[9] )		       
                              NumBotBnc[ irr, k, ird, irtheta, isd ] = bounces
       A            = A[            :,0:maxnarr,:,:,:]
       delay        = delay[        :,0:maxnarr,:,:,:]
       SrcDeclAngle = SrcDeclAngle[ :,0:maxnarr,:,:,:]
       SrcAzimAngle = SrcAzimAngle[ :,0:maxnarr,:,:,:]
       RcvrDeclAngle= RcvrDeclAngle[:,0:maxnarr,:,:,:]
       RcvrAzimAngle= RcvrAzimAngle[:,0:maxnarr,:,:,:]
       NumTopBnc = NumTopBnc[:,0:maxnarr,:,:,:]
       NumBotBnc = NumBotBnc[:,0:maxnarr,:,:,:]
       Pos = {'freq':freq,'sourcex':sourcex,'sourcey':sourcey,'sourcez':sourcez,'receiver_depths':receiver_depths,'receiver_ranges':receiver_ranges,'receiver_thetas':receiver_thetas}
       Arr = {'Narr':Narr,'A':A,'delay':delay,'SrcDeclAngle':SrcDeclAngle,'SrcAzimAngle':SrcAzimAngle,'RcvrDeclAngle':RcvrDeclAngle,'RcvrAzimAngle':RcvrAzimAngle,'NumTopBnc':NumTopBnc,'NumBotBnc':NumBotBnc}

    fid.close()
    return Arr,Pos
