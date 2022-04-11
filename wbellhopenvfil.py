def wbellhopenvfil(filename=None, thetitle=None, source_info=None, surface_info=None, ssp_info=None, bathymetry_info=None, options=None ):
    # Writes BELLHOP env file
    #
    # SYNTAX: wbellhopenvfil(filename, thetitle, source_info, surface_info, ssp_info, bottom_info, options)

    #*******************************************************************************
    # Faro, Qua Mai 26 19:06:06 WEST 2021
    # Written by Orlando Camargo Rodriguez
    #*******************************************************************************
    
    envfil = filename + '.env'
    atifil = filename + '.ati' 
    btyfil = filename + '.bty'
    sspfil = filename + '.ssp' 
    
    #*******************************************************************************
    # Get source data: 

    box      = source_info["box"]
    zs       = source_info["zs"]
    freq     = source_info["f"]
    thetas   = source_info["thetas"]
    beamp    = source_info["p"]
    beamc    = source_info["comp"]

    nthetas = thetas[0]
    theta1  = thetas[1]
    theta2  = thetas[2]

    #*******************************************************************************
    # Get surface data: 

    sitype = surface_info["itype"]
    xati   = surface_info["x"]
    aunits = surface_info["units"]
    pati   = surface_info["p"]
 
    #*******************************************************************************
    # Get sound speed data:
   
    c = ssp_info["c"]
    r = ssp_info["r"]     
    z = ssp_info["z"]
   
    csize = c.shape
   
    Dmax = z[-1]
    if len( r ) > 0:
       nr = r.size
    else:
       nr = 0   
   
    #*******************************************************************************
    # Get bathymetry data:

    bitype = bathymetry_info["itype"]
    bunits = bathymetry_info["units"]
    xbty   = bathymetry_info["x"]
    pbty   = bathymetry_info["p"]

    #*******************************************************************************  
    # Get options: 
    options1 = options["options1"]
    options2 = options["options2"]
    options3 = options["options3"]
    options4 = options["options4"]
    rarray   = options["rarray"]
    zarray   = options["zarray"]
    
    #*******************************************************************************
    # Write the ENVFIL:
    fid = open(envfil, 'w')
    fid.write('\'');fid.write(thetitle);fid.write('\'\n')
    fid.write(str(freq))
    fid.write("\n")
    fid.write('1\n')
    fid.write('\'');fid.write(options1);fid.write('\'\n')
    
    s = options1[0]

    if s == "Q":
       fidssp = open(sspfil, 'w')
       fidssp.write(str(nr))
       fidssp.write("\n")
       fidssp.write(str(r))
       fidssp.write("\n")
       for i in range(c.shape[0]):
           for j in range(c.shape[1]):
               fidssp.write(str(c[i,j]))
           fidssp.write("\n")  
       fidssp.close()

    s = options1[1]

    if s == "A":
       fid.write( str( pati[0] ) );fid.write(" ")
       fid.write( str( pati[1] ) );fid.write(" ")
       fid.write( str( pati[2] ) );fid.write(" ")
       fid.write( str( pati[3] ) );fid.write(" ")
       fid.write( str( pati[4] ) );
       fid.write("/ \n")
       
    if len( options1 ) > 5:
       option5 = options1[5]
    else:
       option5 = " "

    if option5 == "*":
       nati = xati[0,:].size     
       fidati = open(atifil, 'w')
       fidati.write("\'");fidati.write(aitype);fidati.write("\'\n")
       fidati.write(str(nati));fidati.write("\n")
       for i in range(nati):
           fidati.write(str(xati[0,i]));fidati.write(" ")
           fidati.write(str(xati[1,i]));fidati.write("\n")
       fidati.close()
    elif option5 == "~":
       nati = xati[0,:].size
       fidati = open(atifil, 'w')
       fidati.write("\'");fidati.write(bitype);fidati.write("\'\n")
       fidati.write(str(nati));fidati.write("\n")
       for i in range(nati):
           fidati.write(str(xati[0,i]));fidati.write(" ")
           fidati.write(str(xati[1,i]));fidati.write(" ")
           fidati.write(str(xati[2,i]));fidati.write(" ")
           fidati.write(str(xati[3,i]));fidati.write(" ")
           fidati.write(str(xati[4,i]));fidati.write(" ")
           fidati.write(str(xati[5,i]));fidati.write(" ")
           fidati.write(str(xati[6,i]));fidati.write("\n")	   	   
       fidati.close()
    else:
       print("No ati file needed...")    

    fid.write('51 0.0 ' ) # Dummy parameters
    fid.write(str(Dmax) )
    fid.write("\n")
    
    nz = z.size
    
    if nr == 0:
      for i in range(nz):
         fid.write( str( z[i] ) ); fid.write(" ") 
         fid.write( str( c[i] ) ); fid.write(" /\n")
    else:
      for i in range(nz):
         fid.write( str( z[i]   ) ); fid.write(" ") 
         fid.write( str( c[i,0] ) ); fid.write("\n")
    
    fid.write('\'');fid.write( options2 );fid.write('\'')
    fid.write(' 0.0 ' )
    fid.write("\n")

    s = options2[0]

    if s == "A":
       fid.write( str(Dmax) ); fid.write(" ")
       fid.write( str( pbty[0] ) );fid.write(" ")
       fid.write( str( pbty[1] ) );fid.write(" ")
       fid.write( str( pbty[2] ) );fid.write(" ")
       fid.write( str( pbty[3] ) );fid.write(" ")
       fid.write( str( pbty[4] ) );    
       fid.write(" / \n")
      
    s = options2[1]

    if s == "*": 
       nbty = xbty[0,:].size      
       fidbty = open(btyfil, 'w')
       fidbty.write("\'");fidbty.write(bitype);fidbty.write("\'\n")
       fidbty.write(str(nbty));fidbty.write("\n")
       for i in range(nbty):
           fidbty.write(str(xbty[0,i]));fidbty.write(" ")
           fidbty.write(str(xbty[1,i]));fidbty.write("\n")
       fidbty.close()
    elif s == "~":
       nbty = xbty[0,:].size
       fidbty = open(btyfil, 'w')
       fidbty.write("\'");fidbty.write(bitype);fidbty.write("\'\n")
       fidbty.write(str(nbty));fidbty.write("\n")
       for i in range(nbty):
           fidbty.write(str(xbty[0,i]));fidbty.write(" ")
           fidbty.write(str(xbty[1,i]));fidbty.write(" ")
           fidbty.write(str(xbty[2,i]));fidbty.write(" ")
           fidbty.write(str(xbty[3,i]));fidbty.write(" ")
           fidbty.write(str(xbty[4,i]));fidbty.write(" ")
           fidbty.write(str(xbty[5,i]));fidbty.write(" ")
           fidbty.write(str(xbty[6,i]));fidbty.write("\n")	   	   
       fidbty.close()
    else:
       print("No bty file needed...")    
    
    nzs = zs.size
    fid.write( str( nzs ) )
    fid.write("\n")
   
    if nzs == 1:
      fid.write( str( zs[0] ) ); fid.write(" / \n")
    else:
      fid.write( str( zs[0]  ) )
      fid.write( str( zs[-1] ) )
      fid.write(" / \n")  

    nra = rarray.size
    nza = zarray.size

    fid.write( str( nza ) )
    fid.write("\n")

    if nza > 1: 
       fid.write( str( zarray[0]  ) );fid.write(" ")
       fid.write( str( zarray[-1] ) )
       fid.write(" / \n")
    else: 
       fid.write( str( zarray[0] ) )
       fid.write(" / \n")
      
    fid.write( str(nra) )
    fid.write("\n")
    
    if nra > 1: 
       fid.write( str( rarray[0]  ) );fid.write(" ")
       fid.write( str( rarray[-1] ) )
       fid.write(" / \n")
    else: 
       fid.write( str( rarray[0] ) )
       fid.write(" / \n")

    fid.write( options3 ); fid.write("\n")
    fid.write( str( int( nthetas) ) ); fid.write("\n")
    fid.write( str(theta1)  ); fid.write(" ")
    fid.write( str(theta2)  ); fid.write(" / \n")
    fid.write( str(box[0])  ); fid.write(" ")
    fid.write( str(box[1])  ); fid.write(" ")    
    fid.write( str(box[2])  ); fid.write("\n")
    
    if len( options4 ) > 0:
       fid.write( options4 ); fid.write("\n")
       fid.write( beamp[0:2] ); fid.write("\n")
       fid.write( beamp[3:4] )
       fid.write( beamc )
      
    fid.write("\n")      

    fid.close()
