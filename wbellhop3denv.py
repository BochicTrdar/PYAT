def wbellhop3denv(filename=None, thetitle=None, source_info=None, surface_info=None, ssp_info=None, bathymetry_info=None, options=None ):
    # Writes BELLHOP3D env file
    #
    # SYNTAX: wbellhop3denv(filename, thetitle, source_info, surface_info, ssp_info, bottom_info, options)

    #*******************************************************************************
    # Faro, Seg 11 Abr 2022 13:38:52 WEST 
    # Written by Orlando Camargo Rodriguez
    #*******************************************************************************
    
    envfil = filename + '.env'
    atifil = filename + '.ati' 
    btyfil = filename + '.bty'
    sspfil = filename + '.ssp' 
    
    #*******************************************************************************
    # Get source data: 

    box      = source_info["box"]
    xs       = source_info["xs"]    
    ys       = source_info["ys"]
    zs       = source_info["zs"]
    freq     = source_info["f"]
    thetas   = source_info["thetas"]
    phi      = source_info["phi"]
    bearings = source_info["bearings"]   
    beamp    = source_info["p"]
    beamc    = source_info["comp"]
    
    #*******************************************************************************
    # Get surface data: 

    sitype = surface_info["itype"]
    xati   = surface_info["xati"]
    yati   = surface_info["yati"]
    zati   = surface_info["yati"]        
    aunits = surface_info["units"]
    pati   = surface_info["p"]
 
    #*******************************************************************************
    # Get sound speed data:
   
    c = ssp_info["c"]
    r = ssp_info["r"]     
    z = ssp_info["z"]
   
    Dmax = z[-1]
   
    #*******************************************************************************
    # Get bathymetry data:

    bitype = bathymetry_info["itype"]
    bunits = bathymetry_info["units"]
    xbty   = bathymetry_info["xbty"]
    ybty   = bathymetry_info["ybty"]
    zbty   = bathymetry_info["zbty"]    
    pbty   = bathymetry_info["p"]

    #*******************************************************************************  
    # Get options: 
    options1 = options["options1"]
    options2 = options["options2"]
    options3 = options["options3"]
    options4 = options["options4"]  
    opnrd    = options["nrd"]
    opnr     = options["nr"]
    oprd     = options["rd"]
    opr      = options["r"]

    #*******************************************************************************
    # Write the ENVFIL:
    fid = open(envfil, 'w')
    fid.write('\'');fid.write(thetitle);fid.write('\'\n')
    fid.write(str(freq))
    fid.write("\n")
    fid.write('1\n') # Dummy parameter
    fid.write('\'');fid.write(options1);fid.write('\'\n')

    s = options1[2]
    
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
       
    if option5 != " ": 
       nxati = xati.size
       nyati = yati.size      
       fidati = open(atifil, 'w')
       fidati.write('\'')
       fidati.write(sitype)
       fidati.write('\'\n')
       fidati.write(str(nxati));fidati.write("\n")
       for i in range(nxati):
           fidati.write(str(xati[i]));fidati.write(" ")
       fidati.write("\n")
       fidati.write(str(nyati));fidati.write("\n")
       for i in range(nyati):
           fidati.write(str(yati[i]));fidati.write(" ")
       fidati.write("\n")
       for i in range(nyati):
           for j in range(nxati):
               fidati.write(str(zati[j,i]));fidati.write(" ")
           fidati.write("\n")
       fidati.close()

    fid.write('51 0.0 ' ) # Dummy parameters
    fid.write(str(Dmax) )
    fid.write("\n")
    
    nr = r.size
    nz = z.size
    
    if nr == 1:
      for i in range(nz):
         fid.write( str( z[i] ) ); fid.write(" ") 
         fid.write( str( c[i] ) ); fid.write(" /\n")
	 
    fid.write('\'')
    fid.write(options2)
    fid.write('\'')
    fid.write(" 0.0\n") # Roughness (not considered)
    
    s = options2[0]

    if s == "A":
       fid.write( str(Dmax) ); fid.write(" ")
       fid.write( str( pbty[0] ) );fid.write(" ")
       fid.write( str( pbty[1] ) );fid.write(" ")
       fid.write( str( pbty[2] ) );fid.write(" ")
       fid.write( str( pbty[3] ) );fid.write(" ")    
       fid.write(" / \n")
      
    s = options2[1]

    if s == "*":        
       nxbty = xbty.size
       nybty = ybty.size      
       fidbty = open(btyfil, 'w')
       fidbty.write('\'')
       fidbty.write(bitype)
       fidbty.write('\'\n')       
       fidbty.write(str(nxbty));fidbty.write("\n")
       for i in range(nxbty):
           fidbty.write(str(xbty[i]));fidbty.write(" ")
       fidbty.write("\n")
       fidbty.write(str(nybty));fidbty.write("\n")
       for i in range(nybty):
           fidbty.write(str(ybty[i]));fidbty.write(" ")
       fidbty.write("\n")
       for i in range(nybty):
           for j in range(nxbty):
               fidbty.write(str(zbty[j,i]));fidbty.write(" ")
           fidbty.write("\n")
       fidbty.close()

    nxs = xs.size
    fid.write( str( nxs ) );fid.write("\n")
    if nxs == 1:
       fid.write( str( xs[0] ) );fid.write(" /\n")
    else:
       fid.write( str( xs[0] ) );fid.write(" ") 
       fid.write( str( xs[1] ) );fid.write(" /\n")

    nys = ys.size
    fid.write( str( nys ) );fid.write("\n")
    if nys == 1:
       fid.write( str( ys[0] ) );fid.write(" /\n")
    else:
       fid.write( str( ys[0] ) );fid.write(" ") 
       fid.write( str( ys[1] ) );fid.write(" /\n")
    
    nzs = zs.size
    fid.write( str( nzs ) );fid.write("\n")
    if nzs == 1:
       fid.write( str( zs[0] ) );fid.write(" /\n")
    else:
       fid.write( str( zs[0] ) );fid.write(" ") 
       fid.write( str( zs[1] ) );fid.write(" /\n")
    
    fid.write( str(opnrd) ); fid.write("\n")
    if opnrd == 1:
       fid.write( str( oprd[0] ) );fid.write(" /\n")
    else:
       fid.write( str( oprd[0] ) );fid.write(" ")    
       fid.write( str( oprd[1] ) );fid.write(" /\n")
    
    fid.write( str(opnr) ); fid.write("\n")
    if opnr == 1:
       fid.write( str( opr[0] ) );fid.write(" /\n")
    else:
       fid.write( str( opr[0] ) );fid.write(" ")    
       fid.write( str( opr[1] ) );fid.write(" /\n")
          
    fid.write( str(int(bearings[0])) ); fid.write("\n")
    fid.write( str(bearings[1]) ); fid.write(" ")
    fid.write( str(bearings[2]) ); fid.write(' /\n')
  
    fid.write("\'");fid.write(options3);fid.write("\'\n")
    
    fid.write( str( int( thetas[0]) ) );fid.write(' ')
    fid.write( str( int( thetas[1]) ) );fid.write('\n')
    fid.write( str( thetas[2] ) );fid.write(' ')
    fid.write( str( thetas[3] ) );fid.write(' /\n')

    fid.write( str( int( phi[0]) ) );fid.write(' ')
    fid.write( str( int( phi[1]) ) );fid.write('\n')
    fid.write( str( phi[2] ) );fid.write(' ')
    fid.write( str( phi[3] ) );fid.write(' /\n')

    fid.write( str( box[0] ) );fid.write(' ')
    fid.write( str( box[1] ) );fid.write(' ')
    fid.write( str( box[2] ) );fid.write(' ')
    fid.write( str( box[3] ) );fid.write('\n')

    fid.close()
