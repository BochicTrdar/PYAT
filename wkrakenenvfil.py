def wkrakenenvfil(filename=None, thetitle=None, source_info=None, surface_info=None, scatter_info=None, ssp_info=None, bottom_info=None, field_info=None ):
    # Writes KRAKEN env file
    #
    # SYNTAX: wkrakenenvfil(filename, thetitle, source_info, surface_info, scatter_info, ssp_info, bottom_info, field_info)

    #*******************************************************************************
    # Faro, Qua Mai 26 19:05:37 WEST 2021
    # Written by Orlando Camargo Rodriguez
    #*******************************************************************************
    
    envfil = filename + '.env'
    brcfil = filename + '.brc'
    trcfil = filename + '.trc'
    ircfil = filename + '.irc'
        
    #*******************************************************************************
    # Get source data: 

    freq = source_info["f"]
    zs   = source_info["zs"]
    
    nzs = zs.size
    
    #*******************************************************************************
    # Get surface data: 

    top_boundary_condition = surface_info["bc"]
    top_properties         = surface_info["properties"]
    top_reflection_coeff   = surface_info["reflection"]
    
    #*******************************************************************************    
    # Get scatter data: 

    bumden = scatter_info["bumden"]
    eta    = scatter_info["eta"]
    xi     = scatter_info["xi"]
    
    #*******************************************************************************
    # Get sound speed data: 

    ssp_data  = ssp_info["cdata"]
    ssp_type  = ssp_info["type"]
    citype    = ssp_info["itype"]
    nmesh     = ssp_info["nmesh"]
    csigma    = ssp_info["sigma"]
    clow      = ssp_info["clow"]
    chigh     = ssp_info["chigh"]
    zbottom   = ssp_info["zbottom"]

    z = ssp_data[0,]; Dmax = max( z ); nz = z.size

    #*******************************************************************************
    # Get bottom data:

    nlayers		      = bottom_info["n"]
    attenuation_units	      = bottom_info["units"]
    bottom_boundary_condition = bottom_info["bc"]
    bottom_properties	      = bottom_info["properties"]
    bsigma		      = bottom_info["sigma"]
    layer_properties	      = bottom_info["layerp"]
    layer_type  	      = bottom_info["layert"]
    layer_data  	      = bottom_info["bdata"]

    if nlayers >= 20:
       print('Max # of layers = 20!!!!!!')
       stop()

    #*******************************************************************************
    # Get field data:

    thorpe	= field_info["thorpe"]
    finder	= field_info["finder"]
    rmax	= field_info["rmax"]
    nrd 	= field_info["nrd"]
    nrr 	= field_info["nrr"]
    rd  	= field_info["rd"]
    rr  	= field_info["rr"]
    nmodes	= field_info["m"]
    source_type = field_info["stype"]
    nprofiles	= field_info["np"]
    rprofiles	= field_info["rp"]
    dr  	= field_info["dr"]
    range_dependent_modes = field_info["rmodes"]

    #*******************************************************************************
    # Construct the options:

    options1 = citype + top_boundary_condition + attenuation_units + thorpe + finder

    options2 = source_type + range_dependent_modes

    #*******************************************************************************  
    # Write the ENVFIL
    fid = open(envfil, 'w')
    fid.write('\'');fid.write(thetitle);fid.write('\'\n')
    fid.write(str(freq))
    fid.write("\n")
    fid.write(str(nlayers))
    fid.write("\n") 
    fid.write('\'');fid.write(options1);fid.write('\'\n')
    
    if top_boundary_condition == 'A':
      fid.write(str(top_properties[0]));fid.write(" ")
      fid.write(str(top_properties[1]));fid.write(" ")
      fid.write(str(top_properties[2]));fid.write(" ")
      fid.write(str(top_properties[3]));fid.write(" ")
      fid.write(str(top_properties[4]));fid.write(" ")
      fid.write(str(top_properties[5]));fid.write(" /")
      fid.write("\n")

    if top_boundary_condition == 'F':
       nthetas = surface_info["nthetas"]
       angle_data = surface_info["angle_data"]  
       fidtrc = fopen(trcfil,'w')
       fidtrc.write(str(nthetas));fid.write("\n")
       fidtrc.write(str(angle_data));fid.write("\n")
       fclose(fidtrc)

    if ( top_boundary_condition == 'F' )|( top_boundary_condition == 'I' ):
       fid.write(str(bumden));fid.write(" ")
       fid.write(str(eta));fid.write(" ")
       fid.write(str(xi));
       fid.write("\n")
       
    fid.write(str(nmesh));fid.write(" ")
    fid.write(str(csigma));fid.write(" ")
    fid.write(str(zbottom))
    fid.write("\n")
    
    if ( citype != 'A' ):

       nz = ssp_data[0,].size
       
       if ( ssp_type == 'H' ):

         fid.write(str(ssp_data[0,0]));fid.write(" ")
         fid.write(str(ssp_data[1,0]));fid.write(" ")
         fid.write(str(ssp_data[2,0]));fid.write(" ")
         fid.write(str(ssp_data[3,0]));fid.write(" ")
         fid.write(str(ssp_data[4,0]));fid.write(" ")
         fid.write(str(ssp_data[5,0]));
         fid.write(" /\n");
	 
         for i in range(nz-1):
             fid.write(str(ssp_data[0,i+1]));fid.write(" ")
             fid.write(str(ssp_data[1,i+1]));fid.write(" /")
             fid.write("\n")
     
       else:
         
         for i in range(nz):
             fid.write(str(ssp_data[0,i]));fid.write(" ")
             fid.write(str(ssp_data[1,i]));fid.write(" ")
             fid.write(str(ssp_data[2,i]));fid.write(" ")
             fid.write(str(ssp_data[3,i]));fid.write(" ")
             fid.write(str(ssp_data[4,i]));fid.write(" ")
             fid.write("\n")

    for i in range(nlayers-1):
        fid.write(str(int(layer_properties[i,0])));fid.write(" ")
        fid.write(str(layer_properties[i,1]));fid.write(" ")
        fid.write(str(layer_properties[i,2]));fid.write("\n")
    
        if ( layer_type[i] == 'H' ):
          fid.write(str(layer_data[i,0,0]));fid.write(" ")
          fid.write(str(layer_data[i,0,1]));fid.write(" ")
          fid.write(str(layer_data[i,0,2]));fid.write(" ")
          fid.write(str(layer_data[i,0,3]));fid.write(" ")
          fid.write(str(layer_data[i,0,4]));fid.write(" ")
          fid.write(str(layer_data[i,0,5]));fid.write("\n")
          fid.write(str(layer_data[i,1,0]));fid.write(" ")
          fid.write(str(layer_data[i,1,1]));fid.write(" /\n")
	    
        else:
          fid.write(str(layer_data[i,0,0]));fid.write(" ")
          fid.write(str(layer_data[i,0,1]));fid.write(" ")
          fid.write(str(layer_data[i,0,2]));fid.write(" ")
          fid.write(str(layer_data[i,0,3]));fid.write(" ")
          fid.write(str(layer_data[i,0,4]));fid.write(" ")
          fid.write(str(layer_data[i,0,5]));fid.write("\n")
          fid.write(str(layer_data[i,1,0]));fid.write(" ")
          fid.write(str(layer_data[i,1,1]));fid.write(" ")
          fid.write(str(layer_data[i,1,2]));fid.write(" ")
          fid.write(str(layer_data[i,1,3]));fid.write(" ")
          fid.write(str(layer_data[i,1,4]));fid.write(" ")
          fid.write(str(layer_data[i,1,5]));fid.write("\n")
    
    fid.write("\'")
    fid.write(bottom_boundary_condition);fid.write("\' ")
    fid.write(str(bsigma))
    fid.write("\n")
    if bottom_boundary_condition == 'A':
      fid.write(str(bottom_properties[0]));fid.write(" ")
      fid.write(str(bottom_properties[1]));fid.write(" ")
      fid.write(str(bottom_properties[2]));fid.write(" ")
      fid.write(str(bottom_properties[3]));fid.write(" ")
      fid.write(str(bottom_properties[4]));fid.write(" ")
      fid.write(str(bottom_properties[5]));fid.write(" /")
      fid.write("\n")
    fid.write(str(clow));fid.write(" ")
    fid.write(str(chigh))
    fid.write("\n")
    fid.write(str(rmax))
    fid.write("\n")
    fid.write(str(nzs))
    fid.write("\n")
    
    if nzs == 1: 
       fid.write(str(zs[0]))
       fid.write(" /\n")
    else:
        fid.write(str(zs[0]));fid.write(" ") 
        fid.write(str(zs[-1]));
        fid.write(" /\n")
	
    fid.write(str(nrd))
    fid.write("\n")
    
    if nrd == 1: 
       fid.write(str(rd[0]))
       fid.write("\n")
    else:
        fid.write(str(rd[0]));fid.write(" ") 
        fid.write(str(rd[-1]));
        fid.write(" /\n")
	
    fid.close()

    #*******************************************************************************
    # Write the FLPFIL:

    fid = open('field.flp','w');
    fid.write(thetitle);
    fid.write("\n")
    fid.write(options2)
    fid.write("\n")    
    fid.write(str(nmodes))
    fid.write("\n")  
    fid.write(str(nprofiles))
    fid.write("\n") 
    fid.write(str(rprofiles))
    fid.write("\n")
    fid.write(str(nrr))
    fid.write("\n")
    
    if nrr == 1:
       fid.write(str(rr[0]))
       fid.write("\n")
    else:
       fid.write(str(rr[0]));fid.write(" ")
       fid.write(str(rr[-1]))
       fid.write(" /\n")
    
    fid.write(str(nzs))
    fid.write("\n")
    
    if nzs == 1: 
       fid.write(str(zs[0]))
       fid.write(" /\n")
    else:
        fid.write(str(zs[0]));fid.write(" ") 
        fid.write(str(zs[-1]));
        fid.write(" /\n")
	
    fid.write(str(nrd))
    fid.write("\n")
    
    if nrd == 1: 
       fid.write(str(rd[0]))
       fid.write(" /\n")
    else:
        fid.write(str(rd[0]));fid.write(" ") 
        fid.write(str(rd[-1]));
        fid.write(" /\n")

# Yes, this is ugly... but it works!!!
	
    fid.write(str(nrd))
    fid.write("\n")
	
    if nrd == 1: 
       fid.write(str(dr[0]))
       fid.write(" /\n")
    else:
        fid.write(str(dr[0]));fid.write(" ") 
        fid.write(str(dr[-1]));
        fid.write(" /\n")
    
    fid.close()	
