def wbounceenvfil(filename=None, thetitle=None, freq=None, nmedia=None, options=None, water_info=None, layers_info=None,
                  clow=None, chigh=None, rmax=None):
    # Writes BOUNCE env file
    #
    # SYNTAX: wbounceenvfil(filename, thetitle, freq, nmedia, options, water_info, layers_info, rmax, clow, chigh)

    #*******************************************************************************
    # Faro, Fri Jun 21 11:27:39 PM WEST 2024
    # Written by Tordar
    #*******************************************************************************
    
    envfil = filename + '.env'
    
    #*******************************************************************************
    # Write the ENVFIL:
    fid = open(envfil, 'w')
    fid.write('\'');fid.write(thetitle);fid.write('\'\n')
    fid.write(str(freq));fid.write('\n')
    fid.write(str(nmedia));fid.write('\n')
    fid.write('\'');fid.write(options);fid.write('\'\n')
    fid.write(str(water_info[0]));fid.write(' ')
    fid.write(str(water_info[1]));fid.write(' ')
    fid.write(str(water_info[2]));fid.write(' ')
    fid.write(str(water_info[3]));fid.write(' ')
    fid.write(str(water_info[4]));fid.write(' ')
    fid.write(str(water_info[5]));fid.write('\n')
    for i in range(nmedia):    
        fid.write(str('0 0.0'));fid.write(' ')
        fid.write(str(layers_info[i,0]));fid.write('\n')
        fid.write(str(water_info[0]));fid.write(' ')
        fid.write(str(layers_info[i,1]));fid.write(' ')
        fid.write(str(layers_info[i,2]));fid.write(' ')
        fid.write(str(layers_info[i,3]));fid.write(' ')
        fid.write(str(layers_info[i,4]));fid.write(' ')
        fid.write(str(layers_info[i,5]));fid.write('\n')
        fid.write(str(layers_info[i,0]));fid.write(' ')
        fid.write(str(layers_info[i,1]));fid.write(' /\n')
    fid.write('\'A\'');fid.write(str(' 0.0'));fid.write('\n')
    fid.write(str(layers_info[nmedia,0]));fid.write(' ')
    fid.write(str(layers_info[nmedia,1]));fid.write(' ')
    fid.write(str(layers_info[nmedia,2]));fid.write(' ')
    fid.write(str(layers_info[nmedia,3]));fid.write(' ')
    fid.write(str(layers_info[nmedia,4]));fid.write(' ')
    fid.write(str(layers_info[nmedia,5]));fid.write('\n')
    fid.write(str(clow));fid.write(' ') 
    fid.write(str(chigh));fid.write('\n')
    fid.write(str(rmax));fid.write('\n')
    fid.close()
