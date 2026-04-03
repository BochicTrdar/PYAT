import numpy as np

def read_modes_asc(filename, modes=None):
    """
    Read the modes produced by KRAKEN (ASCII format)
    
    Parameters:
    -----------
    filename : str
        Filename including extension
    modes : array_like, optional
        Vector of mode indices
    
    Returns:
    --------
    Modes : dict
        Dictionary containing mode data
    """
    
    Modes = {}
    
    with open(filename, 'r') as fid:
        # Read header
        lrecl = int(fid.readline().strip())
        Modes['pltitl'] = fid.readline().strip()
        
        # Read parameters
        temp = list(map(float, fid.readline().split()))
        Modes['freq'] = temp[0]
        Modes['Nmedia'] = int(temp[1])
        Modes['ntot'] = int(temp[2])
        Modes['nmat'] = int(temp[3])
        Modes['M'] = int(temp[4])
        
        # Skip media info
        for ii in range(Modes['Nmedia']):
            fid.readline()
        
        # Skip halfspace properties
        fid.readline()  # top halfspace
        fid.readline()  # bot halfspace
        fid.readline()  # blank line
        
        # Read depths
        line = ''
        while len(line.split()) < Modes['ntot']:
            line += fid.readline()
        Modes['z'] = np.array(list(map(float, line.split())))
        
        # Read wavenumbers
        line = ''
        while len(line.split()) < 2 * Modes['M']:
            line += fid.readline()
        ckt = np.array(list(map(float, line.split()))).reshape(2, Modes['M'], order='F')
        Modes['k'] = ckt[0, :] + 1j * ckt[1, :]
        
        # Select modes
        if modes is None:
            modes = np.arange(1, Modes['M'] + 1)
        else:
            modes = np.array(modes)
        
        # Keep only modes that exist
        modes = modes[modes <= Modes['M']]
        Modes['k'] = Modes['k'][modes - 1]
        
        # Initialize mode shape array
        Modes['phi'] = np.zeros((Modes['ntot'], len(modes)), dtype=complex)
        
        # Read modes
        mode_list = list(range(1, Modes['M'] + 1))
        for mode in mode_list:
            # Skip mode header line
            fid.readline()
            
            # Read mode data
            line = ''
            while len(line.split()) < 2 * Modes['ntot']:
                line += fid.readline()
            phit = np.array(list(map(float, line.split()))).reshape(2, Modes['ntot'], order='F')
            
            # Store if in list
            if mode in modes:
                idx = np.where(modes == mode)[0][0]
                Modes['phi'][:, idx] = phit[0, :] + 1j * phit[1, :]
    
    return Modes
