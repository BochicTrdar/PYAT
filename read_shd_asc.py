import numpy as np

def read_shd_asc(filename):
    """
    Read an ascii shade file
    
    Parameters:
    -----------
    filename : str
        Shade filename
    
    Returns:
    --------
    PlotTitle : str
        Plot title
    PlotType : str
        Plot type
    freqVec : array
        Frequency vector
    freq0 : float
        Reference frequency
    atten : float
        Attenuation
    Pos : dict
        Position information
    pressure : array
        Pressure field
    """
    
    with open(filename, 'r') as fid:
        PlotTitle = fid.readline().strip()
        PlotType = fid.readline().strip()
        
        # Read dimensions
        line = fid.readline()
        while line.strip() == '':
            line = fid.readline()
        Nfreq = int(line.split()[0])
        Ntheta = int(fid.readline().split()[0])
        Nsd = int(fid.readline().split()[0])
        Nrd = int(fid.readline().split()[0])
        Nrr = int(fid.readline().split()[0])
        
        freq0 = float(fid.readline().split()[0])
        atten = float(fid.readline().split()[0])
        
        # Read vectors
        freqVec = np.array(list(map(float, fid.readline().split())))
        Pos = {}
        Pos['theta'] = np.array(list(map(float, fid.readline().split())))
        Pos['s'] = {'z': np.array(list(map(float, fid.readline().split())))}
        Pos['r'] = {
            'z': np.array(list(map(float, fid.readline().split()))),
            'r': np.array(list(map(float, fid.readline().split())))
        }
        
        # Read pressure data
        # Note: This assumes isd = 1 as in the MATLAB code
        isd = 1
        for ii in range(isd):
            # Read data for one source depth
            data = []
            for _ in range(Nrd):
                line = fid.readline()
                while len(line.split()) < 2 * Nrr:
                    line += ' ' + fid.readline().strip()
                data.append(list(map(float, line.split())))
            
            temp1 = np.array(data)  # shape (Nrd, 2*Nrr)
        
        # Combine real and imaginary parts
        pressure = temp1[:, 0::2] + 1j * temp1[:, 1::2]
    
    return PlotTitle, PlotType, freqVec, freq0, atten, Pos, pressure
