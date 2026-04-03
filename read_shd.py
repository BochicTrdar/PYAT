# DeepSeek + retouching (and needing some fixes with file extensions)
# Faro, qui 19 fev 2026 21:06:57 
import numpy as np
from scipy.io import loadmat
import os
from read_shd_asc import *
from read_shd_bin import *

def read_shd(filename, freq=None, xs=None, ys=None):
    """
    Read the shade file
    
    Parameters:
    -----------
    filename : str
        Filename (recommended to include extension)
    freq : float, optional
        Read source at specified frequency
    xs, ys : float, optional
        Read source at specified x, y coordinate (in km)
    
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
        5-D pressure field p(Nfreq, Ntheta, Nsd, Nrd, Nrr)
    """
    
    PlotType = None
    
    # Parse filename
    filepath, filename_only = os.path.split(filename)
    basename, ext = os.path.splitext(filename_only)

    # Determine file type
    if ext == '.mat':
        basename2, ext2 = os.path.splitext(basename)
        if ext2 == '.shd':
            FileType = 'shdmat'
        elif ext2 == '.grn':
            FileType = 'grnmat'
        else:
            FileType = 'mat'
    else:
        if filename_only == 'ASCFIL':
            FileType = 'asc'
        elif filename_only == basename:
            FileType = 'shd'
            filename = filename + '.shd'
        elif filename_only == 'tl.grid':
            FileType = 'RAM'
        else:
            FileType = ext.lower() if ext else ''
       
    # Read based on file type
    if FileType in ['shd', 'grn']:  # binary format
        if xs is not None and ys is not None:
            PlotTitle, PlotType, freqVec, freq0, atten, Pos, pressure = read_shd_bin(filename, xs, ys)
        elif freq is not None:
            PlotTitle, PlotType, freqVec, freq0, atten, Pos, pressure = read_shd_bin(filename, freq)
        else:
            PlotTitle, PlotType, freqVec, freq0, atten, Pos, pressure = read_shd_bin(filename)
    
    elif FileType == 'shdmat':  # Shade function mat file
        data = loadmat(filename)
        
        # Extract data (simplified - adjust based on actual .mat structure)
        PlotTitle = data.get('PlotTitle', [''])[0]
        if isinstance(PlotTitle, np.ndarray):
            PlotTitle = str(PlotTitle)
        
        PlotType = data.get('PlotType', [''])[0]
        freqVec = data.get('freqVec', np.array([])).flatten()
        freq0 = data.get('freq0', 0).flatten()[0]
        atten = data.get('atten', 0).flatten()[0]
        
        Pos = {}
        if 'Pos' in data:
            Pos_data = data['Pos'][0, 0]
            # Convert MATLAB struct to dict (simplified)
            for field in Pos_data.dtype.names:
                Pos[field] = Pos_data[field][0, 0]
        
        pressure = data.get('pressure', np.array([]))
        
        # Handle specific source xs, ys
        if xs is not None and ys is not None and 'Pos' in Pos and 's' in Pos:
            if 'x' in Pos['s'] and 'y' in Pos['s']:
                xdiff = np.abs(np.array(Pos['s']['x']) - xs * 1000)
                idxX = np.argmin(xdiff)
                ydiff = np.abs(np.array(Pos['s']['y']) - ys * 1000)
                idxY = np.argmin(ydiff)
                
                # Extract appropriate source index
                pressure = pressure[idxX, idxY, ...]
        
        # Handle specific frequency
        if freq is not None:
            freqdiff = np.abs(freqVec - freq)
            ifreq = np.argmin(freqdiff)
            pressure = pressure[ifreq, ...]
    
    elif FileType == 'asc':  # ascii format
        PlotTitle, PlotType, freqVec, freq0, atten, Pos, pressure = read_shd_asc(filename)
    
    elif FileType == 'grnmat':  # Green's function mat file
        data = loadmat(filename)
        
        PlotTitle = data.get('PlotTitle', [''])[0]
        PlotType = data.get('PlotType', [''])[0]
        freqVec = data.get('freqVec', np.array([])).flatten()
        freq0 = data.get('freq0', 0).flatten()[0]
        atten = data.get('atten', 0).flatten()[0]
        
        Pos = {}
        if 'Pos' in data:
            Pos_data = data['Pos'][0, 0]
            for field in Pos_data.dtype.names:
                Pos[field] = Pos_data[field][0, 0]
        
        if 'r' in Pos and 'r' in Pos['r']:
            Pos['r']['r'] = Pos['r']['r'].flatten()[:, np.newaxis]
        
        pressure = data.get('pressure', np.array([]))
        
        # Handle specific frequency
        if freq is not None:
            freqdiff = np.abs(freqVec - freq)
            ifreq = np.argmin(freqdiff)
            pressure = pressure[ifreq, ...]
    
    elif FileType == 'RAM':
        # Placeholder for read_ram_tlgrid
        raise NotImplementedError('RAM format not yet implemented')
    
    else:
        raise ValueError('Unrecognized file extension')
    
    return PlotTitle, PlotType, freqVec, freq0, atten, Pos, pressure
