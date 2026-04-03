import numpy as np
import struct

def read_shd_bin(filename, *args):
    """
    Read a binary shade file
    
    Parameters:
    -----------
    filename : str
        Shade filename
    *args : optional
        Either:
            freq : float - read at specific frequency
            xs, ys : float - read at specific source coordinates (in km)
    
    Returns:
    --------
    title : str
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
        4-D pressure field p(Ntheta, Nsd, Nrd, Nrr)
    """
    
    # Parse arguments
    if len(args) == 1:
        freq = args[0]
        xs = ys = np.nan
    elif len(args) >= 2:
        xs, ys = args[0], args[1]
        freq = np.nan
    else:
        freq = xs = ys = np.nan
    
    with open(filename, 'rb') as fid:
        # Read record length
        recl = struct.unpack('i', fid.read(4))[0]
        
        # Read title
        title = struct.unpack('80s', fid.read(80))[0].decode().strip('\x00')
        
        # Reposition to end of first record
        fid.seek(4 * recl, 0)
        
        # Read PlotType
        PlotType = struct.unpack('10s', fid.read(10))[0].decode().strip('\x00')
        
        # Reposition to end of second record
        fid.seek(2 * 4 * recl, 0)
        
        # Read dimensions
        Nfreq = struct.unpack('i', fid.read(4))[0]
        Ntheta = struct.unpack('i', fid.read(4))[0]
        Nsx = struct.unpack('i', fid.read(4))[0]
        Nsy = struct.unpack('i', fid.read(4))[0]
        Nsz = struct.unpack('i', fid.read(4))[0]
        Nrz = struct.unpack('i', fid.read(4))[0]
        Nrr = struct.unpack('i', fid.read(4))[0]
        freq0 = struct.unpack('d', fid.read(8))[0]
        atten = struct.unpack('d', fid.read(8))[0]
        
        # Reposition to end of record 3
        fid.seek(3 * 4 * recl, 0)
        
        # Read frequency vector
        freqVec = np.frombuffer(fid.read(8 * Nfreq), dtype=np.float64)
        
        # Reposition to end of record 4
        fid.seek(4 * 4 * recl, 0)
        
        # Read theta
        theta = np.frombuffer(fid.read(8 * Ntheta), dtype=np.float64)
        
        Pos = {'theta': theta}
        
        # Read source positions
        if PlotType[:2] != 'TL':
            # Reposition to end of record 5
            fid.seek(5 * 4 * recl, 0)
            Pos['s'] = {'x': np.frombuffer(fid.read(8 * Nsx), dtype=np.float64)}
            
            # Reposition to end of record 6
            fid.seek(6 * 4 * recl, 0)
            Pos['s']['y'] = np.frombuffer(fid.read(8 * Nsy), dtype=np.float64)
        else:
            # Compressed format for TL from FIELD3D
            # Reposition to end of record 5
            fid.seek(5 * 4 * recl, 0)
            xlim = np.frombuffer(fid.read(16), dtype=np.float64)
            Pos['s'] = {'x': np.linspace(xlim[0], xlim[1], Nsx)}
            
            # Reposition to end of record 6
            fid.seek(6 * 4 * recl, 0)
            ylim = np.frombuffer(fid.read(16), dtype=np.float64)
            Pos['s']['y'] = np.linspace(ylim[0], ylim[1], Nsy)
        
        # Read source depths
        fid.seek(7 * 4 * recl, 0)
        Pos['s']['z'] = np.frombuffer(fid.read(4 * Nsz), dtype=np.float32)
        
        # Read receiver depths
        fid.seek(8 * 4 * recl, 0)
        Pos['r'] = {'z': np.frombuffer(fid.read(4 * Nrz), dtype=np.float32)}
        
        # Read receiver ranges
        fid.seek(9 * 4 * recl, 0)
        Pos['r']['r'] = np.frombuffer(fid.read(8 * Nrr), dtype=np.float64)
        
        # Determine number of receivers per range
        if PlotType == 'irregular':
            Nrcvrs_per_range = 1
        else:
            Nrcvrs_per_range = Nrz
        
        # Initialize pressure array
        pressure = np.zeros((Ntheta, Nsz, Nrcvrs_per_range, Nrr), dtype=complex)
        
        # Read data
        if np.isnan(xs):  # Read first xs, ys, but all theta, sz, and rz
            # Get frequency index
            ifreq = 0
            if not np.isnan(freq):
                freqdiff = np.abs(freqVec - freq)
                ifreq = np.argmin(freqdiff)
            
            for itheta in range(Ntheta):
                for isz in range(Nsz):
                    for irz in range(Nrcvrs_per_range):
                        recnum = 10 + (ifreq) * Ntheta * Nsz * Nrcvrs_per_range + \
                                       (itheta) * Nsz * Nrcvrs_per_range + \
                                       (isz) * Nrcvrs_per_range + irz
                        
                        fid.seek(recnum * 4 * recl, 0)
                        
                        # Read complex data
                        temp = np.frombuffer(fid.read(8 * Nrr), dtype=np.float32)
                        pressure[itheta, isz, irz, :] = temp[0::2] + 1j * temp[1::2]
        
        else:  # Read for specific source x, y
            xdiff = np.abs(Pos['s']['x'] - xs * 1000)
            idxX = np.argmin(xdiff)
            ydiff = np.abs(Pos['s']['y'] - ys * 1000)
            idxY = np.argmin(ydiff)
            
            for itheta in range(Ntheta):
                for isz in range(Nsz):
                    for irz in range(Nrcvrs_per_range):
                        recnum = 10 + (idxX) * Nsy * Ntheta * Nsz * Nrcvrs_per_range + \
                                       (idxY) * Ntheta * Nsz * Nrcvrs_per_range + \
                                       (itheta) * Nsz * Nrcvrs_per_range + \
                                       (isz) * Nrcvrs_per_range + irz
                        
                        fid.seek(recnum * 4 * recl, 0)
                        
                        # Read complex data
                        temp = np.frombuffer(fid.read(8 * Nrr), dtype=np.float32)
                        pressure[itheta, isz, irz, :] = temp[0::2] + 1j * temp[1::2]
    
    return title, PlotType, freqVec, freq0, atten, Pos, pressure
