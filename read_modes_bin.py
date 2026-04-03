import numpy as np
import struct

# Persistent variables (using function attributes)
_read_modes_bin_state = {
    'fid': None,
    'iRecProfile': 1,
    'lrecl': None,
    'filename': None
}

def read_modes_bin(filename, freq=None, modes=None):
    """
    Read the modes produced by KRAKEN (binary format)
    
    Parameters:
    -----------
    filename : str
        Filename without extension (assumed to be .moA)
    freq : float, optional
        Frequency (involved for broadband runs, use 0 for single frequency)
    modes : array_like, optional
        Vector of mode indices
    
    Returns:
    --------
    Modes : dict
        Dictionary containing mode data
    """
    global _read_modes_bin_state
    
    Modes = {}
    
    # Open file if not already open or filename changed
    if _read_modes_bin_state['fid'] is None or _read_modes_bin_state['filename'] != filename:
        if _read_modes_bin_state['fid'] is not None:
            _read_modes_bin_state['fid'].close()
        
        _read_modes_bin_state['fid'] = open(filename, 'rb')
        _read_modes_bin_state['filename'] = filename
        _read_modes_bin_state['iRecProfile'] = 1
        
        # Read record length
        recl_word = struct.unpack('i', _read_modes_bin_state['fid'].read(4))[0]
        _read_modes_bin_state['lrecl'] = 4 * recl_word  # convert to bytes
    
    fid = _read_modes_bin_state['fid']
    iRecProfile = _read_modes_bin_state['iRecProfile']
    lrecl = _read_modes_bin_state['lrecl']
    
    # Seek to record
    rec = iRecProfile - 1
    fid.seek(rec * lrecl + 4, 0)
    
    # Read header
    Modes['title'] = struct.unpack('80s', fid.read(80))[0].decode().strip('\x00')
    Modes['Nfreq'] = struct.unpack('i', fid.read(4))[0]
    Modes['Nmedia'] = struct.unpack('i', fid.read(4))[0]
    Ntot = struct.unpack('i', fid.read(4))[0]
    NMat = struct.unpack('i', fid.read(4))[0]
    
    if Ntot < 0:
        return Modes
    
    # Read N and Mater
    rec = iRecProfile
    fid.seek(rec * lrecl, 0)
    
    Modes['N'] = np.zeros(Modes['Nmedia'], dtype=int)
    Modes['Mater'] = []
    
    for medium in range(Modes['Nmedia']):
        Modes['N'][medium] = struct.unpack('i', fid.read(4))[0]
        Modes['Mater'].append(struct.unpack('8s', fid.read(8))[0].decode().strip('\x00'))
    
    # Read depth and density
    rec = iRecProfile + 1
    fid.seek(rec * lrecl, 0)
    
    bulk = np.frombuffer(fid.read(8 * Modes['Nmedia']), dtype=np.float32).reshape(2, Modes['Nmedia'], order='F')
    Modes['depth'] = bulk[0, :]
    Modes['rho'] = bulk[1, :]
    
    # Read frequencies
    rec = iRecProfile + 2
    fid.seek(rec * lrecl, 0)
    
    Modes['freqVec'] = np.frombuffer(fid.read(8 * Modes['Nfreq']), dtype=np.float64)
    
    # Read z
    rec = iRecProfile + 3
    fid.seek(rec * lrecl, 0)
    
    Modes['z'] = np.frombuffer(fid.read(4 * Ntot), dtype=np.float32)
    
    # Find frequency index
    if freq is None:
        freq = 0
    freqdiff = np.abs(Modes['freqVec'] - freq)
    freq_index = np.argmin(freqdiff)
    
    # Skip to selected frequency
    iRecProfile += 4
    rec = iRecProfile
    
    for ifreq in range(freq_index):
        fid.seek(rec * lrecl, 0)
        M_modes = struct.unpack('i', fid.read(4))[0]
        
        if ifreq < freq_index - 1:
            iRecProfile += 3 + M_modes + (4 * (2 * M_modes - 1)) // lrecl
            rec = iRecProfile
    
    # Read number of modes
    fid.seek(rec * lrecl, 0)
    Modes['M'] = struct.unpack('i', fid.read(4))[0]
    
    # Select modes
    if modes is None:
        modes = np.arange(1, Modes['M'] + 1)
    else:
        modes = np.array(modes)
    
    modes = modes[modes <= Modes['M']]
    
    # Read top halfspace info
    rec = iRecProfile + 1
    fid.seek(rec * lrecl, 0)
    
    Modes['Top'] = {}
    Modes['Top']['BC'] = struct.unpack('c', fid.read(1))[0].decode()
    
    cp = np.frombuffer(fid.read(8), dtype=np.float32)
    Modes['Top']['cp'] = complex(cp[0], cp[1])
    
    cs = np.frombuffer(fid.read(8), dtype=np.float32)
    Modes['Top']['cs'] = complex(cs[0], cs[1])
    
    Modes['Top']['rho'] = struct.unpack('f', fid.read(4))[0]
    Modes['Top']['depth'] = struct.unpack('f', fid.read(4))[0]
    
    # Read bottom halfspace info
    Modes['Bot'] = {}
    Modes['Bot']['BC'] = struct.unpack('c', fid.read(1))[0].decode()
    
    cp = np.frombuffer(fid.read(8), dtype=np.float32)
    Modes['Bot']['cp'] = complex(cp[0], cp[1])
    
    cs = np.frombuffer(fid.read(8), dtype=np.float32)
    Modes['Bot']['cs'] = complex(cs[0], cs[1])
    
    Modes['Bot']['rho'] = struct.unpack('f', fid.read(4))[0]
    Modes['Bot']['depth'] = struct.unpack('f', fid.read(4))[0]
    
    # Read modes (eigenfunctions and eigenvalues)
    if Modes['M'] == 0:
        Modes['phi'] = np.array([])
        Modes['k'] = np.array([])
    else:
        Modes['phi'] = np.zeros((NMat, len(modes)), dtype=np.complex64)
        
        # Read eigenfunctions
        for ii, mode_idx in enumerate(modes):
            rec = iRecProfile + 1 + mode_idx
            fid.seek(rec * lrecl, 0)
            
            phi_data = np.frombuffer(fid.read(8 * NMat), dtype=np.float32).reshape(2, NMat, order='F')
            Modes['phi'][:, ii] = phi_data[0, :] + 1j * phi_data[1, :]
        
        # Read wavenumbers
        rec = iRecProfile + 2 + Modes['M']
        fid.seek(rec * lrecl, 0)
        
        k_data = np.frombuffer(fid.read(8 * Modes['M']), dtype=np.float32).reshape(2, Modes['M'], order='F')
        Modes['k'] = (k_data[0, :] + 1j * k_data[1, :]).T
        Modes['k'] = Modes['k'][modes - 1]
    
    # Update state
    _read_modes_bin_state['iRecProfile'] = iRecProfile + 4 + Modes['M'] + (4 * (2 * Modes['M'] - 1)) // lrecl
    
    return Modes

def close_modes_bin():
    """Close the binary mode file"""
    global _read_modes_bin_state
    if _read_modes_bin_state['fid'] is not None:
        _read_modes_bin_state['fid'].close()
        _read_modes_bin_state['fid'] = None
        _read_modes_bin_state['filename'] = None
        _read_modes_bin_state['iRecProfile'] = 1
