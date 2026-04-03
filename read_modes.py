import numpy as np
from scipy.io import loadmat
import os
from read_modes_bin import *

def read_modes(modfil, freq=None, modes=None):
    """
    Read the modes produced by KRAKEN
    
    Parameters:
    -----------
    modfil : str
        Filename including extension
    freq : float, optional
        Frequency
    modes : array_like, optional
        Vector of mode indices
    
    Returns:
    --------
    Modes : dict
        Dictionary containing mode data
    """
    
    # Parse filename
    filepath, filename = os.path.split(modfil)
    basename, ext = os.path.splitext(filename)
    
    if not ext:
        ext = '.mod'  # default extension
    else:
        # Handle case where modfil is FileRoot.mod.mat
        if ext == '.mat':
            basename2, ext2 = os.path.splitext(basename)
            if ext2 == '.mod':
                basename = basename2
                ext = '.mod.mat'
    
    modfil_full = os.path.join(filepath, basename + ext)
    
    # Read modal data based on extension
    if ext == '.mod':  # binary format
        if modes is None:
            Modes = read_modes_bin(modfil_full, freq)
        else:
            Modes = read_modes_bin(modfil_full, freq, modes)
    elif ext == '.mod.mat':  # MATLAB format
        data = loadmat(modfil_full)
        Modes = {}
        # Convert MATLAB struct to Python dict (simplified)
        for key in data:
            if not key.startswith('__'):
                Modes[key] = data[key]
    elif ext == '.moa':  # ascii format
        if modes is None:
            Modes = read_modes_asc(modfil_full)
        else:
            Modes = read_modes_asc(modfil_full, modes)
    else:
        raise ValueError('read_modes: Unrecognized file extension')
    
    # Identify index of frequency closest to user-specified value
    if freq is not None:
        freqdiff = np.abs(np.array(Modes['freqVec']) - freq)
        freq_index = np.argmin(freqdiff)
    else:
        freq_index = 0
    
    # Calculate wavenumbers in halfspaces (if there are any modes)
    if Modes.get('M', 0) != 0:
        # Top halfspace
        if Modes['Top']['BC'] == 'A':
            Modes['Top']['k2'] = (2 * np.pi * Modes['freqVec'][0] / Modes['Top']['cp'])**2
            gamma2 = Modes['k']**2 - Modes['Top']['k2']
            Modes['Top']['gamma'] = pekeris_root(gamma2)  # vertical wavenumber
            Modes['Top']['phi'] = Modes['phi'][0, :]  # mode value at halfspace
        else:
            Modes['Top']['rho'] = 1.0
            Modes['Top']['gamma'] = np.zeros_like(Modes['k'])
            Modes['Top']['phi'] = np.zeros(Modes['phi'].shape[1])
        
        # Bottom halfspace
        if Modes['Bot']['BC'] == 'A':
            Modes['Bot']['k2'] = (2 * np.pi * Modes['freqVec'][freq_index] / Modes['Bot']['cp'])**2
            gamma2 = Modes['k']**2 - Modes['Bot']['k2']
            Modes['Bot']['gamma'] = pekeris_root(gamma2)  # vertical wavenumber
            Modes['Bot']['phi'] = Modes['phi'][-1, :]  # mode value at halfspace
        else:
            Modes['Bot']['rho'] = 1.0
            Modes['Bot']['gamma'] = np.zeros_like(Modes['k'])
            Modes['Bot']['phi'] = np.zeros(Modes['phi'].shape[1])
    
    return Modes

def pekeris_root(gamma2):
    """Helper function for Pekeris root calculation"""
    gamma = np.zeros_like(gamma2, dtype=complex)
    pos = gamma2 >= 0
    gamma[pos] = np.sqrt(gamma2[pos])
    gamma[~pos] = 1j * np.sqrt(-gamma2[~pos])
    return gamma
