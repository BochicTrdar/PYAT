"""
fieldsparc.py - Calculate fields from Green's function produced by SPARC 
Translated/retouched from MATLAB fieldsco.m

This module transforms pressure from (time, wavenumber) domain to (time, range) domain
using Hankel/Bessel transforms.
"""

import numpy as np
from scipy.special import hankel1, hankel2
from scipy.interpolate import interp1d
from scipy.io import loadmat, savemat
import warnings
import os
from pathlib import Path

def readvector(fid):
    """
    Reads a row-vector with Fortran-style '/' termination.
    
    Args:
        fid: Open file handle
    
    Returns:
        x: Vector of values
        Nx: Number of values
    """
    # Read Nx, ignoring anything after ! or %
    line = fid.readline()
    while line.strip() == '':
        line = fid.readline()
    
    # Strip comments and whitespace, convert to int
    comment_pos = line.find('!')
    if comment_pos != -1:
        line = line[:comment_pos]
    comment_pos = line.find('%')
    if comment_pos != -1:
        line = line[:comment_pos]
    
    Nx = int(line.strip())
    
    # Read the data line
    tline = fid.readline()
    while tline.strip() == '':
        tline = fid.readline()
    
    if '/' in tline:
        # Extract values before '/'
        data_part = tline.split('/')[0].strip()
        values = data_part.split()
        x_vals = [float(v) for v in values[:2]]
        
        if Nx == 1:
            x = x_vals[0]
        elif Nx > 2:
            if len(x_vals) > 1:
                x = np.linspace(x_vals[0], x_vals[1], Nx)
            else:
                x = x_vals[0] * np.ones(Nx)
        else:
            x = np.array(x_vals[:Nx])
    else:
        # Read all values individually
        x = np.array([float(v) for v in tline.strip().split()])
        # Read more lines if needed
        while len(x) < Nx:
            next_line = fid.readline()
            if not next_line:
                break
            if next_line.strip():
                x = np.concatenate([x, [float(v) for v in next_line.strip().split()]])
        x = x[:Nx]
    
    return x, Nx

def readr(fid_file_handler):
    """
    Read receiver ranges from file.
    
    Args:
        fid_file_handler: Open file handler
    
    Returns:
        Rr: Receiver ranges in meters
    """
    print('\n_______________________')
    
    Rr, NRr = readvector(fid_file_handler)
    
    print(f'\n Number of receiver ranges, NRr = {NRr}')
    print('\n Receiver ranges, Rr (km)')
    
    if NRr < 10:
        print(' '.join([f'{r:8.2f}' for r in Rr]))
    else:
        print(f'{Rr[0]:8.2f} ... {Rr[-1]:8.2f}')
    
    return Rr  # returns in km, will be converted to m in fieldsco


def readpat(sbpfil, SBP):
    """
    Read source beam pattern.
    
    Args:
        sbpfil: Source beam pattern file root name (without extension)
        SBP: Option flag ('*' to read file, 'O' for omni)
    
    Returns:
        SrcBmPat: Source beam pattern as linear amplitude vs angle
    """
    if SBP == '*':
        print('-----------------------------------')
        print('Using source beam pattern file')
        
        sbp_filename = f'{sbpfil}.sbp'
        if not os.path.exists(sbp_filename):
            warnings.warn(f'Source beam pattern file {sbp_filename} not found. Using omni.')
            SBP = 'O'
        else:
            with open(sbp_filename, 'r') as fid:
                NSBPPts_line = fid.readline()
                while NSBPPts_line.strip() == '':
                    NSBPPts_line = fid.readline()
                NSBPPts = int(NSBPPts_line.strip())
                print(f'Number of source beam pattern points = {NSBPPts}')
                
                SrcBmPat = np.zeros((NSBPPts, 2))
                print(' ')
                print(' Angle (degrees)  Power (dB)')
                
                for I in range(NSBPPts):
                    line = fid.readline()
                    while line.strip() == '':
                        line = fid.readline()
                    vals = line.strip().split()
                    SrcBmPat[I, 0] = float(vals[0])
                    SrcBmPat[I, 1] = float(vals[1])
                    print(f' {SrcBmPat[I, 0]:7.2f}         {SrcBmPat[I, 1]:6.2f}')
    
    if SBP != '*':
        # No pattern given, use omni source pattern
        SrcBmPat = np.array([[-180.0, 0.0],
                              [180.0, 0.0]])
    
    # Convert dB to linear scale
    SrcBmPat[:, 1] = 10.0 ** (SrcBmPat[:, 1] / 20.0)
    
    return SrcBmPat


def scalecol(A, mults):
    """
    Scale each column of A by the respective multiplier in mults.
    Equivalent to B = A * diag(mults) without constructing full matrix.
    
    Args:
        A: Input matrix
        mults: Multiplier vector (row or column)
    
    Returns:
        B: Scaled matrix
    """
    mults = np.asarray(mults).flatten()  # Make sure it's 1D
    
    if A.ndim == 1:
        return mults * A
    
    # For 2D arrays
    B = A * mults[np.newaxis, :]  # Broadcast multiplication
    return B


def taper(G, k, Nk, kleft, kright):
    """
    Apply Hanning window tapering to smooth discontinuities at spectrum ends.
    
    Args:
        G: Matrix to taper (Nrd × Nk)
        k: Wavenumber vector
        Nk: Number of wavenumbers
        kleft: Left limit for tapering
        kright: Right limit for tapering
    
    Returns:
        G: Tapered matrix
    """
    # Calculate left window
    if kleft > k[0]:
        Nwinleft = 2 * round((kleft - k[0]) / (k[-1] - k[0]) * Nk) + 1
        Nwinleft = int(Nwinleft)
        if Nwinleft % 2 == 0:
            Nwinleft += 1
        winleft = np.hanning(Nwinleft)
        Nwinleft = (Nwinleft - 1) // 2
        winleft = winleft[:Nwinleft]  # Take just left half
    else:
        Nwinleft = 0
        winleft = np.array([])
    
    # Calculate right window
    if kright < k[-1]:
        Nwinright = 2 * round((k[-1] - kright) / (k[-1] - k[0]) * Nk) + 1
        Nwinright = int(Nwinright)
        if Nwinright % 2 == 0:
            Nwinright += 1
        winright = np.hanning(Nwinright)
        Nwinright = (Nwinright - 1) // 2
        winright = winright[-Nwinright:]  # Take just right half
    else:
        Nwinright = 0
        winright = np.array([])
    
    if Nk < Nwinleft + Nwinright:
        raise ValueError('Phase speed limits for windowing are invalid')
    
    # Create full window
    window = np.concatenate([winleft, np.ones(Nk - Nwinleft - Nwinright), winright])
    
    # Apply window to each row
    G = scalecol(G, window)
    
    return G

def read_shd_bin(filename, freq=None,xs=None,ys=None):
    #print( freq )
    with open(filename, 'rb') as fid:
        # Record 1: Header length and Title
        recl = np.fromfile(fid, dtype=np.int32, count=1)[0]
        record_size = 4 * recl  # bytes per record
        title = fid.read(80).decode('ascii').strip()

        # Record 2: PlotType
        fid.seek(1 * record_size, 0)
        plot_type = fid.read(10).decode('ascii')

        # Record 3: Dimensions
        fid.seek(2 * record_size, 0)
        n_freq  = np.fromfile(fid, dtype=np.int32, count=1)[0]
        n_theta = np.fromfile(fid, dtype=np.int32, count=1)[0]
        n_sx    = np.fromfile(fid, dtype=np.int32, count=1)[0]
        n_sy    = np.fromfile(fid, dtype=np.int32, count=1)[0]
        n_sz    = np.fromfile(fid, dtype=np.int32, count=1)[0]
        n_rz    = np.fromfile(fid, dtype=np.int32, count=1)[0]
        n_rr    = np.fromfile(fid, dtype=np.int32, count=1)[0]
        freq0   = np.fromfile(fid, dtype=np.float64, count=1)[0]
        atten   = np.fromfile(fid, dtype=np.float64, count=1)[0]

        # Record 4: Frequency Vector
        fid.seek(3 * record_size, 0)
        freq_vec = np.fromfile(fid, dtype=np.float64, count=n_freq)
        #print( freq_vec )
        # Record 5: Thetas
        fid.seek(4 * record_size, 0)
        pos_theta = np.fromfile(fid, dtype=np.float64, count=n_theta)

        # Records 6 & 7: Source X and Y
        pos_s_x = []
        pos_s_y = []
        
        if not plot_type.startswith('TL'):
            fid.seek(5 * record_size, 0)
            pos_s_x = np.fromfile(fid, dtype=np.float64, count=n_sx)
            fid.seek(6 * record_size, 0)
            pos_s_y = np.fromfile(fid, dtype=np.float64, count=n_sy)
        else:
            fid.seek(5 * record_size, 0)
            temp_x = np.fromfile(fid, dtype=np.float64, count=2)
            pos_s_x = np.linspace(temp_x[0], temp_x[-1], n_sx)
            fid.seek(6 * record_size, 0)
            temp_y = np.fromfile(fid, dtype=np.float64, count=2)
            pos_s_y = np.linspace(temp_y[0], temp_y[-1], n_sy)

        # Record 8: Source depths (Z)
        fid.seek(7 * record_size, 0)
        pos_s_z = np.fromfile(fid, dtype=np.float32, count=n_sz)

        # Record 9: Receiver depths (Z)
        fid.seek(8 * record_size, 0)
        pos_r_z = np.fromfile(fid, dtype=np.float32, count=n_rz)

        # Record 10: Receiver ranges (R)
        fid.seek(9 * record_size, 0)
        pos_r_r = np.fromfile(fid, dtype=np.float64, count=n_rr)

        pos = {
            'theta': pos_theta,
            's': {'x': pos_s_x, 'y': pos_s_y, 'z': pos_s_z},
            'r': {'z': pos_r_z, 'r': pos_r_r}
        }
        
        # Logic for selection based on PlotType
        if plot_type.strip() == 'irregular':
            n_rcvrs_per_range = 1
        else:
            n_rcvrs_per_range = n_rz

        # Find ifreq index
        ifreq = 0  # Default to first frequency
        if freq is not None:
            ifreq = np.argmin(np.abs(freq_vec - freq))

        # Initialize Pressure Matrix
        pressure = np.zeros((n_theta, n_sz, n_rcvrs_per_range, n_rr), dtype=complex)

        # THE MAIN LOOP (Equivalent to lines 112-136 of .m file)
        for itheta in range(n_theta):
            for isz in range(n_sz):
                for irz in range(n_rcvrs_per_range):
                    
                    # 0-indexed version of MATLAB's recnum calculation
                    recnum = 10 + (ifreq * n_theta * n_sz * n_rcvrs_per_range) + \
                                  (itheta * n_sz * n_rcvrs_per_range) + \
                                  (isz * n_rcvrs_per_range) + \
                                  irz
                    #print( ifreq, recnum )
                    # Equivalent to status = fseek( fid, recnum * 4 * recl, -1 )
                    fid.seek(recnum * record_size, 0)
                    
                    # Read 2 * Nrr float32 values (Real and Imaginary pairs)
                    temp = np.fromfile(fid, dtype=np.float32, count=2 * n_rr)
                    
                    # Reconstruct complex data: real + 1j * imag
                    pressure[itheta, isz, irz, :] = temp[0::2] + 1j * temp[1::2]

    return title, plot_type, freq_vec, freq0, atten, pos, pressure

def fieldsparc(filename, PlotTitle=None, freq=None, atten=None, Pos=None, 
             Gtemp=None, Opt=None, Rminkm=None, Rmaxkm=None, Nrr=None,
             cmin=1e-10, cmax=1e30):
    """
    Calculate field using Green's function from SPARC/SCOOTER.
    
    Transforms pressure from (time, wavenumber) domain to (time, range) domain.
    
    Args:
        filename: Path to .grn.mat file (without extension or with)
        PlotTitle: Optional plot title
        freq: Optional frequency
        atten: Optional attenuation
        Pos: Optional position structure
        Gtemp: Optional Green's function data
        Opt: Optional options string (e.g., 'RPO' for point source, omni)
        Rminkm: Minimum range in km
        Rmaxkm: Maximum range in km
        Nrr: Number of ranges
        cmin: Minimum phase speed for tapering
        cmax: Maximum phase speed for tapering
    
    Returns:
        pressure: Time-domain pressure field (times × source depths × receiver depths × ranges)
        Pos: Updated position structure with ranges
        freqVec: Time vector
    """
    
    fileroot = filename[:-4] if filename.endswith('.grn') else filename
    
    print(f'fieldsparc {filename}')
    
    # If only filename provided, read control info from .flp file
    if PlotTitle is None:
        flp_filename = f'{fileroot}.flp'
        if not os.path.exists(flp_filename):
            raise FileNotFoundError(f'flp file missing: {flp_filename}')
        
        with open(flp_filename, 'r') as fid:
            # Read option line
            Opt_line = fid.readline()
            while Opt_line.strip() == '':
                Opt_line = fid.readline()
            # Strip comments from option line too
            comment_pos = Opt_line.find('!')
            if comment_pos != -1:
               Opt_line = Opt_line[:comment_pos]
            # Extract letters between quotes
            quote_pos = [i for i, c in enumerate(Opt_line) if c == "'"]
            if len(quote_pos) >= 2:
                Opt = Opt_line[quote_pos[0] + 1:quote_pos[1]]
            else:
                Opt = Opt_line.strip()
            
            if len(Opt) <= 3:
                Opt = Opt + 'O'  # Default beam pattern is omni
            
            print(f'\n Fieldsparc option line = {Opt}')
            
            # Read range vector
            Rr_km, _ = readvector(fid)
            Nrr = len(Rr_km)
            Rr = 1000 * Rr_km  # Convert km to m
            Rr[np.abs(Rr) < 1e-12] = 1  # Replace zero range with 1 m
            rarray = Rr
        # Read Green's function metadata
        xs = np.nan
        ys = np.nan
        PlotTitle, PlotType, freqVec, freq0, atten, Pos, _ = read_shd_bin(filename, freq, xs, ys)
    # Set phase speed limits for optional tapering
    # These can be modified by user
    # cmin = 1e-10  # Disable tapering by default
    # cmax = 1e30
    
    # Read source beam pattern
    SBP = Opt[3] if len(Opt) > 3 else 'O'
    SrcBmPat = readpat(fileroot, SBP)
    
    # Get dimensions
    source_info   = Pos['s'] 
    receiver_info = Pos['r']
    zs     =   source_info["z"]
    zarray = receiver_info["z"]
    Nsz   = np.size( zs )        # Number of source   depths
    Nrz   = np.size( zarray )    # Number of receiver depths
    #Nsz = len(Pos.s.z)  # Number of source depths
    #Nrz = len(Pos.r.z)  # Number of receiver depths
    Nfreq = len(freqVec)  # Number of frequencies (times for SPARC)
    
    # Initialize pressure array
    pressure = np.zeros((Nfreq, Nsz, Nrz, Nrr), dtype=complex)
    # Loop over frequencies/times
    for ifreq in range(Nfreq):
        freq_val = freqVec[ifreq]
        omega = 2.0 * np.pi * freq_val
        kleft = omega / cmax  # Left limit for tapering
        kright = omega / cmin  # Right limit for tapering
        
        # Read Green's function for this frequency/time
        #_, _, _, freq0, atten_val, Pos, Gtemp = read_shd(filename, freqVec[ ifreq ], xs, ys )
        _, _, _, freq0, atten_val, Pos, Gtemp = read_shd_bin(filename, freqVec[ifreq], xs, ys)
        # Get phase speeds (stored in Pos.r.r)
        #cVec = Pos.r.r.flatten()
        cVec = receiver_info['r']

        k = 2 * np.pi * freq0 / cVec
        
        Nk = len(k)
        deltak = (k[-1] - k[0]) / (Nk - 1)
        
        # Check wavenumber sampling
        if Rr[-1] * deltak > 10:
            warnings.warn('The wavenumber sampling is too coarse to accurately calculate the field at the largest range')
        
        # Create complex wavenumber
        ck = k + 1j * atten_val
        
        # Prepare transform matrices
        # x is a matrix of size nk × nrr: ck * r
        x = np.outer(ck, np.abs(Rr))  # (Nk, Nrr)
        
        # Select transform type
        transform_type = Opt[0] if len(Opt) > 0 else 'R'
        
        if transform_type == 'X':  # Line source
            X = np.exp(-1j * x)
            X2 = np.exp(1j * x)
            factor1 = 1
            factor2 = deltak / np.sqrt(2 * np.pi)
        elif transform_type == 'R':  # Point source
            X = np.exp(-1j * (x - np.pi / 4))
            X2 = np.exp(1j * (x - np.pi / 4))
            factor1 = np.sqrt(ck)[:, np.newaxis]
            factor2 = deltak / np.sqrt(2 * np.pi * np.abs(Rr))
        elif transform_type == 'S':  # Point source with cylindrical spreading removed
            X = np.exp(-1j * (x - np.pi / 4))
            X2 = np.exp(1j * (x - np.pi / 4))
            factor1 = np.sqrt(ck)[:, np.newaxis]
            factor2 = deltak / np.sqrt(2 * np.pi)
        elif transform_type == 'H':  # Exact Bessel transform
            # Hankel functions - careful with broadcasting
            X = np.zeros_like(x, dtype=complex)
            X2 = np.zeros_like(x, dtype=complex)
            for i in range(x.shape[0]):
                for j in range(x.shape[1]):
                    X[i, j] = hankel1(0, x[i, j])
                    X2[i, j] = hankel2(0, x[i, j])
            factor1 = ck[:, np.newaxis]
            factor2 = deltak / 2
        else:
            print('fieldsparc.py -- Unknown first letter option in second line of .flp file')
            return None
        
        # Loop over source depths
        for isz in range(Nsz):
            # Extract Green's function for this source depth
            if Gtemp.ndim >= 4:
                G = np.squeeze(Gtemp[isz, 0, :, :])
            elif Gtemp.ndim == 3:
                G = np.squeeze(Gtemp[isz, :, :])
            else:
                G = Gtemp
            
            # Ensure G is 2D
            if G.ndim == 1:
                G = G.reshape(1, -1)
            
            # Apply source beam pattern
            if SBP == '*':
                c = 1500  # Reference sound speed
                kz2 = omega**2 / c**2 - k**2
                kz2[kz2 < 0] = 0
                theta = np.degrees(np.arctan2(np.sqrt(kz2), k))
                
                # Interpolate beam pattern
                interp_func = interp1d(SrcBmPat[:, 0], SrcBmPat[:, 1], 
                                       kind='linear', bounds_error=False, 
                                       fill_value=0)
                S = interp_func(theta)
                G = scalecol(G, S)
            
            # Apply tapering
            G = taper(G, k, Nk, kleft, kright)
            G[np.abs(G) < 1e-40] = 0  # Avoid underflows
            
            # Perform transform
            transform_sign = Opt[1] if len(Opt) > 1 else 'P'
            
            # Scale columns of G by factor1
            if factor1.shape[0] == G.shape[1]:
                G_scaled = scalecol(G, factor1.flatten())
            else:
                G_scaled = G
            
            if transform_sign == 'P':  # Positive (outgoing)
                Y = -G_scaled @ X
            elif transform_sign == 'N':  # Negative (incoming)
                Y = -G_scaled @ X2
            elif transform_sign == 'B':  # Both
                Y = -G_scaled @ (X + X2)
            else:
                print('fieldsco.py -- Unknown second letter option in second line of .flp file')
                return None
            
            # Apply cylindrical spreading factor
            if factor2.ndim == 0:
                pressure[ifreq, isz, :, :] = Y * factor2
            elif factor2.shape[0] == Y.shape[1]:
                #pressure[ifreq, isz, :, :] = scalecol(Y.T, factor2).T
                pressure[ifreq, isz, :, :] = scalecol(Y, factor2)
            else:
                pressure[ifreq, isz, :, :] = Y * factor2[np.newaxis, :]
            
            print(f'      Transform completed for source depth {zs[isz]}, frequency (or time) {freq_val}')
    
    # Update position structure with ranges
    #Pos.r.r = Rr
    #atten = 0
    #PlotType = 'rectilin'
    Pos['r']['r'] = rarray
    # Save results
    output_filename = f'{fileroot}.shd.mat'
    savemat(output_filename, {
        'PlotTitle': PlotTitle,
        'PlotType': PlotType,
        'freqVec': freqVec,
        'freq0': freq0,
        'atten': atten,
        'Pos': Pos,
        'pressure': pressure
    })
    
    print(f'\nSaved field to {output_filename}')
    
    return pressure, Pos, freqVec
