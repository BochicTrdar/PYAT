import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.colors import ListedColormap
import scipy.io as sio

def plotmovie(filename):
    """
    Plot a single TL surface in dB or pressure field
    
    Parameters:
    filename : str
        Path to .mat file without extension (or with extension)
    """
    
    dbflag = False  # Set to True for dB scale, False for pressure field
    
    # Load the .mat file (adjust based on your read_shd function's output format)
    # This assumes the .mat file contains variables: PlotTitle, freqVec, Pos, p
    mat_data = sio.loadmat(filename)
    
    # Extract data - adjust these field names based on your actual .mat structure
    # This is a placeholder - you'll need to adapt to your actual data format
    PlotTitle = mat_data.get('PlotTitle', ['Unknown'])[0]
    freqVec   = mat_data.get('freqVec', np.array([0])).flatten()
    Pos       = mat_data.get('Pos', {})
    p         = mat_data.get('pressure', np.array([]))
    # Squeeze p to remove singleton dimensions
    p = np.squeeze(p)
    Pos_struct = Pos[0, 0]
    receiver_info = Pos_struct['r']
    zarray = np.squeeze( receiver_info[0, 0][0] )
    rarray = np.squeeze( receiver_info[0, 0][1] )
    R, Z = np.meshgrid(rarray, zarray)
    
    if p.ndim == 3:
        tlt_first = np.real(p[0, :, :])
    else:
        tlt_first = np.real(p)
    
    if not dbflag:
        tlt_first = 1e6 * tlt_first
        
    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(111)   
    cf = ax.pcolormesh(R, Z, tlt_first, shading='nearest', cmap=plt.cm.jet_r)
    cf.set_edgecolor('face')
    
    ax.set_xlabel('Range (m)')
    ax.set_ylabel('Depth (m)')
    ax.set_xlim([0, 300])
    ax.invert_yaxis()    

    cbar = plt.colorbar(cf, ax=ax)
    Nfreq = len(freqVec)
        
    # Create figure
    #fig = plt.figure(figsize=(10, 6))

    # Loop over frequencies
    for ifreq in range(Nfreq):
        # Extract pressure field for current frequency
        if p.ndim == 3:
            tlt = np.real(p[ifreq, :, :])
        else:
            tlt = np.real(p)
        
        if dbflag:
            # Convert to dB
            tlt = -20 * np.log10(np.abs(tlt) + 1e-16)  # Add small epsilon to avoid log(0)
        else:
            # Scale pressure values (pcolor has issues with very small values)
            tlt = 1e6 * tlt
            
            # Calculate max for colorbar limits
            tlmax = np.max(np.abs(tlt))
            tlmax = 0.2 * max(tlmax, 1e-6)
        # Create pcolormesh plot (equivalent to pcolor with shading interp)
        # Meshgrid for range and depth
        #R, Z = np.meshgrid(Pos_r_r, Pos_r_z)
        R, Z = np.meshgrid(rarray, zarray)        
        
        #cf = ax.pcolormesh(R, Z, tlt, shading='nearest', cmap=ListedColormap(tej))
        cf = ax.pcolormesh(R, Z, tlt, shading='nearest', cmap=plt.cm.jet_r)        
        cf.set_edgecolor('face')
        
        ax.set_xlabel('Range (m)')
        ax.set_ylabel('Depth (m)')
        ax.set_xlim([0, 300])
        ax.set_title(f'{PlotTitle}\nTime = {freqVec[ifreq]:g} s')
        
        # Reverse y-axis (depth increases downward)
        ax.invert_yaxis()
        
        # Add rectangles for air/water layers (optional - uncomment to enable)
        if not dbflag:
            # Example: add semi-transparent rectangles
            # Air layer (top)
            skyblue = [0.7, 0.9, 1.0]
            r_air = [0, 1000, 1000, 0, 0]
            z_air = [-100, -100, 0, 0, -100]
            ax.fill(r_air, z_air, skyblue, alpha=0.2)
            
            # Water layer (bottom)
            seablue = [0.5, 0.5, 1.0]
            r_water = [0, 1000, 1000, 0, 0]
            z_water = [200, 200, 0, 0, 200]
            ax.fill(r_water, z_water, seablue, alpha=0.2)
        
        # Add colorbar
        #plt.colorbar(cf, ax=ax)
        
        # Update display
        plt.draw()
        plt.pause(0.1)  # Small pause to allow GUI to update
        
        # Optional: clear axes for next frame (uncomment if you want to clear each frame)
        if ifreq < Nfreq-1:
           ax.cla()
    
    plt.show()
