#==================================================================
#  
#  SPARC: Isowaveguide
#  Faro, qua 13 mai 2026 11:17:33  
#  Written by Tordar
#  
#==================================================================

from os import system
from numpy import *
from matplotlib.pyplot import *
sys.path.append ("/home/orodrig/FORdoc/TRACEO/Python")
from wsparcenvfil import *
from fieldsparc   import *
from plotmovie    import *

print("SPARC, Isowaveguide calculations...")

# Case title
case_title = "Isovelocity_waveguide"

# Parameters
freq = 100  # Hz
zs =     0  # Source depth (m)

Dmax = 50 # Maximum depth (m)
Rmaxkm = 0.2   # Maximum range (km)
Rmax = Rmaxkm * 1000  # Maximum range (m)
# ==================================================================
# Source properties
# ==================================================================

tx_data = {
'n': 1,   # Number of sources
'zs': zs, # Source depths
'freq': freq
}

# ==================================================================
# Surface properties
# ==================================================================

# Last option:
# S => snapshot 
# R => range stack (horizontal array)
# D => depth stack (vertical array)
top_data = {
'options': "'NVW S'",
'nmesh': 500,   # Number of mesh points
'csigma': 0,  # RMS roughness
'zbottom': Dmax   # Depth at bottom of medium
}

# ==================================================================
# Sound speed profile properties
# ==================================================================
ssp_array = array([[-200,1500],[Dmax,1500]])
zi = ssp_array[:, 0]
ci = ssp_array[:, 1]

# Create arrays for SSP data
cs  = np.zeros_like(zi)   # Sound speed attenuation
rho = np.ones_like( zi)   # Density
apt = np.zeros_like(zi)  # Attenuation parameter (P-wave)
ast = np.zeros_like(zi)  # Attenuation parameter (S-wave)

# Assemble SSP data: 6 rows (zi, ci, cs, rho, apt, ast)
# The original Octave code does: ssp_data.data = [zi;ci;cs;rho;apt;ast]
ssp_data = {
'data': np.array([zi, ci, cs, rho, apt, ast])
}

# ==================================================================
# Bottom properties
# ==================================================================
# Yup, very clumsy, but it works... 
layer_data = array([[Dmax, 2500, 0, 1, 0, 0],
                    [200 , 2500, 0, 1, 0, 0]])
                    
bottom_data = {
'nlayers': 2,
'nmesh': 300,      # Number of mesh points
'csigma': 0,       # RMS roughness
'zbottom': 200.0,  # Depth at bottom of medium
'layer_data': layer_data
}

# ==================================================================
# Field parameters (receiver info)
# ==================================================================

nrd = 101
rx_data = {
'clow': 1000.0,   # Low cutoff frequency (CAN'T BE 0!)
'chigh': 1.0E7,  # High cutoff frequency
'rmaxkm': 1, # Maximum range (km)
'nrd': nrd,
'rd': np.linspace(0, 100, nrd)  # Receiver depths (m)
}

# ==================================================================
# Pulse info
# ==================================================================

# Pulse options strings:
# PULSE(1:1)
# P => Pseudo-Gaussian 
# R => Ricker wavelet
# A => Approximate Ricker Wavelet
# S => Single sine
# H => Hanning weighted four sine
# N => N-wave
# G => Gaussian 
# F => From a Source Time Series (STS) file 
# B => From a Source Time Series (STS) file backwards
# PULSE(2:2)
# H => perform a Hilbert transform of the source
# N => do not perform... 
# Hilbert transforming eliminates the left travelling wave
# PULSE(3:3)
# + => don't flip (recommended) 
# - => flip 
# PULSE(4:4)
# L => low cut filter 
# H => high cut filter
# B => both high and cut filter
# N => no cut
# Fmin => low  cut frequency 
# Fmax => high cut frequency 
# NRR => number of receiver ranges
# RR => receiver ranges

ntout = 60
pulse_data = {
'options': "'NH+N'",  # Pulse options string
'freqs': [0.0, 600.0],# Frequency range (should contain source freq)
'nrr': 1, # Number of receiver ranges
'rrkm': 0.002,   # Receiver ranges (km) - CAN'T BE LARGER THAN RMAX
'ntout': ntout,
'tout': np.linspace(0.002, 0.22, ntout),  # Output times (s)
'pars': [-0.1, 0.98, 0.0, 0.0, 0.0] # TSTART, TMULT, ALPHA, BETA, V
}

# ==================================================================
# Write the env file
# ==================================================================

wsparcenvfil('iso', case_title, tx_data, top_data, ssp_data, bottom_data, rx_data, pulse_data)

# ==================================================================
# Write the flp file (field parameters file)
# ==================================================================

flpoptions = "'RP'"
rminkm = 0.002
rmaxkm = 1.0
nr = 500

with open('iso.flp', 'w') as fid:
     fid.write(f"{flpoptions}\n")  # 'R/X (coord), Pos/Neg/Both'
     fid.write(f"{nr}\n")  # NR
     fid.write(f"{rminkm} {rmaxkm} /\n")  # RMIN, RMAX

system("sparc.exe iso")
fieldsparc('iso.grn')
plotmovie('iso.shd.mat')

print("done.")
