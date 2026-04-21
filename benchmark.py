"""
Test script for Kraken python Class

1. Cristini and Komatitsch, JASA, 2012 broadband figure
2. Munk profile with bathymetry bump TL figure

Ronan Serre
"""

import numpy as np
from kraken import Kraken


def ricker_wavelet(f, x, x0=0.0):

    y = (1.0 - 2 * np.pi ** 2 * (x - x0) ** 2 * f ** 2) * np.exp(- np.pi ** 2 * (x - x0) ** 2 * f ** 2)

    return y


# scenario
#job_title = ['1_Cristini-Komatitsch', '2_Mumpk'][1]
job_title = ['1_Cristini-Komatitsch', '2_Mumpk'][0]

# path towards acoustic toolbox
acoustic_toolbox_main_path = '/home/orodrig/FORdoc/at'

# first benchmark : constant ssp, constant bathymetry, broadband signal, Cristini and Komatitsch JASA 2012
if job_title == '1_Cristini-Komatitsch':
    # range array
    range_min = 0.0
    range_max = 400.0
    range_step = 10.0
    range_n = int((range_max - range_min) / range_step + 1)
    range_array = np.linspace(range_min, range_max, range_n)

    # depth array
    depth_min = 0.0
    depth_max = 400.0
    depth_step = 10.0
    depth_n = int((depth_max - depth_min) / depth_step + 1)
    depth_array = np.linspace(depth_min, depth_max, depth_n)

    # bathymetry values
    bathymetry_values = depth_max

    # constant ssp
    ssp_array = 1500.0

    # bottom geometry
    bottom = {"thickness": np.array([0.0, 400.0]),
              "p_speed": np.array([1500.0, 2400.0]),
              "absorption": np.array([0.5, 0.5]),
              "density": np.array([1.0 , 2.0])}

    # source(s) : define only by depth and frequency
    src_freq = 30.0
    src_depth = np.array([395.0])

    # receiver(s)
    rec_range_array = np.array([400.0])
    rec_range_n = rec_range_array.shape[0]

    rec_depth_array = np.array([380.0])
    rec_depth_n = rec_depth_array.shape[0]

    # broadband input signal
    fs = 200.0
    dt = 1.0 / fs
    tmax = 4.2
    nt = int(tmax * fs)
    t = np.linspace(0, nt - 1, nt) / fs
    py = ricker_wavelet(src_freq, t, 0.05)

    # amplification factor
    py *= 1.0e3

# second benchmark : range dependent TL on Munk profile with bump in bathymety (Mumpk)
elif job_title == '2_Mumpk':
    # range array
    range_min = 0.0
    range_max = 1e5
    range_step = 1000.0
    range_n = int((range_max - range_min) / range_step + 1)
    range_array = np.linspace(range_min, range_max, range_n)

    # depth array
    depth_min = 0.0
    depth_max = 5000.0
    depth_step = 100.0
    depth_n = int((depth_max - depth_min) / depth_step + 1)
    depth_array = np.linspace(depth_min, depth_max, depth_n)

    # bathymetry values
    bathymetry_values = depth_max - 1000.0 * np.exp(- (range_array - 0.5 * range_max) ** 2 / (5000.0 ** 2))

    # range dependent sound speed profil (Munk with small random variations)
    ssp_array = np.zeros((depth_n, range_n), dtype=float)
    ssp_eps = 0.00737
    c0_value_distribution = np.random.normal(1500.0, 10.0, size=range_n)
    c0_depth_distribution = np.random.normal(1300.0, 20.0, size=range_n)
    for r in range(range_n):
        x = 2.0 * (depth_array - c0_depth_distribution[r]) / c0_depth_distribution[r]
        ssp_array[:, r] = cenv = c0_value_distribution[r] * (1.0 + ssp_eps * (x - 1.0 + np.exp(-x)))

    # bottom geometry
    bottom = {"thickness": np.array([0.0, 2000.0]),
              "p_speed": np.array([1800.0, 2400.0]),
              "absorption": np.array([0.11, 0.18]),
              "density": np.array([2.0 , 1.6 ])}

    # source(s) : define only by depth and frequency
    src_freq = 100.0
    src_depth = np.array([1000.0])

    # receiver(s)
    rec_range_min = range_array[0]
    rec_range_max = range_array[-1]
    rec_range_n = 101
    rec_range_array = np.linspace(rec_range_min, rec_range_max, rec_range_n)
    rec_range_step = (rec_range_max - rec_range_min) / (rec_range_n - 1)

    rec_depth_step = 100.0
    rec_depth_min = rec_depth_step
    rec_depth_max = 4000.0
    rec_depth_n = int((rec_depth_max - rec_depth_min) / rec_depth_step + 1)
    rec_depth_array = np.linspace(rec_depth_min, rec_depth_max, rec_depth_n)

# initalise class
kraken = Kraken(job_title, at_path=acoustic_toolbox_main_path)

# assign environment
kraken.environment["range"] = range_array
kraken.environment["depth"] = depth_array
kraken.environment["bathymetry"] = bathymetry_values
kraken.environment["ssp"] = ssp_array
kraken.environment["bottom"] = bottom

# assign source(s)
kraken.src["depth"] = src_depth
kraken.src["frequency"] = src_freq

# assign receiver(s)
kraken.rec["range"] = rec_range_array
kraken.rec["depth"] = rec_depth_array

# assign broadband input signal
if job_title == '1_Cristini-Komatitsch':
    kraken.environment["signal"] = py
    kraken.environment["fs"] = fs

# generate environment
kraken.generate_env()

# run
kraken.run()

# generate figures
kraken.generate_figures()
