import time, copy

import numpy as np
import scipy.signal as ssignal
import matplotlib.pyplot as plt
import scipy.interpolate as sinterp
from mpl_toolkits.axisartist.axislines import AxesZero

import comm as comm

DAC_SR = 100e9

# contruct signal
sig_tx = comm.signal.Signal(n_dims=2)
sig_tx.symbol_rate = 10e9

TX_UPSAMPLE_FACTOR = DAC_SR / sig_tx.symbol_rate[0]

#%% # generate bits
sig_tx.generate_bits(n_bits=2**15)

#%% # set constellation (modulation format)
sig_tx.generate_constellation(format='QAM', order=4)

#%% # create symbols
sig_tx.mapper()

#%% # upsampling and pulseshaping
# sig_tx.pulseshaper(upsampling=TX_UPSAMPLE_FACTOR, pulseshape='rc', roll_off=0.2)
sig_tx.pulseshaper(upsampling=TX_UPSAMPLE_FACTOR, pulseshape='rect')



comm.visualizer.plot_poincare_sphere(sig_tx.samples[0],
                                     sig_tx.samples[1],
                                     decimation=1, 
                                     fNum = 1, 
                                     tit = 'Poincaré sphere simple',
                                     simple_plot=True)

comm.visualizer.plot_poincare_sphere(sig_tx.samples[0],
                                     sig_tx.samples[1],
                                     decimation=1, 
                                     fNum = 1, 
                                     tit = 'Poincaré sphere',
                                     simple_plot=False)