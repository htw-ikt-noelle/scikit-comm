# ESNR estimation - do not push!

import sys
import os
# if not any(os.path.abspath('..') == p for p in sys.path): 
#     print('adding comm module to path...')
#     sys.path.insert(0, os.path.abspath('..'))
import numpy as np
import matplotlib.pyplot as plt
import comm as comm

# global
TX_UPSAMPLE_FACTOR = 5
SNR = 9 #in dB

# construct signal
sig_tx = comm.signal.Signal(n_dims=1)
sig_tx.symbol_rate = 50e6 

# generate bits
sig_tx.generate_bits(n_bits=2**12, seed=1)

# set constellation (modulation format)
sig_tx.generate_constellation(order=4)
sig_tx.modulation_info = 'QPSK'

# create symbols
sig_tx.mapper()

# sig_tx.samples = sig_tx.symbols[0]
# sig_tx.sample_rate = sig_tx.symbol_rate[0]
# upsampling and pulseshaping
ROLL_OFF = 0.2
sig_tx.pulseshaper(upsampling=TX_UPSAMPLE_FACTOR, pulseshape='rrc', roll_off=[ROLL_OFF])

# add complex noise
sig_tx.set_snr(snr_dB=SNR,seed=None)

# RX matched filter
sig_tx.samples = comm.filters.raised_cosine_filter(sig_tx.samples[0],root_raised=True,roll_off=ROLL_OFF,
                                                   symbol_rate=sig_tx.symbol_rate[0],
                                                   sample_rate=sig_tx.sample_rate[0])


# downsample to 1 sps
sig_tx.samples = sig_tx.samples[0][::TX_UPSAMPLE_FACTOR]

sig_tx.plot_constellation()


# TODO: estimate SNR from (electrical) samples

# calc mean
symb_mean = np.mean(np.abs(sig_tx.samples[0]))

# calc variance
symb_var = np.mean((np.abs(sig_tx.samples[0])-symb_mean)**2)


# estimate SNR
snr_est = 10*np.log10((np.abs(symb_mean)**2)/(2*symb_var))

# estimated SNR is higher than actual SNR, the error increases as the actual
# SNR decreases; a modification of the estimation is performed:

if snr_est < 10:
    # introduce alias for sample vector for better legibility
    y = sig_tx.samples[0]
    
    # calc z (Qun & Jian, Eq. 8)
    z = snr_est
    #z = (np.mean(np.real(y)**2)+np.mean(np.imag(y)**2)) / (np.abs(np.mean(np.real(y)))**2+np.abs(np.mean(np.imag(y)))**2)
    
    # calc modified SNR, which is valid instead of snr_est
    #snr_est_mod = 1e4*(-0.041292958452235*(z**5)+2.66418532072905*(z**4)-6.86724072350538*(z**3)+8.84039993634297*(z**2)-5.68658561155135*z+1.464045795143920)
    snr_est_mod = np.sqrt((z-2.5)*39.2)-7

