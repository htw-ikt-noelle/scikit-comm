# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 14:29:12 2021

@author: noelle
"""
import sys
import os
if not any(os.path.abspath('..') == p for p in sys.path): 
    print('adding comm module to path...')
    sys.path.insert(0, os.path.abspath('..'))
import time
import numpy as np
import scipy.signal as ssignal
import matplotlib.pyplot as plt
import scipy.interpolate as sinterp
import comm as comm
#from scipy.signal.signaltools import wiener as wiener



# contruct signal
sig_tx = comm.signal.Signal(n_dims=1)
sig_tx.symbol_rate = 50e6 #50e6

# generate bits
sig_tx.generate_bits(n_bits=2**10)

# generate symols
sig_tx.constellation[0] = np.asarray([0,1])
sig_tx.mapper()


upsam = 64
#samples_up = ssignal.upfirdn(np.asarray([1]), sig_tx.symbols[0], up=upsam, down=1)
sig_tx.pulseshaper(upsampling=upsam, pulseshape='None')
sig_tx.sample_rate[0] = upsam * sig_tx.symbol_rate[0]


# sig = comm.filters.raised_cosine_filter(sig_tx.samples[0], sample_rate=sig_tx.sample_rate[0], 
#                                         symbol_rate=sig_tx.symbol_rate[0], roll_off=0.0, 
#                                         length=-1, root_raised=True, domain='freq')

# sig = comm.filters.moving_average(sig_tx.samples[0], average=63, domain='freq')

# sig = comm.filters.windowed_sinc(sig_tx.samples[0], fc=0.1, order=-1, window=None)

sig = comm.filters.ideal_lp(sig_tx.samples[0], 0.1)['samples_out']


sig = sig[0:1024*upsam]

comm.visualizer.plot_eye(sig, sample_rate=sig_tx.sample_rate[0], 
                         bit_rate=sig_tx.symbol_rate[0], offset = 0, 
                         fNum = 1, tit = 'eye diagramm')

comm.visualizer.plot_spectrum(sig, sample_rate=sig_tx.sample_rate[0],
                              fNum=1, scale='logNorm', tit='spectrum')