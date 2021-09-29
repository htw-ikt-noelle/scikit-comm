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
import copy

# upsampling factor
UP = 5

# contruct signal
sig_tx = comm.signal.Signal(n_dims=1)
sig_tx.symbol_rate = 50e6 #50e6

# generate bits
sig_tx.generate_bits(n_bits=2**12, seed=1)

# set constellation (modualtion format)
sig_tx.generate_constellation(format='QAM',order=16)
sig_tx.modulation_info = '16-QAM'

# create symbols
sig_tx.mapper()
# upsampling (without pulse shaping)
sig_tx.pulseshaper(upsampling=UP,pulseshape='rect')
# sig_tx.samples = sig_tx.symbols

### generate sig_rx
sig_rx = copy.deepcopy(sig_tx)

# add artificial delay 
sps = int(sig_rx.sample_rate[0] / sig_rx.symbol_rate[0])
delay = 5*sps
sig_rx.samples = sig_rx.samples[0][delay:]

# add amplitude noise
sig_rx.set_snr(30)
# add phase ambiguity (phase rotation and optional complex conjugation)
# sig_rx.samples = sig_rx.samples[0] * np.exp(1j*np.pi/3)

sig_rx.plot_constellation()

# add phase noise
# sig_rx.samples = comm.channel.add_phase_noise(sig_rx.samples[0],sig_rx.sample_rate[0],1e-3)['samples']

# downsample to 1sps
sig_rx.samples = sig_rx.samples[0][::UP]
sig_rx.sample_rate = (sig_rx.sample_rate[0] / UP)

# CPE
# cpe_results = comm.rx.carrier_phase_estimation_VV(sig_rx.samples[0], n_taps=31, filter_shape='wiener', mth_power=4, rho=.3)
# sig_rx.samples = cpe_results['rec_symbols']
# est_phase = cpe_results['phi_est']

# delay and phase ambiguity estimation and compensation
sig_rx = comm.rx.symbol_sequence_sync(sig_rx, dimension=-1)

# decide
sig_rx.decision()
# plot constellation after decision
sig_rx.plot_constellation()
# demap
sig_rx.demapper()



# since CPE and artificial delay throw away some samples, BER calculation can 
# only be performed on length of the shorter of the two bit vectors (sig_rx.samples)
length = len(sig_rx.samples[0])

ber_res = comm.rx.count_errors(sig_rx.bits[0][:length], sig_rx.samples[0][:length])
print('BER = {}'.format(ber_res['ber']))