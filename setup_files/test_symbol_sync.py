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


# contruct signal
sig_tx = comm.signal.Signal(n_dims=1)
sig_tx.symbol_rate = 50e6 #50e6
sig_tx.sample_rate = 50e6 #50e6

# generate bits
sig_tx.generate_bits(n_bits=2**14)

# set constellation (modualtion format)
sig_tx.generate_constellation(order=4)
sig_tx.modulation_info = 'QPSK'

# create symbols
sig_tx.mapper()
# pulse shaping
sig_tx.samples = sig_tx.symbols[0]

############# Rx #####################
sig_rx = sig_tx

# get samples from scope (repeat rx sequence)
ext = 20000
ratio_base = ext // sig_rx.samples[0].size
ratio_rem = ext % sig_rx.samples[0].size        
sig_rx.samples[0] = np.concatenate((np.tile(sig_rx.samples[0], ratio_base), sig_rx.samples[0][:ratio_rem]), axis=0)

# artificially delay received samples (cut away delay leading symbols)
delay = 1059
sig_rx.samples = sig_rx.samples[0][delay:]

# introduce ambiguity (phase shift / flip)
sig_rx.samples = sig_rx.samples[0] * np.exp(-1j*1)
# sig_rx.samples = np.conj(sig_rx.samples) * np.exp(-1j*1.58)

########################## compensation algorithm ########################################
# TODO: put into rx function

# determine symbol delay and constellation ambiguity (phase rotation and conj)
corr_len = sig_tx.symbols[0].size
# complex correlation of rx_symbols and tx_smybols (only one symbol block needed)
corr_norm = ssignal.correlate(sig_rx.samples[0][:corr_len], sig_tx.symbols[0], mode='same')
# complex correlation of conjugated rx_symbols and tx_symbols (only one symbol block needed)
corr_conj = ssignal.correlate(np.conj(sig_rx.samples[0][:corr_len]), sig_tx.symbols[0], mode='same')

# decide which correlation is larger and determine delay index and phase
if np.max(np.abs(corr_norm)) > np.max(np.abs(corr_conj)):
    symbols_conj = False
    idx = np.argmax(np.abs(corr_norm))    
    phase_est = np.angle(corr_norm[idx])
else:
    symbols_conj = True
    idx = np.argmax(np.abs(corr_conj))
    phase_est = np.angle(corr_conj[idx])

# determine symbol delay from index depending on the location of the 
# correlation peak (either left or right from the center)
if idx <= corr_len/2:
    symbol_delay_est = int(corr_len/2 - idx)
else:
    symbol_delay_est = int(corr_len - idx + corr_len/2)

print('conjugated:{}, delay={}, phase={}'.format(symbols_conj, symbol_delay_est, phase_est))

# plt.plot(np.abs(corr_norm))
# plt.plot(np.abs(corr_conj))
# plt.show()

# manipulate logical reference symbol sequence in order to compensate for 
# delay and ambiguity
if symbols_conj:
    sig_rx.symbols = np.roll(np.conj(sig_rx.symbols[0]), -int(symbol_delay_est)) * np.exp(-1j*phase_est)
else:
    sig_rx.symbols = np.roll(sig_rx.symbols[0], -int(symbol_delay_est)) * np.exp(1j*phase_est)
# generate reference bit sequence from manipulated symbol sequence by decision and demapping
sig_rx.symbols = comm.rx.decision(sig_rx.symbols[0], sig_rx.constellation[0])
sig_rx.bits = comm.rx.demapper(sig_rx.symbols[0], sig_rx.constellation[0])

########################################################################################

# decision and demapper
sig_rx.decision()
sig_rx.demapper()

# BER counting
ber_res = comm.rx.count_errors(sig_rx.bits[0], sig_rx.samples[0])
print(ber_res['ber'])
# print(np.sum(ber_res['err_idx']))

# plt.plot(ber_res['err_idx'])








