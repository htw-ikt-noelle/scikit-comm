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


sps = 8 # samples per modulation symobol


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
# sig_tx.samples = sig_tx.symbols[0]
sig_tx.pulseshaper(upsampling=sps, pulseshape='rrc', roll_off=0.2)

################ channel ###########################
sig_rx = copy.deepcopy(sig_tx)


# tmp = copy.deepcopy(sig_rx)
# get samples from scope (repeat rx sequence)
ext = 8192*sps + 4000*sps
ratio_base = ext // sig_rx.samples[0].size
ratio_rem = ext % sig_rx.samples[0].size        
sig_rx.samples[0] = np.concatenate((np.tile(sig_rx.samples[0], ratio_base), sig_rx.samples[0][:ratio_rem]), axis=0)

# artificially delay received samples (cut away delay leading symbols)
delay = 10*sps
sig_rx.samples = sig_rx.samples[0][delay:]

# introduce ambiguity (phase shift / flip)
# sig_rx.samples = sig_rx.samples[0] * np.exp(-1j*np.pi/2)
# sig_rx.samples = np.conj(sig_rx.samples) * np.exp(-1j*1.58)

# TODO: implement once link performance is satisfactory without any noise
# AWGN

# phase noise


############# Rx #####################
# matched filter
sig_rx.samples = comm.filters.raised_cosine_filter(sig_rx.samples[0], sample_rate=sig_rx.sample_rate[0],
                                                   symbol_rate=sig_rx.symbol_rate[0], roll_off=0.2,
                                                   root_raised=True)

# artificially delay received samples (cut away delay leading symbols)
# cut_lead = 1000*sps
# cut_trail = 1000*sps
# sig_rx.samples = sig_rx.samples[0][cut_lead:-cut_trail]
crop = 10*sps
if crop != 0:
    sig_rx.samples = sig_rx.samples[0][crop:-crop]
else:
    sig_rx.samples = sig_rx.samples

comm.visualizer.plot_eye(sig_rx.samples[0][:500*sps],sample_rate=sig_rx.sample_rate[0], 
                          bit_rate=sig_rx.symbol_rate[0])

# correct for sampling instant
sig_rx.sampling_phase_adjustment()


comm.visualizer.plot_eye(sig_rx.samples[0][:500*sps],sample_rate=sig_rx.sample_rate[0], 
                         bit_rate=sig_rx.symbol_rate[0])

# sampling to 1 sps
sig_rx.samples = sig_rx.samples[0][::sps]


# CPE


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
    # symbols: only delay compensation is performed, symbols are then
    # independently decided and decided before counting errors against
    # rx.samples
    # sig_rx.symbols = np.roll(np.conj(sig_rx.symbols[0]), -int(symbol_delay_est)) * np.exp(-1j*phase_est)
    sig_rx.symbols = np.roll(np.conj(sig_rx.symbols[0]), -int(symbol_delay_est)) 
    # samples: ambiguity compensation
    # TODO: continue here
    sig_rx.samples = ???
else:
    # symbols: only delay compensation
    # sig_rx.symbols = np.roll(sig_rx.symbols[0], -int(symbol_delay_est)) * np.exp(1j*phase_est)
    sig_rx.symbols = np.roll(sig_rx.symbols[0], -int(symbol_delay_est))
    # samples: ambiguity compensation
    sig_rx.samples = ???
    
# generate reference bit sequence from manipulated symbol sequence by decision and demapping
sig_rx.symbols = comm.rx.decision(sig_rx.symbols[0], sig_rx.constellation[0])
sig_rx.bits = comm.rx.demapper(sig_rx.symbols[0], sig_rx.constellation[0])

########################################################################################

sig_rx.plot_constellation()

# decision and demapper
sig_rx.decision()
sig_rx.demapper()

# BER counting
ber_res = comm.rx.count_errors(sig_rx.bits[0], sig_rx.samples[0])
print(ber_res['ber'])
# print(np.sum(ber_res['err_idx']))

# plt.plot(ber_res['err_idx'])








