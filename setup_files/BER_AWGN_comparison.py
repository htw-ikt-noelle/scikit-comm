### comparison of theoretical BER curves with simulated BER over AWGN channel
### for different modulation formats

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

# contruct signal
sig_tx = comm.signal.Signal(n_dims=1)
sig_tx.symbol_rate = 32e9 

# range of SNR values in dB
snr_range = np.arange(20)
# QAM orders to generate
orders = [4,16,32]
# due to the implementation of generate bits, the number of generated bits must
# be a multiple of the least common multiple (LCM) of log2 of all the QAM orders
# that will be compared 
# ex. for LCM of 2,4,5,6,7,8 (QAM orders 4,16,32,64,128,256): LCM = 840
# TODO: implement a function that calculates the LCM for a list/array of QAM orders
LCM = 840

# generate bits
sig_tx.generate_bits(n_bits=10*LCM, seed=1)

### upsampling factor/sps
UP = 1
sig_tx.sample_rate = sig_tx.symbol_rate[0] * UP

# calc BER vectors per order
ber_order = []
ber_theoretical = []
for order in orders:
    # set constellation (modulation format)
    sig_tx.generate_constellation(format='QAM',order=order)
    # sig_tx.modulation_info = 'QAM'
    
    # create symbols
    sig_tx.mapper()
    # upsampling (without pulse shaping)
    sig_tx.pulseshaper(upsampling=UP,pulseshape='rect')
    
    # calc BER for different SNR values
    ber = []
    for i in snr_range:
    # add amplitude noise
        sig_rx = copy.deepcopy(sig_tx)  
        sig_rx.set_snr(int(i))
        # downsample to 1 sps if not already the case
        sig_rx.samples = sig_rx.samples[0][::UP]
        sig_rx.decision()
        sig_rx.demapper()
        ber.append(comm.rx.count_errors(sig_rx.bits[0], sig_rx.samples[0])['ber'])

    ber_order.append(np.asarray(ber))
    
    # calc theoretical BER
    ber_theoretical.append(comm.utils.ber_awgn(modulation_order=order))

# plot BER curves
fig, ax0 = plt.subplots()
ax0.plot(snr_range,ber_order[0],ber_theoretical[0])
ax0.set_yscale('log')
plt.legend(('simulated BER','theoretical BER'))
plt.xlabel('SNR [dB]')
plt.ylabel('BER')
if np.min(ber_order[0]) < 1e-5:
    plt.ylim((1e-5,0))
plt.grid()
plt.title('simulated vs. theoretical BER, 4-QAM, 1 sps, no MF')

fig, ax1 = plt.subplots()
ax1.plot(snr_range,ber_order[1],ber_theoretical[1])
ax1.set_yscale('log')
plt.legend(('simulated BER','theoretical BER'))
plt.xlabel('SNR [dB]')
plt.ylabel('BER]')
if np.min(ber_order[1]) < 1e-5:
    plt.ylim((1e-5,0))
plt.grid()
plt.title('simulated vs. theoretical BER, 16-QAM, 1 sps, no MF')

fig, ax2 = plt.subplots()
ax2.plot(snr_range,ber_order[2],ber_theoretical[2])
ax2.set_yscale('log')
plt.legend(('simulated BER','theoretical BER'))
plt.xlabel('SNR [dB]')
plt.ylabel('BER')
if np.min(ber_order[2]) < 1e-5:
    plt.ylim((1e-5,0))
plt.grid()
plt.title('simulated vs. theoretical BER, 32-QAM, 1 sps, no MF')
plt.show()

### comparison at 5 sps w/ matched filter

# upsampling factor/sps
UP = 5
sig_tx.sample_rate = sig_tx.symbol_rate[0] * UP

# calc BER vectors per order
ber_order = []
ber_theoretical = []
for order in orders:
    # set constellation (modulation format)
    sig_tx.generate_constellation(format='QAM',order=order)
    # sig_tx.modulation_info = 'QAM'
    
    # create symbols
    sig_tx.mapper()
    # upsampling (without pulse shaping)
    sig_tx.pulseshaper(upsampling=UP,pulseshape='rrc',roll_off=0.2)
    
    # calc BER for different SNR values
    ber = []
    for i in snr_range:
    # add amplitude noise
        sig_rx = copy.deepcopy(sig_tx)  
        sig_rx.set_snr(int(i))
        # matched filter at Rx
        sig_rx.raised_cosine_filter(roll_off=0.2,root_raised=True) 
        # downsample to 1 sps if not already the case
        sig_rx.samples = sig_rx.samples[0][::UP]
        sig_rx.decision()
        sig_rx.demapper()
        ber.append(comm.rx.count_errors(sig_rx.bits[0], sig_rx.samples[0])['ber'])

    ber_order.append(np.asarray(ber))
    
    # calc theoretical BER
    ber_theoretical.append(comm.utils.ber_awgn(modulation_order=order))
    
# plot BER curves
fig, ax3 = plt.subplots()
ax3.plot(snr_range,ber_order[0],ber_theoretical[0])
ax3.set_yscale('log')
plt.legend(('simulated BER','theoretical BER'))
plt.xlabel('SNR [dB]')
plt.ylabel('BER')
if np.min(ber_order[0]) < 1e-5:
    plt.ylim((1e-5,0))
plt.grid()
plt.title('simulated vs. theoretical BER, 4-QAM, 5 sps, MF')

fig, ax4 = plt.subplots()
ax4.plot(snr_range,ber_order[1],ber_theoretical[1])
ax4.set_yscale('log')
plt.legend(('simulated BER','theoretical BER'))
plt.xlabel('SNR [dB]')
plt.ylabel('BER]')
if np.min(ber_order[1]) < 1e-5:
    plt.ylim((1e-5,0))
plt.grid()
plt.title('simulated vs. theoretical BER, 16-QAM, 5 sps, MF')

fig, ax5 = plt.subplots()
ax5.plot(snr_range,ber_order[2],ber_theoretical[2])
ax5.set_yscale('log')
plt.legend(('simulated BER','theoretical BER'))
plt.xlabel('SNR [dB]')
plt.ylabel('BER')
if np.min(ber_order[2]) < 1e-5:
    plt.ylim((1e-5,0))
plt.grid()
plt.title('simulated vs. theoretical BER, 32-QAM, 5 sps, MF')
plt.show()