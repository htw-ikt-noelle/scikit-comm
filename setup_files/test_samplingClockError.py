import time
import numpy as np
import scipy.signal as ssignal
import scipy.interpolate as sinterp
import matplotlib.pyplot as plt
import sys
import comm as comm


###################### Tx ##################################
# signal parameters
TX_UPSAMPLE_FACTOR = 5

# contruct signal
sig_tx = comm.signal.Signal(n_dims=1)
sig_tx.symbol_rate = 50e6 #50e6

# generate bits
sig_tx.generate_bits(n_bits=2**11)

# set constellation (modualtion format)
sig_tx.generate_constellation(order=4)
sig_tx.modulation_info = 'QPSK'

# create symbols
sig_tx.mapper()

# upsampling (to 5*50e6 --> 250 MSa/s)  and pulseshaping
ROLL_OFF = 0.1
sig_tx.pulseshaper(upsampling=TX_UPSAMPLE_FACTOR, pulseshape='rrc', roll_off=[ROLL_OFF])



#############################################################################
######################## Rx #################################################

# contruct rx signal
sig_rx = comm.signal.Signal(n_dims=1)
sig_rx.symbol_rate = sig_tx.symbol_rate
sig_rx.sample_rate = sig_tx.sample_rate[0]
sig_rx.samples[0] = sig_tx.samples[0]

# add artificial sample clock error
ratio = 1.001 # ratio of sampling frequency missmatch 
samples = sig_rx.samples[0]
n_old = np.size(samples, axis=0)
t_old = np.arange(n_old) / sig_rx.sample_rate[0]
n_new = int(np.round(ratio * n_old))
t_new = np.linspace(start=t_old[0], stop=t_old[-1], num=n_new, endpoint=True)
sr_new = 1 / (t_new[1] - t_new[0])
# interpolate signal at different timing / sampling instants
f = sinterp.interp1d(t_old, samples)
samples_new = f(t_new)
sig_rx.samples[0] = samples_new

# "standard" coherent complex baseband signal processing
# Rx matched filter
sig_rx.raised_cosine_filter(roll_off=ROLL_OFF,root_raised=True) 

# sampling adjustment
BLOCK_SIZE = 500 # number of blocks... -1 for only one block
samples = sig_rx.samples[0]

if BLOCK_SIZE == -1:
    # sampling phase adjustment (once per simulation) --> timing recovery
    results = comm.rx.sampling_phase_adjustment(samples, sample_rate=sig_rx.sample_rate[0], symbol_rate=sig_rx.symbol_rate[0])
    samples = results['samples_out']
    shift = results['est_shift']
else:
    # sampling clock adjustment (multiple times per simulation, multiple blocks)
    # --> crude sampling frequency offset correction
    # cut down to multiple of block size
    n = np.size(samples, axis=0)
    n_new = int(np.floor(n / BLOCK_SIZE) * BLOCK_SIZE)
    samples_array = np.reshape(samples[:n_new], (-1, BLOCK_SIZE))

    results = list()
    # run sampling_phase_adjustment multiple times
    for idx, samples in enumerate(samples_array):
        results.append(comm.rx.sampling_phase_adjustment(samples, sample_rate=sig_rx.sample_rate[0], symbol_rate=sig_rx.symbol_rate[0]))
    # generate output samples and shifts as ndarrays
    samples = np.asarray([result['samples_out'] for result in results]).reshape(-1)
    shifts = [result['est_shift'] for result in results]
    # plt.stem(np.asarray(shifts[:]),use_line_collection=True)

sig_rx.samples[0] = samples
                  

# sig_rx.samples[0] = results['samples_out']
# print(results['est_shift'])
# sig_rx.plot_eye()
# sig_rx.sampling_phase_adjustment()
# sig_rx.plot_eye()

# sampling
START_SAMPLE = 0
sps = sig_rx.sample_rate[0] / sig_rx.symbol_rate[0] # CHECK FOR INTEGER SPS!!!
rx_symbols = sig_rx.samples[0][START_SAMPLE::int(sps)]

comm.visualizer.plot_constellation(rx_symbols)