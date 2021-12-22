import sys
import os
if not any(os.path.abspath('..') == p for p in sys.path): 
    print('adding comm module to path...')
    sys.path.insert(0, os.path.abspath('..'))

import comm as comm


############################################################
###################### Tx ##################################
############################################################

# signal parameters
LASER_LINEWIDTH = 0*600 # [Hz]
TX_UPSAMPLE_FACTOR = 5
EXPERIMENT = False
UPLOAD_SAMPLES = False
USE_PREDIST = False
SNR = 200

# contruct signal
sig_tx = comm.signal.Signal(n_dims=1)
sig_tx.symbol_rate = 50e6 

# generate bits
sig_tx.generate_bits(n_bits=2**14, seed=1)

# set constellation (modulation format)
# sig_tx.generate_constellation(format='PAM', order=2)
sig_tx.generate_constellation(format='PSK', order=2)

# create symbols
sig_tx.mapper()

sig_tx.samples = sig_tx.symbols[0]
sig_tx.sample_rate = sig_tx.symbol_rate[0]

sig_tx.plot_constellation(tit='test', axMax=2, hist=True, nBins=100)

sig_tx.set_snr(snr_dB=10)

sig_tx.plot_constellation(tit='test', axMax=2, hist=True, nBins=80)

