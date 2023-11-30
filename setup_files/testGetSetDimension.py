import time, copy, sys, pathlib

import numpy as np
import scipy.signal as ssignal
import matplotlib.pyplot as plt
import scipy.interpolate as sinterp

mod_path = str(pathlib.Path(__file__).parent.parent)
if not mod_path in sys.path:
    sys.path.append(mod_path)
import skcomm as skc

plt.ion()

DAC_SR = 16e9

# contruct signal
sig_tx = skc.signal.Signal(n_dims=3)
sig_tx.symbol_rate = 12.8e9

TX_UPSAMPLE_FACTOR = DAC_SR / sig_tx.symbol_rate[0]

# generate bits
sig_tx.generate_bits(n_bits=2**10*6, seed=[1,2,3])

# set constellation (modulation format)
sig_tx.generate_constellation(format='QAM', order=[4,16,64])

# create symbols
sig_tx.mapper()

# upsampling and pulseshaping
ROLL_OFF = 0.1
sig_tx.pulseshaper(upsampling=TX_UPSAMPLE_FACTOR, pulseshape='rrc', roll_off=ROLL_OFF)

# test get method
# sig_new = sig_tx.get_dimensions([0,2])
# sig_new = sig_tx.get_dimensions([0,1,2,0,1,2,0])
# sig_new = sig_tx.get_dimensions([0,0,0,0,0,0,0])

# # test add method
# sig_new = sig_tx.get_dimensions([0])
# sig_tx.add_dimension(sig_new,dim=3)

# sig_test = skc.signal.Signal()
# sig_tx.add_dimension(sig_test)

# test set method
sig_new = skc.signal.Signal(n_dims=2)
# sig_new = sig_tx.get_dimensions(dims=[2,0])
sig_tx.set_dimensions(sig_new, dims=[0,1])

print('1')