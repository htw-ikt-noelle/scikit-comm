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

value = np.arange(15)

ber = comm.utils.ber_awgn(value=value, type='SNR', modulation_format='QAM', modulation_order=16, symbol_rate=32e9, PDM=False)

print(ber)

plt.semilogy(value, ber)
plt.grid()