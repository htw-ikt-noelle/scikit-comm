# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 17:04:51 2021

@author: UET-Labor
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


samples_in = comm.pre_distortion.generate_wn_probesignal(f_max=1, n_samples=2**17)

comm.visualizer.plot_spectrum(samples_in)

# FILTER = np.array([[0, 0.001,0.1, 0.5], [0, 0, -30, -10]],dtype='float').transpose()
# # FILTER = np.array([[0, 0.001,0.2], [0, 0, -30]],dtype='float').transpose()
# samples_out = comm.filters.filter_arbitrary(samples_in, FILTER, sample_rate = 1)
# samples_out = comm.pre_distortion.generate_wn_probesignal(sample_rate=250e6, n_samples=2**15, f_max=80e6)
samples_out = comm.filters.raised_cosine_filter(samples_in,sample_rate=1,symbol_rate=0.5)
samples_out = np.real(np.fft.ifft(np.fft.fft(samples_out, n=samples_in.size*1)))

comm.visualizer.plot_spectrum(samples_out)

results = comm.pre_distortion.estimate_tf_welch(samples_in, 1, samples_out, 1, f_max=0.5,visualize=226)

comm.visualizer.plot_signal(20*np.log10(results['tf']))
comm.visualizer.plot_signal(20*np.log10(results['tf_inv']))






