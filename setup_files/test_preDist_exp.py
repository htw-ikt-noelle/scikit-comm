# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 14:10:51 2021

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

sample_rate = 250e6
IP_AWG_AG33522A = '192.168.1.45'
AWG_AMPLpp = 4
f_max=70e6/(sample_rate/2)

IP_MSO_YK_DLM3024 = '192.168.1.14'

# generate noise samples
samples = comm.pre_distortion.generate_wn_probesignal(n_samples=2**17, f_max=f_max)
comm.visualizer.plot_spectrum(samples, sample_rate=sample_rate)

# upload to AWG
comm.instrument_control.write_samples_AWG33522A(np.real(samples), ip_address=IP_AWG_AG33522A, sample_rate=[sample_rate], offset=[0.0], amp_pp=[AWG_AMPLpp], channels=[1], out_filter=['normal'])
time.sleep(0.5)

comm.instrument_control.write_samples_AWG33522A(np.real(samples), ip_address=IP_AWG_AG33522A, sample_rate=[sample_rate], offset=[0.0], amp_pp=[AWG_AMPLpp], channels=[2], out_filter=['normal'])
time.sleep(0.5)

# get samples from scope
sr, samples_out = comm.instrument_control.get_samples_DLM2034(channels=[1, 2], address=IP_MSO_YK_DLM3024)
comm.visualizer.plot_spectrum(samples_out[0,:], sample_rate=sr)

# estimate transfer function
transfer_funtion = comm.pre_distortion.estimate_tf_welch(np.real(samples), sample_rate, samples_out[0,:], sr, f_max*sample_rate/2, nperseg=64, visualize=1)
tf = transfer_funtion['tf']
tf_inv = transfer_funtion['tf_inv']
tf_freq = transfer_funtion['freq']

plt.plot(20*np.log10(tf))
plt.show()
plt.plot(20*np.log10(tf_inv))
plt.axis([0,tf_inv.size, -20, 10])
plt.show()

# generate filter for saving and for the use with arbitrary filter
filtershape = np.column_stack((tf_freq, 20*np.log10(tf_inv)))

np.save('preDistFilter.npy', filtershape)
