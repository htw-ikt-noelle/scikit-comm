# Used modules
import sys
import os
if not any(os.path.abspath('..') == p for p in sys.path): 
    print('adding comm module to path...')
    sys.path.insert(0, os.path.abspath('.'))
    print(os.path.abspath('.'))
import numpy as np  # Version 1.20.1
import matplotlib.pyplot as plt # Version 3.3.4
import comm as comm
import time

# # Setting AWG for test signal
# n = 1000
# phi = np.arange(n)/n * np.pi*2
# samples = np.array([np.cos(phi),np.sin(phi)])

# # write samples to AWG
# comm.instrument_control.write_samples_Tektronix_AWG70002B(samples, ip_address='192.168.1.21',
#                                                 sample_rate=[500e6],
#                                                 amp_pp=[0.5,0.3], channels=[1,2], 
#                                                 out_filter=['normal'],log_mode=True)

# Reading data form scope (only one channel)
sample_rate, wfm = comm.instrument_control.get_samples_Tektronix_MSO6B(word_length = 1)

print(sample_rate)
plt.plot(wfm[0])
plt.show()


# Reading data form scope (two channels)
sample_rate, wfm = comm.instrument_control.get_samples_Tektronix_MSO6B(channels = [1,2], word_length = 2)

print(sample_rate)
for idx,item in enumerate(wfm):
    plt.figure(idx)
    plt.plot(item)

plt.show()

# Reading data form scope (three channels)
sample_rate, wfm = comm.instrument_control.get_samples_Tektronix_MSO6B(channels = [1,2,3], word_length = 2)

print(sample_rate)
for idx,item in enumerate(wfm):
    plt.figure(idx)
    plt.plot(item)

plt.show()

# Reading data form scope (four channels)
sample_rate, wfm = comm.instrument_control.get_samples_Tektronix_MSO6B(channels = [1,2,3,4], word_length = 2)

print(sample_rate)
for idx,item in enumerate(wfm):
    plt.figure(idx)
    plt.plot(item)

plt.show()

# Reading data form scope (Onyl channel 1 and 3)
sample_rate, wfm = comm.instrument_control.get_samples_Tektronix_MSO6B(channels = [1,4], word_length = 2)

print(sample_rate)
for idx,item in enumerate(wfm):
    plt.figure(idx)
    plt.plot(item)

plt.show()