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


# Test with one signal
n = 1000
phi = np.arange(n)/n * np.pi*2
samples = np.array([np.cos(phi)])


# write samples to AWG
comm.instrument_control.write_samples_Tektronix_AWG70002B(samples, ip_address='192.168.1.21',
                                                sample_rate=[1e6],
                                                amp_pp=[0.5], channels=[1], 
                                                out_filter=['normal'],log_mode=True)

'''
# Test with two signals
n = 1000
phi = np.arange(n)/n * np.pi*2
samples = np.array([np.cos(phi),np.sin(phi)])

# write samples to AWG
comm.instrument_control.write_samples_Tektronix_AWG70002B(samples, ip_address='192.168.1.21',
                                                sample_rate=[500e6],
                                                amp_pp=[0.5,0.3], channels=[1,2], 
                                                out_filter=['normal'],log_mode=True)
'''
'''
# Test with NaN in sample vector
n = 1000
phi = np.arange(n)/n * np.pi*2
samples = np.array([np.cos(phi)])
samples = np.array([np.append(samples,np.nan)])

# write samples to AWG
comm.instrument_control.write_samples_Tektronix_AWG70002B(samples, ip_address='192.168.1.21',
                                                sample_rate=[500e6],
                                                amp_pp=[0.5], channels=[1], 
                                                out_filter=['normal'],log_mode=True)

# Result: The NaN will be played by the AWG with maximum amplitude.
#         Also the AWG will give a warning: "Waveform contains invalid sample values" at his HOME screen.
#         The driver should prevent this by not accepting NaN values

# Test with Inf in sample vector
n = 1000
phi = np.arange(n)/n * np.pi*2
samples = np.array([np.sin(phi)])
samples = np.array([np.append(samples,np.Inf)])



# write samples to AWG
comm.instrument_control.write_samples_Tektronix_AWG70002B(samples, ip_address='192.168.1.21',
                                                sample_rate=[500e6],
                                                amp_pp=[0.5], channels=[1], 
                                                out_filter=['normal'],log_mode=True)

# Result: The inf will be played by the AWG with maximum amplitude.
#         The driver should prevent this by not accepting inf values
#         Positive infinity and negative infinity will be catched by np.isinf()


# Test with on signal, maximal length detection  
# Maximum Bytes for one command : 999,999,999 Bytes          
# Bytes per float = 4 Bytes
# Total length of signal for one command : 249999999 Samples                       
start = time.time()
print("START")
n = 150000000
#n = 10000000
phi = np.arange(n)/n * np.pi*2
samples = np.array([np.ones(n)])
end = time.time()
print("END")
print(end - start)

start = time.time()
print("START")
# write samples to AWG
comm.instrument_control.write_samples_Tektronix_AWG70002B(samples, ip_address='192.168.1.21',
                                                sample_rate=[500e6],
                                                amp_pp=[0.5], channels=[1], 
                                                out_filter=['normal'],log_mode=True)
end = time.time()
print("END")
print(end - start)
'''

