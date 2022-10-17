import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import comm as comm
import copy 
import warnings
import scipy.signal as ssignal
import scipy.interpolate as sinterp

def add_subsample_time_delay(samples, sample_rate, subsample_delay=0.2/1):
    t_orig = np.arange(0, (len(samples)/sample_rate),1/sample_rate)
    t_off = t_orig+subsample_delay
    f = sinterp.interp1d(t_off, samples, kind='cubic', fill_value="extrapolate")
    r_samples = f(t_orig) 

    return r_samples

## options for test
amount_symbols_wanted = 6 #must be higher than 4 for cropping
timedelay_in_sampels = 0
symbol_rate = 2
upsampling = 2
phase_offset = 0
#phase_offset = np.pi/3

#### GENERATE SIGNALS
#generate some signal
mod_format = 'QAM'
mod_order = 4
sig1 = comm.signal.Signal(n_dims=int(1))
sig1.symbol_rate[0] = symbol_rate
bits = (amount_symbols_wanted*2)
sig1.generate_bits(n_bits=bits,seed=1)
sig1.generate_constellation(format=mod_format,order=mod_order)
sig1.mapper()
sig1.pulseshaper(upsampling=upsampling,pulseshape='rrc', roll_off=0.2)
sig2 = copy.deepcopy(sig1)

sig2.samples[0] = add_subsample_time_delay(sig2.samples[0], sig2.sample_rate[0], 0.2/sig2.symbol_rate[0])

results = comm.rx.sampling_phase_adjustment(sig2.samples[0], sample_rate=sig2.sample_rate[0], symbol_rate=sig2.symbol_rate[0], shift_dir='both')
samples_out = results['samples_out']
print(results['est_shift'])

t_orig = np.arange(0, (len(sig1.samples[0])/sig1.sample_rate[0]),1/sig1.sample_rate[0])

plt.figure()
plt.plot(t_orig, sig1.samples[0].real, label="sig 1 ohne subsampling offset real")
plt.plot(t_orig, sig1.samples[0].imag, label="sig 1 ohne subsampling offset imag")
plt.plot(t_orig, sig2.samples[0].real, label="sig 2 mit subsampling offset real")
plt.plot(t_orig, sig2.samples[0].imag, label="sig 2 mit subsampling offset imag")
plt.legend()
plt.show()

#samples_out = comm.filters.time_shift(sig2.samples[0], 20, -0.1)

plt.figure()
plt.plot(t_orig, sig1.samples[0].real, label="sig 1 ohne subsampling offset real")
plt.plot(t_orig, sig1.samples[0].imag, label="sig 1 ohne subsampling offset imag")
plt.plot(t_orig, samples_out.real, label="sig 2 mit subsampling offset real")
plt.plot(t_orig, samples_out.imag, label="sig 2 mit subsampling offset imag")
plt.legend()
plt.show()

print("end")