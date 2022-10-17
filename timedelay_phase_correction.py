import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import comm as comm
import copy 
import warnings
import scipy.signal as ssignal
import scipy.interpolate as sinterp

## options for test
amount_symbols_wanted = 20 #must be higher than 4 for cropping
timedelay_in_sampels = 6
subsample_relative_delay = 0.2
symbol_rate = 1
upsampling = 2
#phase_offset = 0
phase_offset = np.pi/3

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

#add phase offset to sig 2
sig2.samples[0] = sig2.samples[0]*np.exp(1j*phase_offset)

# roll 
sig2.samples[0] = np.roll(sig2.samples[0],timedelay_in_sampels)

# add subsample time delay
def add_subsample_time_delay(samples, sample_rate, subsample_delay=0.2/1):
    t_orig = np.arange(0, (len(samples)/sample_rate),1/sample_rate)
    t_off = t_orig+subsample_delay
    f = sinterp.interp1d(t_off, samples, kind='cubic', fill_value="extrapolate")
    r_samples = f(t_orig) 

    return r_samples

sig2.samples[0] = add_subsample_time_delay(sig2.samples[0], sig2.sample_rate[0], subsample_delay=subsample_relative_delay/sig2.sample_rate[0])

#sig2.samples[0] = sig2.samples[0][int((amount_symbols_wanted/2)*upsampling):-int((amount_symbols_wanted/2)*upsampling)]
#sig1.samples[0] = sig1.samples[0][int((amount_symbols_wanted/2)*upsampling):-int((amount_symbols_wanted/2)*upsampling)]

t1 = np.arange(0, (len(sig1.samples[0])/sig1.sample_rate[0]),1/sig1.sample_rate[0])
t2 = np.arange(0, (len(sig2.samples[0])/sig1.sample_rate[0]),1/sig1.sample_rate[0])

tit = "Signals before comepensating time + phase"
plt.figure()
plt.title(tit+str(", constellation"))
plt.scatter(sig1.samples[0].real, sig1.samples[0].imag, label="Sig 1")
plt.scatter(sig2.samples[0].real, sig2.samples[0].imag, label="Sig 2")
plt.legend()
plt.grid()
plt.show()

plt.figure()
plt.title(tit+str(", timebase"))
plt.plot(t1,sig1.samples[0].real, label="sig1 real")
plt.plot(t1,sig1.samples[0].imag, label="sig1 imag",ls="--", marker="o")
plt.plot(t2,sig2.samples[0].real, label="sig2 real")
plt.plot(t2,sig2.samples[0].imag, label="sig2 imag",ls="--", marker="*")
plt.ylabel("Amplitude")
plt.xlabel("time [s]")
plt.legend()
plt.show()

#### COMPENSATE FOR TIME DELAY
sig1.samples[0], sig2.samples[0], lag = comm.rx.comb_timedelay_compensation(sig1.samples[0], sig2.samples[0], method="crop", xcorr="abs")

t = np.arange(0, (len(sig1.samples[0])/sig1.sample_rate[0]),1/sig1.sample_rate[0])


tit = "(zeropad) Signals after time compensation"
plt.figure()
plt.title(tit+str(", constellation"))
plt.scatter(sig1.samples[0].real, sig1.samples[0].imag, label="Sig 1")
plt.scatter(sig2.samples[0].real, sig2.samples[0].imag, label="Sig 2")
plt.legend()
plt.grid()
plt.show()

plt.figure()
plt.title(tit+str(", timebase"))
plt.plot(t,sig1.samples[0].real, label="sig1 real")
plt.plot(t,sig1.samples[0].imag, label="sig1 imag",ls="--", marker="o")
plt.plot(t,sig2.samples[0].real, label="sig2 real")
plt.plot(t,sig2.samples[0].imag, label="sig2 imag",ls="--", marker="*")
plt.ylabel("Amplitude")
plt.xlabel("time [s]")
plt.legend()
plt.show()

#### COMPENSATE FOR SUBSAMPLE DELAY
results = comm.rx.sampling_phase_adjustment(sig2.samples[0], sample_rate=sig2.sample_rate[0], symbol_rate=sig2.symbol_rate[0], shift_dir='both')
sig2.samples[0] = results['samples_out']
print(results['est_shift'])

tit = "(zeropad) Signals after subsample time compensation"
plt.figure()
plt.title(tit+str(", constellation"))
plt.scatter(sig1.samples[0].real, sig1.samples[0].imag, label="Sig 1")
plt.scatter(sig2.samples[0].real, sig2.samples[0].imag, label="Sig 2")
plt.legend()
plt.grid()
plt.show()

plt.figure()
plt.title(tit+str(", timebase"))
plt.plot(t,sig1.samples[0].real, label="sig1 real")
plt.plot(t,sig1.samples[0].imag, label="sig1 imag",ls="--", marker="o")
plt.plot(t,sig2.samples[0].real, label="sig2 real")
plt.plot(t,sig2.samples[0].imag, label="sig2 imag",ls="--", marker="*")
plt.ylabel("Amplitude")
plt.xlabel("time [s]")
plt.legend()
plt.show()

# ### COMPENSATE FOR PHASE
#print("should be the phase: "+str(np.rad2deg(phase_offset)))
sig1.samples[0], sig2.samples[0], phase_est = comm.rx.comb_phase_compensation(sig1.samples[0], sig2.samples[0])

tit = "Signals after time + phase compensation"
plt.figure()
plt.title(tit+str(", constellation"))
plt.scatter(sig1.samples[0].real, sig1.samples[0].imag, label="Sig 1")
plt.scatter(sig2.samples[0].real, sig2.samples[0].imag, label="Sig 2", marker="*")
plt.legend()
plt.grid()
plt.show()

plt.figure()
plt.title(tit+str(", timebase"))
plt.plot(t,sig1.samples[0].real, label="sig1 real")
plt.plot(t,sig1.samples[0].imag, label="sig1 imag",ls="--", marker="o")
plt.plot(t,sig2.samples[0].real, label="sig2 real")
plt.plot(t,sig2.samples[0].imag, label="sig2 imag",ls="--", marker="*")
plt.ylabel("Amplitude")
plt.xlabel("time [s]")
plt.legend()
plt.show()

print("[sample] sample delay est: {:.1f}, sample delay soll: {:.1f}".format(lag,timedelay_in_sampels))
print("[seconds????] sub sample delay est: {:.4f}, sub sample delay soll: {:.4f}".format(results['est_shift'],subsample_relative_delay/sig1.symbol_rate[0]))
print("[grad] phase est: {:.2f}, phase soll: {:.2f}".format(np.rad2deg(phase_est),np.rad2deg(phase_offset)))