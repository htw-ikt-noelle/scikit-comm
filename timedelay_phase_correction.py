import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import comm as comm
import copy 
import warnings

## options for test
amount_symbols_wanted = 20 #must be higher than 4 for cropping
timedelay_in_sampels = 3
symbol_rate = 1
upsampling = 1
phase_offset = np.pi/3

#### GENERATE SIGNALS
#generate some signal
mod_format = 'QAM'
mod_order = 4
sig1 = comm.signal.Signal(n_dims=int(1))
sig1.symbol_rate[0] = symbol_rate
bits = (amount_symbols_wanted*2)*2
sig1.generate_bits(n_bits=bits,seed=1)
sig1.generate_constellation(format=mod_format,order=mod_order)
sig1.mapper()
sig1.pulseshaper(upsampling=upsampling,pulseshape='rrc', roll_off=0.2)
sig2 = copy.deepcopy(sig1)

#add phase offset to sig 2
sig2.samples[0] = sig2.samples[0]*np.exp(1j*phase_offset)

# roll and cut sampels (cut -> same legnth but some samples different for reality)
sig2.samples[0] = np.roll(sig2.samples[0],timedelay_in_sampels)
sig2.samples[0] = sig2.samples[0][int((amount_symbols_wanted/2)*upsampling):-int((amount_symbols_wanted/2)*upsampling)]
sig1.samples[0] = sig1.samples[0][int((amount_symbols_wanted/2)*upsampling):-int((amount_symbols_wanted/2)*upsampling)]

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
sig1.samples[0], sig2.samples[0], _ = comm.rx.comb_timedelay_compensation(sig1.samples[0], sig2.samples[0], sr=sig1.sample_rate[0], method="crop")
t = np.arange(0, (len(sig1.samples[0])/sig1.sample_rate[0]),1/sig1.sample_rate[0])

tit = "(cropped) Signals after time compensation"
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
sig1.samples[0], sig2.samples[0] = comm.rx.comb_phase_compensation(sig1.samples[0], sig2.samples[0])

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