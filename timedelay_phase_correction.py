from time import time
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import comm as comm
import copy 

## options
amount_symbols_wanted = 10
timedelay_in_sampels = 0
symbol_rate = 1
upsampling = 2

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
#sig1.pulseshaper(upsampling=upsampling, pulseshape="None")

#fix samples to something
#sig1.samples[0] = np.array([1+1j, 1-1j, -1+1j, -1-1j])

# crop/roll and clone sampels
sig2 = copy.deepcopy(sig1)
#add phase offset to sig 2
phase_offset = np.pi/5
sig2.samples[0] = sig2.samples[0]*np.exp(1j*phase_offset)

sig2.samples[0] = np.roll(sig2.samples[0],timedelay_in_sampels)
#sig2.samples[0] = sig2.samples[0][int((amount_symbols_wanted/2)*upsampling):-int((amount_symbols_wanted/2)*upsampling)]
#sig1.samples[0] = sig1.samples[0][int((amount_symbols_wanted/2)*upsampling):-int((amount_symbols_wanted/2)*upsampling)]

plt.figure()
plt.scatter(sig1.samples[0].real, sig1.samples[0].imag, label="Sig 1")
plt.scatter(sig2.samples[0].real, sig2.samples[0].imag, label="Sig 2")
plt.legend()
plt.grid()
plt.show()

t1 = np.arange(0, (len(sig1.samples[0])/sig1.sample_rate[0]),1/sig1.sample_rate[0])
t2 = np.arange(0, (len(sig2.samples[0])/sig1.sample_rate[0]),1/sig1.sample_rate[0])

plt.figure()
plt.title("Signals BEFORE time compensation")
plt.plot(t1,sig1.samples[0].real, label="sig1 real")
plt.plot(t1,sig1.samples[0].imag, label="sig1 imag",ls="--", marker="o")
plt.plot(t2,sig2.samples[0].real, label="sig2 real")
plt.plot(t2,sig2.samples[0].imag, label="sig2 imag",ls="--", marker="*")
plt.ylabel("Amplitude")
plt.xlabel("time [s]")
plt.legend()
plt.show()

correlation = signal.correlate(sig1.samples[0], sig2.samples[0], mode="full")
lags = signal.correlation_lags(sig1.samples[0].size, sig2.samples[0].size, mode="full")
lag = lags[np.argmax(correlation)]

time_offset = lag/sig1.sample_rate[0]

sig2.samples[0] = np.roll(sig2.samples[0], lag)

if lag >= 0:
    sig2.samples[0]= sig2.samples[0][lag:]
    sig1.samples[0] = sig1.samples[0][lag:]
    t1_crop = t1[lag:]
    t2_crop = t2[lag:]
elif lag < 0:
    sig2.samples[0] = sig2.samples[0][:lag]
    sig1.samples[0] = sig1.samples[0][:lag]
    t1_crop = t1[:lag]
    t2_crop = t2[:lag]

plt.figure()
plt.scatter(sig1.samples[0].real, sig1.samples[0].imag, label="Sig 1")
plt.scatter(sig2.samples[0].real, sig2.samples[0].imag, label="Sig 2")
plt.legend()
plt.grid()
plt.show()

plt.figure()
plt.title("Signals AFTER time compensation after cutting")
plt.plot(t1_crop,sig1.samples[0].real, label="sig1 real")
plt.plot(t1_crop,sig1.samples[0].imag, label="sig1 imag",ls="--", marker="o")
plt.plot(t2_crop,sig2.samples[0].real, label="sig2 real")
plt.plot(t2_crop,sig2.samples[0].imag, label="sig2 imag",ls="--", marker="*")
plt.ylabel("Amplitude")
plt.xlabel("time [s]")
plt.legend()
plt.show()

C = sum(sig1.samples[0]*np.conj(sig2.samples[0]))

print(np.angle(C))
print(phase_offset)

sig2.samples[0] = sig2.samples[0]*np.exp(1j*np.angle(C))

plt.figure()
plt.scatter(sig1.samples[0].real, sig1.samples[0].imag, label="Sig 1")
plt.scatter(sig2.samples[0].real, sig2.samples[0].imag, label="Sig 2", marker="*")
plt.legend()
plt.grid()
plt.show()

plt.figure()
plt.title("Signals AFTER time compensation/cutting/phase_compensation")
plt.plot(t1_crop,sig1.samples[0].real, label="sig1 real")
plt.plot(t1_crop,sig1.samples[0].imag, label="sig1 imag",ls="--", marker="o")
plt.plot(t2_crop,sig2.samples[0].real, label="sig2 real")
plt.plot(t2_crop,sig2.samples[0].imag, label="sig2 imag",ls="--", marker="*")
plt.ylabel("Amplitude")
plt.xlabel("time [s]")
plt.legend()
plt.show()

# #add noise
# sig2.samples[0] = comm.channel.set_snr(sig2.samples[0],snr_dB=5,sps=upsampling)
# sig1.samples[0] = comm.channel.set_snr(sig1.samples[0],snr_dB=5,sps=upsampling)

# t1 = np.arange(0, (len(sig1.samples[0])/sig1.sample_rate[0]),1/sig1.sample_rate[0])
# t2 = np.arange(0, (len(sig1.samples[0])/sig1.sample_rate[0]),1/sig2.sample_rate[0])

# plt.figure()
# plt.title("Signals BEFORE time compensation")
# plt.plot(t1,sig1.samples[0].real, label="sig1 real")
# plt.plot(t1,sig1.samples[0].imag, label="sig1 imag",ls="--")
# plt.plot(t2,sig2.samples[0].real, label="sig2 real")
# plt.plot(t2,sig2.samples[0].imag, label="sig2 imag",ls="--")
# plt.ylabel("Amplitude")
# plt.xlabel("time [s]")
# plt.legend()
# plt.show()

# correlation = signal.correlate(sig1.samples[0], sig2.samples[0], mode="full")
# lags = signal.correlation_lags(sig1.samples[0].size, sig2.samples[0].size, mode="full")
# lag = lags[np.argmax(correlation)]

# print("time offset: {}s".format(lag/sig1.sample_rate[0]))

# sig2.samples[0] = np.roll(sig2.samples[0], lag)

# plt.figure()
# plt.title("Signals AFTER time compensation")
# plt.plot(t1,sig1.samples[0].real, label="sig1 real")
# plt.plot(t1,sig1.samples[0].imag, label="sig1 imag",ls="--", marker="o")
# plt.plot(t2,sig2.samples[0].real, label="sig2 real")
# plt.plot(t2,sig2.samples[0].imag, label="sig2 imag",ls="--", marker="*")
# plt.ylabel("Amplitude")
# plt.xlabel("time [s]")
# plt.legend()
# plt.show()

# if lag >= 0:
#     sig2.samples[0] = sig2.samples[0][lag:]
#     t1 = t1[lag:]
#     sig1.samples[0] = sig1.samples[0][lag:]
#     t2 = t2[lag:]
# elif lag < 0:
#     sig2.samples[0] = sig2.samples[0][:lag]
#     t2 = t2[:lag]
#     sig1.samples[0] = sig1.samples[0][:lag]
#     t1 = t1[:lag]

# plt.figure()
# plt.title("Signals AFTER time compensation after cutting")
# plt.plot(t1,sig1.samples[0].real, label="sig1 real")
# plt.plot(t1,sig1.samples[0].imag, label="sig1 imag",ls="--", marker="o")
# plt.plot(t2,sig2.samples[0].real, label="sig2 real")
# plt.plot(t2,sig2.samples[0].imag, label="sig2 imag",ls="--", marker="*")
# plt.ylabel("Amplitude")
# plt.xlabel("time [s]")
# plt.legend()
# plt.show()
