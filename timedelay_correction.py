from time import time
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import comm as comm
import copy 


def comp_timedelay(samples1, samples2, samplerate):
    """
    Compensate the time delay of different time shifted signals.
    """
    t1 = np.arange(0, (len(samples1)/samplerate),1/samplerate)
    t2 = np.arange(0, (len(samples1)/samplerate),1/samplerate)

    correlation = signal.correlate(samples1, samples2, mode="full")
    lags = signal.correlation_lags(samples1.size, samples2.size, mode="full")
    lag = lags[np.argmax(correlation)]

    time_offset = lag/samplerate

    samples2_comp = np.roll(samples2, lag)

    if lag >= 0:
        samples2_crop = samples2_comp[lag:]
        samples1_crop = samples1[lag:]
        t1_crop = t1[lag:]
        t2_crop = t2[lag:]
    elif lag < 0:
        samples2_crop = samples2_comp[:lag]
        samples1_crop = samples1[:lag]
        t1_crop = t1[:lag]
        t2_crop = t2[:lag]
 
    
    result_array = {
        "samples1": samples1,
        "samples2": samples2,
        "time 1": t1,
        "time 2": t2,
        "samples2_comp": samples2_comp,
        "samples2_crop": samples2_crop,
        "samples1_crop": samples1_crop,
        "t1_crop": t1_crop,
        "t2_crop": t2_crop,
        "time_offset": time_offset
    }

    return result_array

## options
amount_symbols_wanted = 1024
timedelay_in_sampels = 10
symbol_rate = 23
upsampling = 2

#generate some signal
mod_format = 'QAM'
mod_order = 4
sig1 = comm.signal.Signal(n_dims=int(1))
sig1.symbol_rate[0] = symbol_rate
bits = (amount_symbols_wanted*2)*2
sig1.generate_bits(n_bits=bits,seed=None)
sig1.generate_constellation(format=mod_format,order=mod_order)
sig1.mapper()
sig1.pulseshaper(upsampling=upsampling,pulseshape='rrc', roll_off=0.2)

# crop/roll and clone sampels
sig2 = copy.deepcopy(sig1)
sig2.samples[0] = np.roll(sig2.samples[0],timedelay_in_sampels)
sig2.samples[0] = sig2.samples[0][int((amount_symbols_wanted/2)*upsampling):-int((amount_symbols_wanted/2)*upsampling)]
sig1.samples[0] = sig1.samples[0][int((amount_symbols_wanted/2)*upsampling):-int((amount_symbols_wanted/2)*upsampling)]

timed_dict = comp_timedelay(sig1.samples[0], sig2.samples[0], sig1.sample_rate[0])

print(timed_dict["time_offset"])

plt.figure()
plt.title("Signals AFTER time compensation after cutting")
plt.plot(timed_dict["t1_crop"],timed_dict["samples1_crop"].real, label="sig1 real")
plt.plot(timed_dict["t1_crop"],timed_dict["samples1_crop"].imag, label="sig1 imag",ls="--", marker="o")
plt.plot(timed_dict["t2_crop"],timed_dict["samples2_crop"].real, label="sig2 real")
plt.plot(timed_dict["t2_crop"],timed_dict["samples2_crop"].imag, label="sig2 imag",ls="--", marker="*")
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
