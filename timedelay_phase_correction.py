import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import comm as comm
import copy 
import warnings

def comb_timedelay_compensation(x, y, sr=1, method="zeropad"):
    """
    Compensate the delay on sample basis between x and y. This is done by using
    crosscorrelation to estimate the "lag" between both numpy arrays (which can be complex).
    Specifing samplerate sr correctly results in calculated delay time between both.  
    There are two modes of timedelay compensation results:

        - "CROP" means, only overlapping values will be returned, values which are not present in 
        both (the samples in the lag time regime) are cropped.

        - "ZEROPAD" works with zeropadding and add zeros to the lag samples of the different vector.
    
    If the returned time value is positive, y is delayed to x and vice versa if lag is negative.
    
    Please note: This function is designed to work with arrays of complex values. (Unscaled) real valued arrays 
    leading to a less precise correlation and therefore possible errors.
    
    Parameters
    ----------
    x : np.array
        First vector
    y : np.array
        Second vector
    sr : int, optional
        sample rate of the vectors. Default is 1.
    method : string, optional
        method for return arrays. See docstring for help. Default is "zeropad". 
    

    Returns
    -------
    x : np.array
        First vector, timeshifted by method
    y : np.array
        Second vector, timeshifted by method
    time: float
        delay time in [s] between both vectors
    """

    if x.dtype != np.complex128 or y.dtype != np.complex128:
        warnings.warn("Combining timedelay compensation: arrays are not complex -> less precise xcorr if unscaled arrays.")

    correlation = signal.correlate(y, x, mode="full")
    lags = signal.correlation_lags(x.size, y.size, mode="full")
    lag = lags[np.argmax(correlation)]
    # if lag is positive, y is delayed to x and vice versa if lag is negative

    # use crop - cut away not matching data
    if method == "crop":
        #apply estimated time in samples
        y = np.roll(y, -lag)

        #cut signal, only return matching area wihtout lag samples
        if lag < 0:
            y = y[abs(lag):]
            x = x[abs(lag):]
        elif lag > 0:
            y = y[:-lag]
            x = x[:-lag]

    # use zeropadding - dont throw away any data, shift to compensate lag and fill zeros
    elif method == "zeropad":

        #build up some placeholders
        x_temp = np.zeros((len(x)+abs(lag)),dtype=x.dtype)
        y_temp = np.zeros((len(y)+abs(lag)),dtype=y.dtype)

        #fill placeholders as matching
        if lag > 0:
            y_temp[0:-lag] = y
            x_temp[lag:] = x
        elif lag < 0:
            y_temp[abs(lag):] = y
            x_temp[0:-abs(lag)] = x
        elif lag == 0:
            x_temp = x
            y_temp = y

        x = x_temp
        y = y_temp

    else:
        raise Exception("timedelay compensation method must be crop or zeropad!")

    return x, y, lag/sr

def comb_phase_compensation(x, y):
    """
    Compensate the phase on np.array y in respect to np.array x.
    Build up a complex correlation factor where the argument represent the phase offset between x and y.

    Sources:
    [1] Rao et. al., Toward Practical Digital Phase Alignment for Coherent Beam Combining in Multi-Aperture Free Space Coherent Optical Receivers, 2020
    [2] Geisler et. al., Multi-aperture digital coherent combining for free-space optical communication receivers, 2016

    Parameters
    ----------
    x : np.array
        First vector
    y : np.array
        Second vector
    
    Returns
    -------
    x : np.array
        First vector
    y : np.array
        Second vector, with compensated phase in respect to x
    """

    C = sum(x*np.conj(y))
    y = y*np.exp(1j*np.angle(C))
    
    return x, y

## options for test
amount_symbols_wanted = 6 #must be higher than 4 for cropping
timedelay_in_sampels = -1
symbol_rate = 1
upsampling = 20

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
#sig1.pulseshaper(upsampling=upsampling, pulseshape="None")

#fix samples to something
#sig1.samples[0] = np.array([1+1j, 1-1j, -1+1j, -1-1j])

# crop/roll and clone sampels
sig2 = copy.deepcopy(sig1)

#add phase offset to sig 2
phase_offset = np.pi/7
sig2.samples[0] = sig2.samples[0]*np.exp(1j*phase_offset)

#if we have different phase offset in signal
#sig2.samples[0][0:5] = sig2.samples[0][0:5]*np.exp(1j*phase_offset)
#sig2.samples[0][5:] = sig2.samples[0][5:]*np.exp(1j*(phase_offset+np.pi/3))

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

#### COMPENSATE FOR TIME DELAY (CROPPING METHOD)
#do a cross correlation, get back assumed lag in samples between both signals
# correlation = signal.correlate(sig1.samples[0], sig2.samples[0], mode="full")
# lags = signal.correlation_lags(sig1.samples[0].size, sig2.samples[0].size, mode="full")
# lag = lags[np.argmax(correlation)]

# time_offset = lag/sig1.sample_rate[0]

# #apply estimated time in samples
# sig2.samples[0] = np.roll(sig2.samples[0], lag)

# #cut signal, only show matching area
# if lag >= 0:
#     sig2.samples[0]= sig2.samples[0][lag:]
#     sig1.samples[0] = sig1.samples[0][lag:]
#     t1_crop = t1[lag:]
#     t2_crop = t2[lag:]
# elif lag < 0:
#     sig2.samples[0] = sig2.samples[0][:lag]
#     sig1.samples[0] = sig1.samples[0][:lag]
#     t1_crop = t1[:lag]
#     t2_crop = t2[:lag]

# tit = "(cropped) Signals after time compensation"
# plt.figure()
# plt.title(tit+str(", constellation"))
# plt.scatter(sig1.samples[0].real, sig1.samples[0].imag, label="Sig 1")
# plt.scatter(sig2.samples[0].real, sig2.samples[0].imag, label="Sig 2")
# plt.legend()
# plt.grid()
# plt.show()

# plt.figure()
# plt.title(tit+str(", timebase"))
# plt.plot(t1_crop,sig1.samples[0].real, label="sig1 real")
# plt.plot(t1_crop,sig1.samples[0].imag, label="sig1 imag",ls="--", marker="o")
# plt.plot(t2_crop,sig2.samples[0].real, label="sig2 real")
# plt.plot(t2_crop,sig2.samples[0].imag, label="sig2 imag",ls="--", marker="*")
# plt.ylabel("Amplitude")
# plt.xlabel("time [s]")
# plt.legend()
# plt.show()

# #### COMPENSATE FOR TIME DELAY (ZERO PADDING METHOD)
# #do a cross correlation, get back assumed lag in samples between both signals
# correlation = signal.correlate(sig2.samples[0], sig1.samples[0], mode="full")
# lags = signal.correlation_lags(sig1.samples[0].size, sig2.samples[0].size, mode="full")
# lag = lags[np.argmax(correlation)]

# time_offset = lag/sig1.sample_rate[0]

# samples_1_temp = np.zeros((len(sig1.samples[0])+abs(lag)),dtype=complex)
# samples_2_temp = np.zeros((len(sig2.samples[0])+abs(lag)),dtype=complex)

# #cut signal, only show matching area
# if lag >= 0:
#     samples_2_temp[0:-lag] = sig2.samples[0]
#     samples_1_temp[lag:] = sig1.samples[0]
# elif lag < 0:
#     samples_2_temp[abs(lag):] = sig2.samples[0]
#     samples_1_temp[0:-abs(lag)] = sig1.samples[0]

# sig1.samples[0] = samples_1_temp
# sig2.samples[0] = samples_2_temp

sig1.samples[0], sig2.samples[0], _ = comb_timedelay_compensation(sig1.samples[0], sig2.samples[0], sr=sig1.sample_rate[0], method="crop")

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
# C = sum(sig1.samples[0]*np.conj(sig2.samples[0]))

# print(np.angle(C))
# print(phase_offset)

# sig2.samples[0] = sig2.samples[0]*np.exp(1j*np.angle(C))

sig1.samples[0], sig2.samples[0] = comb_phase_compensation(sig1.samples[0], sig2.samples[0])

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