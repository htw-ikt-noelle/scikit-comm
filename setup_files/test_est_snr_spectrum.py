import math

import numpy as np
from scipy import optimize
import scipy.interpolate as sinter
import scipy.signal as ssignal
import scipy.special as sspecial
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial

import comm



sig = comm.signal.Signal(n_dims=1)
sig.symbol_rate = 12.8e9

TX_UPSAMPLE_FACTOR = 2#DAC_SR / sig_tx.symbol_rate[0]

sig.generate_bits(n_bits=2**15, seed=1)

sig.generate_constellation(format='QAM', order=4)

sig.mapper()

ROLL_OFF = 0.1
sig.pulseshaper(upsampling=TX_UPSAMPLE_FACTOR, pulseshape='rrc', roll_off=[ROLL_OFF])

# sig.plot_spectrum(tit='noise free')

sig.set_snr(10)

# sig.plot_spectrum(tit='noisy')

y = abs(np.fft.fftshift(np.fft.fft(sig.samples[0])))**2 
x = np.fft.fftshift(np.fft.fftfreq(y.size, d=1/sig.sample_rate[0])) 

# spectrogram_in = ssignal.welch(sig.samples[0], fs=sig.sample_rate[0], window='hann', 
#                               nperseg=128, noverlap=None, nfft=None, detrend=False, 
#                               return_onesided=False, scaling='spectrum', axis=- 1, average='mean')

# # shifting the zero-frequency component to the center of the spectrum
# x = np.fft.fftshift(spectrogram_in[0]) 
# # magnitude spectrum of input signal
# y = np.fft.fftshift((spectrogram_in[1]))

# dx = 1e9
# x = np.arange(193.0e12, 193.5e12+dx, dx)

noise_bw = sig.symbol_rate[0]

# # linear!!!!
# rng = np.random.default_rng()
# y = ssignal.windows.gaussian(x.size, 30) + 0*np.linspace(0, 0.2, x.size) + 0.2 + rng.normal(scale=0.02,size=x.size)

sig_range = np.asarray([-7e9, 7e9])
noise_range = np.asarray([-12e9, -8e9, 8e9, 12e9])

order = 3
##############################################

comm.utils.estimate_snr_spectrum(x, y, sig_range, noise_range, order=order, noise_bw=noise_bw,scaling='lin',plotting=True)