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

order = 1
##############################################

sig_range.sort()
sig_range = np.expand_dims(sig_range, axis=-1)
sig_range_idx = np.argmin(np.abs(x-sig_range), axis=-1)
sig_range=np.squeeze(sig_range)

noise_range.sort()

noise_range = np.expand_dims(noise_range, axis=-1)
noise_range_idx = np.argmin(np.abs(x-noise_range), axis=-1)
noise_range = np.squeeze(noise_range)

# noise_range_lft = noise_range[:2]
# noise_range_lft = np.expand_dims(noise_range_lft, axis=-1)
# noise_range_lft_idx = np.argmin(np.abs(x-noise_range_lft), axis=-1)
# noise_range_lft = np.squeeze(noise_range_lft)

# noise_range_rgt = noise_range[2:]
# noise_range_rgt = np.expand_dims(noise_range_rgt, axis=-1)
# noise_range_rgt_idx = np.argmin(np.abs(x-noise_range_rgt), axis=-1)
# noise_range_rgt = np.squeeze(noise_range_rgt)

y_sig_n2 = y[sig_range_idx[0]:sig_range_idx[1]]
# dx_sig = np.diff(x[sig_range_idx[0]:sig_range_idx[1]])
# p_sig_n2 = np.sum(y_sig_n2[:-1]*dx_sig)
p_sig_n2 = np.trapz(y_sig_n2, x=x[sig_range_idx[0]:sig_range_idx[1]])


x_n = np.append(x[noise_range_idx[0]:noise_range_idx[1]], x[noise_range_idx[2]:noise_range_idx[3]])
y_n = np.append(y[noise_range_idx[0]:noise_range_idx[1]], y[noise_range_idx[2]:noise_range_idx[3]])
c = Polynomial.fit(x_n, y_n, order)
xx_n, yy_n = c.linspace(n=x.size, domain=[x[0], x[-1]])

C = c.integ()

x_sig_mean = np.mean(sig_range)
p_n1 = C(x_sig_mean + noise_bw/2) - C(x_sig_mean - noise_bw/2)

p_n2 = C(x[sig_range_idx[1]]) - C(x[sig_range_idx[0]])

# print(p_n2/p_n1)

snr = (p_sig_n2 - p_n2) / p_n1
snr_db = 10 * np.log10(snr)

print('est. SNR = {:.1f} dB'.format(snr_db))



# plt.plot(x_n,y_n,'o')
# plt.plot(xx_n,yy_n,'r-')
# plt.show()

plt.plot(x,y)
plt.plot(x[sig_range_idx], y[sig_range_idx], 'o')
plt.plot(x[noise_range_idx[0]:noise_range_idx[1]], y[noise_range_idx[0]:noise_range_idx[1]], 'r--')
plt.plot(x[noise_range_idx[2]:noise_range_idx[3]], y[noise_range_idx[2]:noise_range_idx[3]], 'r--')
plt.plot(xx_n,yy_n,'g-')

plt.figure()
plt.plot(x,10*np.log10(y))
plt.plot(x[sig_range_idx], 10*np.log10(y[sig_range_idx]), 'o')
plt.plot(x[noise_range_idx[0]:noise_range_idx[1]], 10*np.log10(y[noise_range_idx[0]:noise_range_idx[1]]), 'r--')
plt.plot(x[noise_range_idx[2]:noise_range_idx[3]], 10*np.log10(y[noise_range_idx[2]:noise_range_idx[3]]), 'r--')
plt.plot(xx_n,10*np.log10(yy_n),'g-')
