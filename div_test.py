import comm 
import numpy as np
import matplotlib.pyplot as plt

#### diversity gain test
sig = comm.signal.Signal(n_dims=1)
sig.symbol_rate = 5e9
sig.sample_rate = 15e9
x_pol_snr = np.random.randint(0,20)
y_pol_snr = np.random.randint(0,20)
n_bits = 2**16
df = sig.sample_rate[0]/n_bits
#### TX
sig.generate_bits(n_bits=n_bits,seed=None)
sig.generate_constellation(format='QAM',order=4)
sig.mapper()
sig.pulseshaper(upsampling=sig.sample_rate[0]/sig.symbol_rate[0],pulseshape='rrc',roll_off=.2)
#### CH
x_pol = comm.channel.set_snr(sig.samples[0],snr_dB=x_pol_snr,seed=1)
y_pol = comm.channel.set_snr(sig.samples[0],snr_dB=y_pol_snr,seed=1)
#### RX
# pre-combining SNR estimation
x_snr_est = comm.utils.estimate_snr_spectrum(np.linspace(-sig.sample_rate[0]/2+df/2,sig.sample_rate[0]/2-df/2,sig.samples[0].size), 
                                             np.abs(np.fft.fftshift(np.fft.fft(x_pol)))**2, 
                                             sig_range=np.array([-3e9,3e9]), 
                                             noise_range=np.array([-4e9,-3e9,3e9,4e9]),
                                             order=1,
                                             noise_bw=sig.sample_rate[0],
                                             plotting=False)

y_snr_est = comm.utils.estimate_snr_spectrum(np.linspace(-sig.sample_rate[0]/2+df/2,sig.sample_rate[0]/2-df/2,sig.samples[0].size), 
                                             np.abs(np.fft.fftshift(np.fft.fft(y_pol)))**2, 
                                             sig_range=np.array([-3e9,3e9]), 
                                             noise_range=np.array([-4e9,-3e9,3e9,4e9]),
                                             order=1,
                                             noise_bw=sig.sample_rate[0],
                                             plotting=False)
# EGC combining
sig_EGC = x_pol + y_pol
EGC_snr_est = comm.utils.estimate_snr_spectrum(np.linspace(-sig.sample_rate[0]/2+df/2,sig.sample_rate[0]/2-df/2,sig.samples[0].size), 
                                             np.abs(np.fft.fftshift(np.fft.fft(sig_EGC)))**2, 
                                             sig_range=np.array([-3e9,3e9]), 
                                             noise_range=np.array([-4e9,-3e9,3e9,4e9]),
                                             order=1,
                                             noise_bw=sig.sample_rate[0],
                                             plotting=False)
# MRC combining
sig_MRC = (10**(x_snr_est/10))*x_pol + (10**(y_snr_est/10))*y_pol
MRC_snr_est = comm.utils.estimate_snr_spectrum(np.linspace(-sig.sample_rate[0]/2+df/2,sig.sample_rate[0]/2-df/2,sig.samples[0].size), 
                                             np.abs(np.fft.fftshift(np.fft.fft(sig_MRC)))**2, 
                                             sig_range=np.array([-3e9,3e9]), 
                                             noise_range=np.array([-4e9,-3e9,3e9,4e9]),
                                             order=1,
                                             noise_bw=sig.sample_rate[0],
                                             plotting=False)

print('Sum of individual channel SNRs = {:.2f} dB.'.format(x_snr_est+y_snr_est))
print('MRC SNR = {:.2f} dB.'.format(MRC_snr_est))
print('MRC SNR gain = {:.2f} dB.'.format(MRC_snr_est-EGC_snr_est))