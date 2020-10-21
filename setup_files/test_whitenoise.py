# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 12:14:04 2020

@author: UET-Labor
"""

import numpy as np

import matplotlib.pyplot as plt
import comm as comm
import numpy as np
import scipy.signal as scisig
from scipy.interpolate import interp1d


############################# test parameter ######################################
SAMPLE_RATE= 250e6 # AWG sample rate in Hz
N = 20000 # Number of complex noise samples
NperSEG = 64 # segment size for Welch spectral estimation

############################# Test signal erzeugen ################################
noise = np.random.randn(2*N).view(np.complex128)
# noise = np.random.uniform(low=-1, high=1, size=N*2).view(np.complex128)

plt.figure(1)
plt.plot(np.real(noise), np.imag(noise),'.'); plt.axis('equal'); plt.show()
#comm.visualizer.plot_spectrum(noise,SAMPLE_RATE)

spectrogram = scisig.welch(noise, fs=SAMPLE_RATE, window='hann', 
                                 nperseg=NperSEG, noverlap=None, nfft=None, detrend='constant', 
                                 return_onesided=False, scaling='spectrum', axis=- 1, average='mean')
freq = np.fft.fftshift(spectrogram[0])
PSD_IN = np.fft.fftshift(spectrogram[1])
PSD_IN_dB = 10*np.log10(PSD_IN/10)
plt.figure(2)
plt.plot(freq,PSD_IN_dB,'og-'); #plt.show()
############################# Ende Test signal erzeugen ################################


############################
# pre-filter (test)
FILTER = np.array([[0,35e6,55e6,65e6], [0,3.5,10.5,-10]],dtype='float').transpose()
#FILTER = np.array([[0,30e6,56e6,65e6], [0,4,12.5,-20]],dtype='float').transpose()
#plt.figure(1); plt.plot(FILTER[:,0],FILTER[:,1],'r-o'),plt.xlim((-sr/2,sr/2)); plt.show()
noise = comm.filters.filter_arbitrary(noise,FILTER,sample_rate=SAMPLE_RATE)


spectrogram = scisig.welch(noise, fs=SAMPLE_RATE, window='hann', 
                                 nperseg=NperSEG, noverlap=None, nfft=None, detrend='constant', 
                                 return_onesided=False, scaling='spectrum', axis=- 1, average='mean')

freq = np.fft.fftshift(spectrogram[0])
PSD_OUT = np.fft.fftshift(spectrogram[1])
PSD_OUT_dB = 10*np.log10(PSD_OUT/10)

plt.plot(freq,PSD_OUT_dB,'or-'); plt.show()


TF_dB = PSD_OUT_dB - PSD_IN_dB
plt.figure(2)
plt.plot(FILTER[:,0],FILTER[:,1], 'xr-')
plt.plot(freq,TF_dB,'ok-'); plt.show()


############################### send to AWG ########################################
# noise = noise[...,np.newaxis].transpose()
# noise = np.concatenate((np.real(noise), np.imag(noise)))

# plt.figure(2)
# plt.plot(noise[0],noise[1],'.r'); plt.axis('equal'); plt.show()

# comm.instrument_control.write_samples_AWG33522A(noise, ip_address='192.168.1.44',
#                                                         sample_rate=[SAMPLE_RATE]*2,
#                                                         offset=[0.0, 0.0], amp_pp=[3.0]*2, channels=[1,2], 
#                                                         out_filter=['normal']*2)



# comm.visualizer.plot_spectrum(np.real(noise))
# comm.visualizer.plot_spectrum(np.imag(noise))
# sr = 100
# FILTER = np.array([[-0.4*sr,-0.3*sr,0*sr,0.1*sr,0.2*sr], [0,-100,10,-150,-50]],dtype='float').transpose()

# plt.figure(1)
# plt.plot(FILTER[:,0],FILTER[:,1],'r-o'),plt.xlim((-sr/2,sr/2)); plt.show()


# samples_out = comm.filters.filter_arbitrary(noise,FILTER,sample_rate=sr)

# plt.figure(2)
# comm.visualizer.plot_spectrum(samples_out)



 