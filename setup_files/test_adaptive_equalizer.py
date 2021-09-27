import sys
import os
if not any(os.path.abspath('..') == p for p in sys.path): 
    print('adding comm module to path...')
    sys.path.insert(0, os.path.abspath('..'))
import time
import numpy as np
import scipy.signal as ssignal
import matplotlib.pyplot as plt
import scipy.interpolate as sinterp
import comm as comm
import copy
import time
#from scipy.signal.signaltools import wiener as wiener



############################################################
###################### Tx ##################################
############################################################

# signal parameters
LASER_LINEWIDTH = 0*5e2 # [Hz]
TX_UPSAMPLE_FACTOR = 5
SNR = 30

# contruct signal
sig_tx = comm.signal.Signal(n_dims=1)
sig_tx.symbol_rate = 50e6 #50e6

# generate bits
sig_tx.generate_bits(n_bits=2**12, seed=1)

# set constellation (modualtion format)
sig_tx.generate_constellation(order=16)
sig_tx.modulation_info = 'QAM'

# create symbols
sig_tx.mapper()

# upsampling and pulseshaping
ROLL_OFF = 0.1
sig_tx.pulseshaper(upsampling=TX_UPSAMPLE_FACTOR, pulseshape='rrc', roll_off=[ROLL_OFF])

# sig_tx.plot_eye()

# plot for checks
# sig_tx.plot_constellation(0)
#sig_tx.plot_eye(0)
#comm.visualizer.plot_signal(sig_tx.samples[0], sample_rate=sig_tx.sample_rate[0])

# generate DAC samples (analytical signalg at IF)
f_IF_nom = 1*30e6 #30e6
f_granularity = 1 / sig_tx.samples[0].size * sig_tx.sample_rate[0]
f_if = round(f_IF_nom / f_granularity) * f_granularity
print('intermediate frequency: {} MHz'.format(f_if/1e6))
t = np.arange(0, np.size(sig_tx.samples[0])) / sig_tx.sample_rate
# sig_tx.plot_spectrum(0)

# upmixing to IF
sig_tx.samples[0] = sig_tx.samples[0] * np.exp(1j * 2 * np.pi * f_if * t)
sig_tx.center_frequency = f_if

# TODO: equalization of cosine MZM transfer function

# sig_tx.plot_spectrum(0)

# format samples so that driver can handle them (range +-1)
maxVal = np.max(np.abs(np.concatenate((np.real(sig_tx.samples), np.imag(sig_tx.samples)))))
samples = np.asarray(sig_tx.samples) / maxVal
samples = np.concatenate((np.real(samples), np.imag(samples)))


############################################################################
################## Link ####################################################
############################################################################


samples = samples[0] + 1j*samples[1] # build ideal complex signal from Tx samples (no ampl. and phase noise)

# =============================================================================
sps = int(sig_tx.sample_rate[0] / sig_tx.symbol_rate[0])

# get samples from scope (repeat rx sequence)
ext = 80000*sps + 4000*sps
ratio_base = ext // samples.size
ratio_rem = ext % samples.size        
samples = np.concatenate((np.tile(samples, ratio_base), samples[:ratio_rem]), axis=0)

# add artificial delay 
delay = 10*sps
samples = samples[delay:]

## add noise
samples = comm.channel.set_snr(samples, snr_dB=SNR, sps=int(sig_tx.sample_rate[0]/sig_tx.symbol_rate[0]), seed=1)

## phase noise emulation
samples = comm.channel.add_phase_noise(samples ,sig_tx.sample_rate[0] , LASER_LINEWIDTH, seed=1)['samples']
sr = sig_tx.sample_rate[0]
# plt.figure(1); plt.plot(phaseAcc); plt.show()

# add artificial sample clock error
ratio = 1.0 # ratio of sampling frequency missmatch     
n_old = np.size(samples, axis=0)
t_old = np.arange(n_old) / sr
n_new = int(np.round(ratio * n_old))
t_new = np.linspace(start=t_old[0], stop=t_old[-1], num=n_new, endpoint=True)
sr_new = 1 / (t_new[1] - t_new[0])
# interpolate signal at different timing / sampling instants
f = sinterp.interp1d(t_old, samples, kind='cubic')
samples = f(t_new)

# after heterodyne detection and balanced detection
samples = np.real(samples)



#############################################################################
######################## Rx #################################################
#############################################################################
    
# contruct rx signal structure
sig_rx = copy.deepcopy(sig_tx)
sig_rx.samples = samples
sig_rx.sample_rate = sr

#comm.visualizer.plot_spectrum(rx_samples, sample_rate=sr_rx)
# # comm.visualizer.plot_signal(samples, sample_rate=sr)

# resampling to the same sample rate as at the transmitter
sr_dsp = sig_tx.sample_rate[0]

# # watch out, that this is really an integer, otherwise the samplerate is asynchronous with the data afterwards!!!
len_dsp = sr_dsp / sig_rx.sample_rate[0] * np.size(samples)
if len_dsp % 1:
    raise ValueError('DSP samplerate results in asynchronous sampling of the data symbols')
sig_rx.samples = ssignal.resample(sig_rx.samples[0], num=int(len_dsp), window=None)
sig_rx.sample_rate = sr_dsp
# #comm.visualizer.plot_spectrum(rx_samples, sample_rate=sr)
# sig_rx.plot_spectrum()


# IQ-Downmixing and (ideal) lowpass filtering
# ...either real signal processing
# t = np.arange(0, np.size(sig_rx.samples[0])) / sig_rx.sample_rate[0]
t = comm.utils.create_time_axis(sig_rx.sample_rate[0], np.size(sig_rx.samples[0]))
samples_r = sig_rx.samples[0] *  np.cos(2 * np.pi * f_if * t)
# comm.visualizer.plot_spectrum(samples_r, sample_rate=sr)
fc = sig_tx.symbol_rate[0] / 2 * (1 + ROLL_OFF) * 1.1 # cuttoff frequency of filter
fc = fc/(sig_rx.sample_rate[0]/2) # normalization to the sampling frequency
tmp = comm.filters.ideal_lp(samples_r, fc)
samples_r = tmp['samples_out']
# comm.visualizer.plot_spectrum(samples_r, sample_rate=sr)
samples_i = sig_rx.samples[0] *  np.sin(2 * np.pi * f_if * t)
# # comm.visualizer.plot_spectrum(samples_i, sample_rate=sr)
tmp = comm.filters.ideal_lp(samples_i, fc)
samples_i = tmp['samples_out']
# # comm.visualizer.plot_spectrum(samples_i, sample_rate=sr)
sig_rx.samples[0] = samples_r - 1j * samples_i

# ... OR complex singal processing
# samples_bb = samples *  np.exp(-1j*2*np.pi*(f_if+1e4*0)*t)
# sig_rx.samples[0] = samples_bb

#sig_rx.plot_spectrum()
#sig_rx.plot_constellation()

############# From here: "standard" coherent complex baseband signal processing ############
# resample to 2 sps
sps_new = 2
sps = sig_rx.sample_rate[0]/sig_rx.symbol_rate[0]
new_length = int(sig_rx.samples[0].size/sps*sps_new)
sig_rx.samples = ssignal.resample(sig_rx.samples[0], new_length, window='boxcar')
sig_rx.sample_rate = sps_new*sig_rx.symbol_rate[0]

# add artificial low pass filter
filtershape = np.asarray([[0, 0.0], [25e6, -10.0], [50e6, 2-0.0], [75e6, -30.0]])
sig_rx.samples = comm.filters.filter_arbitrary(sig_rx.samples[0], filtershape, sample_rate=sig_rx.sample_rate[0])

# copy first to second dim
sig_rx.n_dims = 1
sig_rx.samples = sig_rx.samples[0]
sig_rx.sample_rate = sig_rx.sample_rate[0]
sig_rx.symbol_rate = sig_rx.symbol_rate[0]
sig_rx.constellation = sig_rx.constellation[0]


start = time.time()
results = comm.rx.blind_adaptive_equalizer(sig_rx, n_taps=31, mu_cma=5e-3, mu_rde=5e-3, mu_dde=1e-3, 
                                            decimate=False, return_info=True, 
                                            stop_adapting=-1, start_rde=20000, start_dde=20000)
# results = comm.rx.blind_adaptive_equalizer(sig_rx, n_taps=[111, 55], mu=[1e-2, 5e-3], 
#                                             decimate=[True]*2, return_info=[True, True], 
#                                             stop_adapting=[-1, -1])
end = time.time()
print('equalizer took {:1.1f} s'.format(end - start))

sig_rx = results['sig']
h = results['h'][0]
eps = results['eps'][0]

plt.plot(np.abs(eps))
plt.show()

plt.plot(np.abs(np.fft.fftshift(np.fft.fft(h[-1,:]))))
plt.show()     

# h = results['h'][1]
# eps = results['eps'][1]     

# plt.plot(np.abs(eps))
# plt.show()

# plt.plot(np.abs(np.fft.fftshift(np.fft.fft(h[-1,:]))))
# plt.show() 

# take only one dimension for further processing
pick = 0
sig_rx.n_dims = 1
sig_rx.samples = sig_rx.samples[pick]
sig_rx.sample_rate = sig_rx.sample_rate[pick]
sig_rx.symbol_rate = sig_rx.symbol_rate[pick] 
sig_rx.constellation = sig_rx.constellation[pick]  
  
sps = int(sig_rx.sample_rate[0]/sig_rx.symbol_rate[0])
cut = 30000
# cut away init symbols and sample signal
sig_rx.samples = sig_rx.samples[0][int(cut)*sps::sps]
sig_rx.sample_rate = sig_rx.symbol_rate[0]

sig_rx.plot_constellation()
  
# CPE
cpe_results = comm.rx.carrier_phase_estimation_VV(sig_rx.samples[0], n_taps=31, filter_shape='wiener', mth_power=4, rho=.3)
sig_rx.samples = cpe_results['rec_symbols']
est_phase = cpe_results['phi_est']

# sig_rx.plot_constellation()

# delay and phase ambiguity estimation and compensation
sig_rx = comm.rx.symbol_sequence_sync(sig_rx, dimension=-1)
    
# plot constellation and calc BER
# sig_rx.plot_constellation()

# calc EVM
evm = comm.rx.calc_evm(sig_rx.samples[0], sig_rx.constellation[0], norm='max')
print("EVM: {:2.2%}".format(evm))

# decision and demapper
sig_rx.decision()
sig_rx.demapper()

# BER counting
ber_res = comm.rx.count_errors(sig_rx.bits[0], sig_rx.samples[0])
print('BER = {}'.format(ber_res['ber']))

# plt.plot(ber_res['err_idx'])