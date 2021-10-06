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
LASER_LINEWIDTH = 600 # [Hz]
TX_UPSAMPLE_FACTOR = 5
SNR = 30 
bps = 5

# contruct signal
sig_tx = comm.signal.Signal(n_dims=1)
sig_tx.symbol_rate = 50e6 #50e6

# generate bits
sig_tx.generate_bits(n_bits=2**12*bps, seed=1)

# generate constellation
sig_tx.generate_constellation(format='QAM', order=2**bps)

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
ext = 40000*sps + 4000*sps
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

# # add artificial low pass filter
# filtershape = np.asarray([[0, 0.0], [25e6, -10.0], [50e6, 2-0.0], [75e6, -30.0]])
# sig_rx.samples = comm.filters.filter_arbitrary(sig_rx.samples[0], filtershape, sample_rate=sig_rx.sample_rate[0])

# copy first to second dim
sig_rx.n_dims = 1
sig_rx.samples = sig_rx.samples[0]
sig_rx.sample_rate = sig_rx.sample_rate[0]
sig_rx.symbol_rate = sig_rx.symbol_rate[0]
sig_rx.constellation = sig_rx.constellation[0]


# matched filter
sig_rx.raised_cosine_filter(roll_off=ROLL_OFF,root_raised=True) 


sps = int(sig_rx.sample_rate[0]/sig_rx.symbol_rate[0])
cut = 500
# cut away boundary symbols and sample signal to 1 sps
sig_rx.samples = sig_rx.samples[0][int(cut)*sps:-(int(cut)*sps):sps]
sig_rx.sample_rate = sig_rx.symbol_rate[0]


# take only one dimension for further processing
pick = 0
sig_rx.n_dims = 1
sig_rx.samples = sig_rx.samples[pick]
sig_rx.sample_rate = sig_rx.sample_rate[pick]
sig_rx.symbol_rate = sig_rx.symbol_rate[pick] 
sig_rx.constellation = sig_rx.constellation[pick]  
  

sig_rx.plot_constellation()

  
# # VV CPE
# cpe_results = comm.rx.carrier_phase_estimation_VV(sig_rx.samples[0], n_taps=31, filter_shape='wiener', mth_power=4, rho=.3)
# sig_rx.samples = cpe_results['rec_symbols']
# est_phase = cpe_results['phi_est']

# sig_rx.samples[0] = sig_rx.samples[0] * np.exp(1j * np.pi/8)

# BPS CPE
cpe_results = comm.rx.carrier_phase_estimation_bps(sig_rx.samples[0], n_taps=15, n_test_phases=35, constellation=sig_rx.constellation[0])
samples_outs = cpe_results['samples_out']
est_phase_noise = cpe_results['est_phase_noise']
samples_corrected = cpe_results['samples_corrected']

sig_rx.samples[0] = samples_corrected

    
sig_rx.plot_constellation()

plt.hist2d(sig_rx.samples[0].real, sig_rx.samples[0].imag, bins=200, cmap=plt.get_cmap('jet'))
# plt.xlim(-5,5)
# plt.ylim(-5,5)
plt.axis('equal')
plt.show()


comm.visualizer.plot_signal(est_phase_noise)

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