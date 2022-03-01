import time
import copy

import numpy as np
import scipy.signal as ssignal
import matplotlib.pyplot as plt
import scipy.interpolate as sinterp

import comm as comm


#%% Transmitter

#%%# signal parameters
LASER_LINEWIDTH = 1*600 # [Hz]
TX_UPSAMPLE_FACTOR = 5
EXPERIMENT = False
UPLOAD_SAMPLES = False
USE_PREDIST = True
SNR = 20

#%%# construct signal
sig_tx = comm.signal.Signal(n_dims=1)
sig_tx.symbol_rate = 50e6 

#%%# generate bits
sig_tx.generate_bits(n_bits=2**12, seed=1)

#%%# set constellation (modulation format)
sig_tx.generate_constellation(format='QAM', order=4)

#%%# create symbols
sig_tx.mapper()

#%%# upsampling and pulseshaping
ROLL_OFF = 0.1
sig_tx.pulseshaper(upsampling=TX_UPSAMPLE_FACTOR, pulseshape='rrc', roll_off=[ROLL_OFF])

#%%# generate DAC samples (analytical signalg at IF)
f_IF_nom = 1*30e6 #30e6
f_granularity = 1 / sig_tx.samples[0].size * sig_tx.sample_rate[0]
f_if = round(f_IF_nom / f_granularity) * f_granularity
print('intermediate frequency: {} MHz'.format(f_if/1e6))
t = np.arange(0, np.size(sig_tx.samples[0])) / sig_tx.sample_rate

#%%# upmixing to IF
sig_tx.samples[0] = sig_tx.samples[0] * np.exp(1j * 2 * np.pi * f_if * t)
sig_tx.center_frequency = f_if

# TODO: equalization of cosine MZM transfer function

#%%# pre-equalization of AWG frequency response
if USE_PREDIST:
    filtershape = np.load('setup_files/preDistFilter.npy')
    sig_tx.samples[0] = comm.filters.filter_arbitrary(sig_tx.samples[0], filtershape, sample_rate=sig_tx.symbol_rate[0]*TX_UPSAMPLE_FACTOR)

# format samples so that driver can handle them (range +-1)
maxVal = np.max(np.abs(np.concatenate((np.real(sig_tx.samples), np.imag(sig_tx.samples)))))
samples = np.asarray(sig_tx.samples) / maxVal
samples = np.concatenate((np.real(samples), np.imag(samples)))


#%% Link

#%%# Experiment 
if EXPERIMENT:
    if UPLOAD_SAMPLES:                    
        #%%## write samples to AWG
        comm.instrument_control.write_samples_AWG33522A(samples, ip_address='192.168.1.45',
                                                        sample_rate=[sig_tx.symbol_rate[0]*TX_UPSAMPLE_FACTOR]*2,
                                                        offset=[0.0, 0.0], amp_pp=[3.0]*2, channels=[1,2], 
                                                        out_filter=['normal']*2)
        time.sleep(0.3)
    #%%## get samples from scope
    sr, samples = comm.instrument_control.get_samples_DLM2034(channels=[1, 2], address='192.168.1.13')
    
    #%%## subtration of pos. and neg. detector
    samples = samples[0] - samples[1]
    
    # remove mean of signal
    samples = samples - np.mean(samples)

#%%# Simulation
else:
    samples = samples[0] + 1j*samples[1] # build ideal complex signal from Tx samples (no ampl. and phase noise)

    # =============================================================================
    sps = int(sig_tx.sample_rate[0] / sig_tx.symbol_rate[0])
    
    #%%## get samples from scope (repeat rx sequence)
    ext = 40000*sps + 4000*sps
    ratio_base = ext // samples.size
    ratio_rem = ext % samples.size        
    samples = np.concatenate((np.tile(samples, ratio_base), samples[:ratio_rem]), axis=0)
    
    #%%## add artificial delay 
    delay = 10*sps
    samples = samples[delay:]
    
    #%%## add phase ambiguity (phase rotation and conmplex conjugation)
    ### w/o conjugate
    # samples = samples * np.exp(1j*np.pi/3)
    ### w/ conjugate
    # ATTENTION: if conj is applied before linear phase rotation, sign of the
    # additional phase is flipped and subsequently "misinterpreted" (but compensated
    # correctly) by ambiguity compensation
    # samples = np.conj(samples * np.exp(-1j*np.pi/3))
    # =============================================================================
    
    #%%## add amplitude noise
    samples = comm.channel.set_snr(samples, snr_dB=SNR, sps=int(sig_tx.sample_rate[0]/sig_tx.symbol_rate[0]), seed=None)

    #%%## phase noise emulation
    samples = comm.channel.add_phase_noise(samples ,sig_tx.sample_rate[0] , LASER_LINEWIDTH, seed=1)['samples']
    sr = sig_tx.sample_rate[0]
        
    #%%## add artificial sample clock error
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

#%% Receiver

# contruct rx signal structure
sig_rx = copy.deepcopy(sig_tx)
sig_rx.samples = samples
sig_rx.sample_rate = sr

#%%# resampling to the same sample rate as at the transmitter
sr_dsp = sig_tx.sample_rate[0]

# # watch out, that this is really an integer, otherwise the samplerate is asynchronous with the data afterwards!!!
len_dsp = sr_dsp / sig_rx.sample_rate[0] * np.size(samples)
if len_dsp % 1:
    raise ValueError('DSP samplerate results in asynchronous sampling of the data symbols')
sig_rx.samples = ssignal.resample(sig_rx.samples[0], num=int(len_dsp), window=None)
sig_rx.sample_rate = sr_dsp
sig_rx.plot_spectrum(tit='received spectrum before IF downmixing')

#%%# IQ-Downmixing and (ideal) lowpass filtering
# ...either real signal processing
# t = np.arange(0, np.size(sig_rx.samples[0])) / sig_rx.sample_rate[0]
t = comm.utils.create_time_axis(sig_rx.sample_rate[0], np.size(sig_rx.samples[0]))
samples_r = sig_rx.samples[0] *  np.cos(2 * np.pi * f_if * t)

fc = sig_tx.symbol_rate[0] / 2 * (1 + ROLL_OFF) * 1.1 # cuttoff frequency of filter
fc = fc/(sig_rx.sample_rate[0]/2) # normalization to the sampling frequency
tmp = comm.filters.ideal_lp(samples_r, fc)
samples_r = tmp['samples_out']

samples_i = sig_rx.samples[0] *  np.sin(2 * np.pi * f_if * t)
# # comm.visualizer.plot_spectrum(samples_i, sample_rate=sr)
tmp = comm.filters.ideal_lp(samples_i, fc)
samples_i = tmp['samples_out']
sig_rx.samples[0] = samples_r - 1j * samples_i

# ... OR complex singal processing
# samples_bb = samples *  np.exp(-1j*2*np.pi*(f_if+1e4*0)*t)
# sig_rx.samples[0] = samples_bb

## From here: "standard" coherent complex baseband signal processing ############
#%%# resample to 2 sps
sps_new = 2
sps = sig_rx.sample_rate[0]/sig_rx.symbol_rate[0]
new_length = int(sig_rx.samples[0].size/sps*sps_new)
sig_rx.samples = ssignal.resample(sig_rx.samples[0], new_length, window='boxcar')
sig_rx.sample_rate = sps_new*sig_rx.symbol_rate[0]

#%%# normalize samples to mean magnitude of original constellation
mag_const = np.mean(abs(sig_rx.constellation[0]))
mag_samples = np.mean(abs(sig_rx.samples[0]))
sig_rx.samples = sig_rx.samples[0] * mag_const / mag_samples

sig_rx.plot_constellation(hist=True, tit='constellation before EQ')

adaptive_filter = True
#%%# equalizer
# either blind adaptive filter....
if adaptive_filter == True:    
    results = comm.rx.blind_adaptive_equalizer(sig_rx, n_taps=31, mu_cma=1e-3, 
                                               mu_rde=5e-3, mu_dde=0.5, decimate=False, 
                                               return_info=True, stop_adapting=-1, 
                                               start_rde=5000*0, start_dde=5000*0)
    
    sig_rx = results['sig']
    h = results['h'][0]
    eps = results['eps'][0]
    # plot error evolution
    plt.plot(np.abs(eps))
    plt.title('evolution of equalizer error')
    plt.xlabel('time / symbols')
    plt.ylabel('error /a.u.')
    plt.show()
    # plot last filter frequency response
    plt.plot(np.abs(np.fft.fftshift(np.fft.fft(h[-1,:]))))
    plt.show()            
    # plot evolution of filters frequency response
    plt.figure()
    ax = plt.subplot(projection='3d')
    f = np.fft.fftshift(np.fft.fftfreq(h[0,:].size, d=1/sig_rx.sample_rate[0]))
    outsymbs = [0, 1000, 5000, 10000, 20000, 30000, h[:,0].size-1]    
    for outsymb in outsymbs:
        plt.plot(f, np.ones(f.size)*outsymb, np.abs(np.fft.fftshift(np.fft.fft(h[int(outsymb),:]))))
    plt.title('evolution of equalizer frequency response')
    plt.xlabel('frequency / Hz')
    plt.ylabel('time / symbols')    
    plt.show()       
        
    # cut away init symbols
    sps = int(sig_rx.sample_rate[0]/sig_rx.symbol_rate[0])
    cut = 10000
    sig_rx.samples = sig_rx.samples[0][int(cut)*sps:]

# ... or matched filtering
else:
    # Rx matched filter
    sig_rx.raised_cosine_filter(roll_off=ROLL_OFF,root_raised=True) 
    # sig_rx.plot_eye()
    
    # crop samples here, if necessary
    sps = int(sig_rx.sample_rate[0] / sig_rx.symbol_rate[0])
    crop = 10*sps
    if crop != 0:
        sig_rx.samples = sig_rx.samples[0][crop:-crop]
    else:
        sig_rx.samples = sig_rx.samples[0]
    
    # sampling phase / clock adjustment
    BLOCK_SIZE = 1000 # size of one block in SYMBOLS... -1 for only one block
    sig_rx.sampling_clock_adjustment(BLOCK_SIZE)
    
#%%# sampling (if necessary)
START_SAMPLE = 0
sps = sig_rx.sample_rate[0] / sig_rx.symbol_rate[0] # CHECK FOR INTEGER SPS!!!
sig_rx.samples = sig_rx.samples[0][START_SAMPLE::int(sps)]
sig_rx.plot_constellation(0, hist=True, tit='constellation after EQ')

#%%# CPE
viterbi = False
# ...either VV
if viterbi:
    cpe_results = comm.rx.carrier_phase_estimation_VV(sig_rx.samples[0], n_taps=51, 
                                                      filter_shape='wiener', mth_power=4, 
                                                      rho=.05)
    sig_rx.samples = cpe_results['rec_symbols']
    est_phase = cpe_results['phi_est'].real
## ...or BPS
else:
    cpe_results = comm.rx.carrier_phase_estimation_bps(sig_rx.samples[0], sig_rx.constellation[0], 
                                               n_taps=15, n_test_phases=15, const_symmetry=np.pi/2)
    sig_rx.samples = cpe_results['samples_corrected']
    est_phase = cpe_results['est_phase_noise']
    
comm.visualizer.plot_signal(est_phase, tit='estimated phase noise')
sig_rx.plot_constellation(hist=True, tit='constellation after CPE')

#%%# delay and phase ambiguity estimation and compensation
sig_rx = comm.rx.symbol_sequence_sync(sig_rx, dimension=-1)
    
#%%# calc EVM
evm = comm.rx.calc_evm(sig_rx.samples[0], sig_rx.constellation[0], norm='max')
print("EVM: {:2.2%}".format(evm))

#%%# decision and demapper
sig_rx.decision()
sig_rx.demapper()

#%%# BER counting
ber_res = comm.rx.count_errors(sig_rx.bits[0], sig_rx.samples[0])
print('BER = {}'.format(ber_res['ber']))