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
#from scipy.signal.signaltools import wiener as wiener



############################################################
###################### Tx ##################################
############################################################

# signal parameters
LASER_LINEWIDTH = 0*1e3 # [Hz]
TX_UPSAMPLE_FACTOR = 5
EXPERIMENT = True
UPLOAD_SAMPLES = False
USE_PREDIST = True
SNR = 200

# contruct signal
sig_tx = comm.signal.Signal(n_dims=1)
sig_tx.symbol_rate = 50e6 #50e6

# generate bits
sig_tx.generate_bits(n_bits=2**12, seed=1)

# set constellation (modualtion format)
sig_tx.generate_constellation(order=4)
sig_tx.modulation_info = 'QPSK'

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

# TODO: pre-equalization of AWG frequency response
if USE_PREDIST:
    filtershape = np.load('preDistFilter.npy')
    sig_tx.samples[0] = comm.filters.filter_arbitrary(sig_tx.samples[0], filtershape, sample_rate=sig_tx.symbol_rate[0]*TX_UPSAMPLE_FACTOR)

# sig_tx.plot_spectrum(0)

# format samples so that driver can handle them (range +-1)
maxVal = np.max(np.abs(np.concatenate((np.real(sig_tx.samples), np.imag(sig_tx.samples)))))
samples = np.asarray(sig_tx.samples) / maxVal
samples = np.concatenate((np.real(samples), np.imag(samples)))


############################################################################
################## Link ####################################################
############################################################################

##################### Experiment ###########################################
if EXPERIMENT:
    if UPLOAD_SAMPLES:                    
        # write samples to AWG
        comm.instrument_control.write_samples_AWG33522A(samples, ip_address='192.168.1.45',
                                                        sample_rate=[sig_tx.symbol_rate[0]*TX_UPSAMPLE_FACTOR]*2,
                                                        offset=[0.0, 0.0], amp_pp=[3.0]*2, channels=[1,2], 
                                                        out_filter=['normal']*2)
        time.sleep(0.3)
    # get samples from scope
    sr, samples = comm.instrument_control.get_samples_DLM2034(channels=[1, 2], address='192.168.1.13')
    
    # subtration of pos. and neg. detector
    samples = samples[0] - samples[1]
    
    #TODO: maybe remove mean of signal
    samples = samples - np.mean(samples)

###################### Simulation ###########################################
else:
    samples = samples[0] + 1j*samples[1] # build ideal complex signal from Tx samples (no ampl. and phase noise)

    # =============================================================================
    sps = int(sig_tx.sample_rate[0] / sig_tx.symbol_rate[0])
    
    # get samples from scope (repeat rx sequence)
    ext = 8192*sps + 4000*sps
    ratio_base = ext // samples.size
    ratio_rem = ext % samples.size        
    samples = np.concatenate((np.tile(samples, ratio_base), samples[:ratio_rem]), axis=0)
    
    # add artificial delay 
    delay = 10*sps
    samples = samples[delay:]
    
    # add phase ambiguity (phase rotation and conmplex conjugation)
    ### w/o conjugate
    # samples = samples * np.exp(1j*np.pi/3)
    ### w/ conjugate
    # ATTENTION: if conj is applied before linear phase rotation, sign of the
    # additional phase is flipped and subsequently "misinterpreted" (but compensated
    # correctly) by ambiguity compensation
    # samples = np.conj(samples * np.exp(-1j*np.pi/3))
    # =============================================================================
    
    ## add amplitude noise
    samples = comm.channel.set_snr(samples, snr_dB=SNR, sps=int(sig_tx.sample_rate[0]/sig_tx.symbol_rate[0]), seed=None)

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
sig_rx.plot_spectrum()

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


adaptive_filter = True
# blind adaptive filter....
if adaptive_filter == True:
    # blind adaptive equalizer
    # see [1] D. Godard, “Self-recovering equalization and carrier tracking in twodimensional data communication systems,” IEEE Trans. Commun., vol. 28, no. 11, pp. 1867–1875, Nov. 1980.
    # and [2] S. Savory, "Digital Coherent Optical Receivers: Algorithms and Subsystems", IEEE STQE, vol 16, no. 5, 2010
    
    # length of filter impuse resoponse
    n_taps = 111 # has to be odd
    sps = int(sig_rx.sample_rate[0] / sig_rx.symbol_rate[0])
    # step size for stochastic gradient method
    mu = 2e-2
    # init equalizer impulse response to delta
    h = np.zeros(n_taps, dtype=np.complex128)
    h[n_taps//2] = 1.0
    h_tmp = [h]
    eps_tmp = []
    
    samples_in = sig_rx.samples[0]
    samples_out = np.full(samples_in.size-n_taps, np.nan, dtype=np.complex128)
    
    # desired modulus for p=2, see [1], eq. (28) + 1
    r = np.mean(np.abs(samples_in)**4) / np.mean(np.abs(samples_in)**2)
    
    # comm.visualizer.plot_eye(samples_in[-500:], sample_rate = sig_rx.sample_rate[0], bit_rate = sig_rx.symbol_rate[0])
    cut = 10e3
    # comm.visualizer.plot_constellation(samples_in[cut*sps:-cut*sps:sps])
    
    for sample in range(0, samples_out.size, 1):
        
        # filter the signal for each output sample (convolution)
        # see [1], eq. (5)
        samples_out[sample] = np.sum(h * samples_in[n_taps+sample:sample:-1])
        
        # for each symbol, calculate error signal and update impulse response
        if (sample % sps == 0):
            # print('calc error \n')
            # see [1], eq. (26)
            eps = samples_out[sample] * (np.abs(samples_out[sample])**2 - r)
            eps_tmp.append(eps)
            # see [1], eq (28)
            h -= mu * np.conj(samples_in[n_taps+sample:sample:-1]) * eps
            h_tmp.append(h)        
        
        
    plt.plot(np.abs(np.asarray(eps_tmp)))
    plt.show()
    
    plt.plot(np.abs(np.fft.fftshift(np.fft.fft(np.asarray(h_tmp[-1])))))
    plt.show()
            
    # comm.visualizer.plot_eye(samples_out[-500:], sample_rate = sig_rx.sample_rate[0], bit_rate = sig_rx.symbol_rate[0]) 
    # comm.visualizer.plot_constellation(samples_out[cut*sps:-cut*sps:sps])    

    # cut away first samples
    sig_rx.samples = samples_out[int(cut)*sps:]
    

# matched filtering
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
    BLOCK_SIZE = -1 # size of one block in SYMBOLS... -1 for only one block
    sig_rx.sampling_clock_adjustment(BLOCK_SIZE)
    
# sampling
START_SAMPLE = 0
sps = sig_rx.sample_rate[0] / sig_rx.symbol_rate[0] # CHECK FOR INTEGER SPS!!!
sig_rx.samples = sig_rx.samples[0][START_SAMPLE::int(sps)]
sig_rx.plot_constellation(0)


# CPE
cpe_results = comm.rx.carrier_phase_estimation_VV(sig_rx.samples[0], n_taps=31, filter_shape='wiener', mth_power=4, rho=.3)
sig_rx.samples = cpe_results['rec_symbols']
est_phase = cpe_results['phi_est']

sig_rx.plot_constellation()

# delay and phase ambiguity estimation and compensation
sig_rx = comm.rx.symbol_sequence_sync(sig_rx, dimension=-1)
    
# plot constellation and calc BER
sig_rx.plot_constellation()

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