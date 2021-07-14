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
#from scipy.signal.signaltools import wiener as wiener



############################################################
###################### Tx ##################################
############################################################

# signal parameters
LASER_LINEWIDTH = 0*1e3 # [Hz]
TX_UPSAMPLE_FACTOR = 5
EXPERIMENT = False
UPLOAD_SAMPLES = True
USE_PREDIST = False
SNR = 20

# contruct signal
sig_tx = comm.signal.Signal(n_dims=1)
sig_tx.symbol_rate = 50e6 #50e6

# generate bits
sig_tx.generate_bits(n_bits=2**14)

# set constellation (modualtion format)
sig_tx.generate_constellation(order=4)
sig_tx.modulation_info = 'QPSK'

# create symbols
sig_tx.mapper()

# upsampling and pulseshaping
ROLL_OFF = 0.1
sig_tx.pulseshaper(upsampling=TX_UPSAMPLE_FACTOR, pulseshape='rrc', roll_off=[ROLL_OFF])

# sig_tx.plot_eye()
# TODO: compensate for the group delay of RRC filter??

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

sig_tx.plot_spectrum(0)
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

###################### Simulation ###########################################
else:
    samples = samples[0] + 1j*samples[1] # build ideal complex signal from Tx samples (no ampl. and phase noise)

    ## add amplitude noise
    samples = comm.channel.set_snr(samples, snr_dB=SNR, sps=int(sig_tx.sample_rate[0]/sig_tx.symbol_rate[0]), seed=None)

    ## phase noise emulation
    samples, phaseAcc, varPN = comm.channel.add_phase_noise(samples ,sig_tx.sample_rate[0] , LASER_LINEWIDTH, seed=1)
    sr = sig_tx.sample_rate[0]
    # plt.figure(1); plt.plot(phaseAcc); plt.show()
    
    # add artificial sample clock error
    ratio = 1.00 # ratio of sampling frequency missmatch     
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

#comm.visualizer.plot_spectrum(rx_samples, sample_rate=sr_rx)
# # comm.visualizer.plot_signal(samples, sample_rate=sr)

# resampling to the same sample rate as at the transmitter
sr_dsp = sig_tx.symbol_rate[0] * TX_UPSAMPLE_FACTOR

# # watch out, that this is really an integer, otherwise the samplerate is asynchronous with the data afterwards!!!
len_dsp = sr_dsp / sr * np.size(samples)
if len_dsp % 1:
    raise ValueError('DSP samplerate results in asynchronous sampling of the data symbols')
samples = ssignal.resample(samples, num=int(len_dsp), window=None)
sr = sr_dsp
# #comm.visualizer.plot_spectrum(rx_samples, sample_rate=sr)

# contruct rx signal structure
sig_rx = comm.signal.Signal(n_dims=1)
sig_rx.symbol_rate = sig_tx.symbol_rate
sig_rx.sample_rate = sr

# IQ-Downmixing and (ideal) lowpass filtering
# ...either real signal processing
t = np.arange(0, np.size(samples)) / sr
t = comm.utils.create_time_axis(sr, np.size(samples))
samples_r = samples *  np.cos(2 * np.pi * f_if * t)
# comm.visualizer.plot_spectrum(samples_r, sample_rate=sr)
fc = sig_tx.symbol_rate[0] / 2 * (1 + ROLL_OFF) * 1.1 # cuttoff frequency of filter
fc = fc/(sr/2) # normalization to the sampling frequency
tmp = comm.filters.ideal_lp(samples_r, fc)
samples_r = tmp['samples_out']
# comm.visualizer.plot_spectrum(samples_r, sample_rate=sr)
samples_i = samples *  np.sin(2 * np.pi * f_if * t)
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
# Rx matched filter
sig_rx.raised_cosine_filter(roll_off=ROLL_OFF,root_raised=True) 
# TODO: compensate for the group delay of the filter???
# sig_rx.plot_eye()

# sampling phase / clock adjustment
BLOCK_SIZE = -1 # size of one block in SYMBOLS... -1 for only one block
sig_rx.sampling_clock_adjustment(BLOCK_SIZE)
# samples = sig_rx.samples[0]
# results = comm.rx.sampling_clock_adjustment(samples, sample_rate=sr, 
#                                             symbol_rate=sig_tx.symbol_rate[0], 
#                                             block_size=BLOCK_SIZE)
# samples = results['samples_out']
# shifts = results['est_shift']
# # plt.stem(shifts),plt.show()
# sig_rx.samples[0] = samples

# sampling
START_SAMPLE = 0
sps = sig_rx.sample_rate[0] / sig_rx.symbol_rate[0] # CHECK FOR INTEGER SPS!!!
rx_symbols = sig_rx.samples[0][START_SAMPLE::int(sps)]
comm.visualizer.plot_constellation(rx_symbols)

# CPE
cpe_results = comm.rx.carrier_phase_estimation_VV(rx_symbols, n_taps=21, filter_shape='wiener', mth_power=4, rho=.3)
rx_symbols = cpe_results['rec_symbols']
est_phase = cpe_results['phi_est']

# calc EVM
evm = comm.rx.calc_evm(rx_symbols, sig_tx.constellation[0], norm='max')
print("EVM: {:2.2%}".format(evm))


comm.visualizer.plot_signal(est_phase)
comm.visualizer.plot_constellation(rx_symbols)
# comm.visualizer.plot_signal(abs(rx_symbols))