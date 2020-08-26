import time
import numpy as np
import scipy.signal as ssignal
import matplotlib.pyplot as plt
import sys
sys.path.append("..\\comm")
import comm as comm
#from scipy.signal.signaltools import wiener as wiener




###################### Tx ##################################
# signal parameters
LASER_LINEWIDTH = 0*200e3 # [Hz]
TX_UPSAMPLE_FACTOR = 5
EXPERIMENT = False
UPLOAD_SAMPLES = False
SNR = 200


# contruct signal
sig_tx = comm.signal.Signal(n_dims=1)
sig_tx.symbol_rate = 50e6 #50e6

# generate bits
sig_tx.generate_bits(n_bits=2**10)

# set constellation (modualtion format)
sig_tx.generate_constellation(order=4)
sig_tx.modulation_info = 'QPSK'

# create symbols
sig_tx.mapper()

# upsampling (to 5*50e6 --> 250 MSa/s)  and pulseshaping
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
# t = comm.utils.create_time_axis(sig_tx.sample_rate, np.size(sig_tx.samples[0]))

# sig_tx.plot_spectrum(0)
# upmixing to IF
sig_tx.samples[0] = sig_tx.samples[0] * np.exp(1j * 2 * np.pi * f_if * t)
sig_tx.center_frequency = f_if

# TODO: equalization of cosine MZM transfer function

# TODO: pre-equalization of AWG frequency response

# sig_tx.plot_spectrum(0)

# format samples so that driver can handle them (range: +-32768, integer, list)
maxVal = np.max(np.abs(np.concatenate((np.real(sig_tx.samples), np.imag(sig_tx.samples)))))
samples = np.asarray(sig_tx.samples) / maxVal
samples = np.concatenate((np.real(samples), np.imag(samples)))


############################################################################
################## Link ####################################################

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
    samples, phaseAcc, varPN = comm.channel.add_phase_noise(samples ,sig_tx.sample_rate[0] , LASER_LINEWIDTH, seed=None)
    sr = sig_tx.sample_rate[0]
    # plt.figure(1); plt.plot(phaseAcc); plt.show()
    
    # after heterodyne detection and balanced detection
    samples = np.real(samples)



#############################################################################
######################## Rx #################################################

#comm.visualizer.plot_spectrum(rx_samples, sample_rate=sr_rx)
# # comm.visualizer.plot_signal(samples, sample_rate=sr)

# resampling to N x symbolrate
sr_dsp = sig_tx.symbol_rate[0] * TX_UPSAMPLE_FACTOR
# # watch out, that this is really an integer, otherwise the samplerate is asynchronous with the data afterwards!!!
len_dsp = sr_dsp / sr * np.size(samples)
if len_dsp % 1:
    raise ValueError('DSP samplerate results in asynchronous sampling of the data symbols')
samples = ssignal.resample(samples, num=int(len_dsp), window=None)
sr = sr_dsp
# #comm.visualizer.plot_spectrum(rx_samples, sample_rate=sr)

# contruct rx signal
sig_rx = comm.signal.Signal(n_dims=1)
sig_rx.symbol_rate = sig_tx.symbol_rate
sig_rx.sample_rate = sr

# IQ-Downmixing and (ideal) lowpass filtering (real signal processing)
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

# IQ-Downmixing (ideal) (complex singal processing)
# samples_bb = samples *  np.exp(-1j*2*np.pi*(f_if+1e4*0)*t)
# sig_rx.samples[0] = samples_bb

#sig_rx.plot_spectrum()
#sig_rx.plot_constellation()

# "standard" coherent complex baseband signal processing
# Rx matched filter
sig_rx.raised_cosine_filter(roll_off=ROLL_OFF,root_raised=True) 
# TODO: compensate for the group delay of the filter???
# sig_rx.plot_eye()


# sampling phase adjustment
# results = comm.rx.sampling_phase_adjustment(sig_rx.samples[0], sample_rate=sig_rx.sample_rate[0], symbol_rate=sig_rx.symbol_rate[0])
# sig_rx.samples[0] = results['samples_out']
# sig_rx.plot_eye()
sig_rx.sampling_phase_adjustment()
# sig_rx.plot_eye()

# sampling
START_SAMPLE = 0
sps = sig_rx.sample_rate[0] / sig_rx.symbol_rate[0] # CHECK FOR INTEGER SPS!!!
rx_symbols = sig_rx.samples[0][START_SAMPLE::int(sps)]
#comm.visualizer.plot_constellation(rx_symbols)
# comm.visualizer.plot_constellation(rx_symbols)
# sig_rx.samples[0] = sig_rx.samples[0][START_SAMPLE::int(sps)]
# sig_rx.plot_constellation()


## implementing viterbi-viterbi for phase estimation
N_CPE = 21 # must be an odd number for Wiener filter
  # number of CPE_taps
mth_Power = 4 # 4=QAM&QPSK, 2=BPSK,...
phi_est = np.roll(comm.filters.moving_average(np.unwrap(mth_Power*np.angle(rx_symbols))/mth_Power, N_CPE, domain='freq'), -N_CPE//2)


## Wiener filter (see Ip2007, Acuna2013)
r = .3 # r = sigma²_phi / sigma²_ampl;  r>0; r is the ratio between the magnitude of the phase noise variance sigma²_phi and the additive noise variance sigma²
a = 1+r/2-np.sqrt((1+r/2)**2-1) # alpha
h_Wiener = a*r/(1-a**2) * a**np.arange(N_CPE//2+1) # postive half
h_Wiener = np.concatenate((np.flip(h_Wiener[1:]), h_Wiener)) # make symmetric
h_Wiener = h_Wiener / np.sum(h_Wiener) # normalize to unit sum (make unbiased estimator)
# plt.figure(2); plt.stem(np.arange(2*(N_CPE//2)+1)-N_CPE//2,h_Wiener,basefmt='C2-',use_line_collection=True); plt.show();
phi_est = np.roll(comm.filters.filter_samples(np.unwrap(mth_Power*np.angle(rx_symbols))/mth_Power, h_Wiener, domain='freq'), -N_CPE//2+1)
       
rx_symbols = rx_symbols * np.exp(-1j*(phi_est + np.pi/4))
rx_symbols = rx_symbols[1*N_CPE+1:-N_CPE*1] # crop start and end
comm.visualizer.plot_constellation(rx_symbols)