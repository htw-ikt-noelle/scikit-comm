import time, copy

import numpy as np
import scipy.signal as ssignal
import matplotlib.pyplot as plt
import scipy.interpolate as sinterp
import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '..'))
import comm as comm

#########################################################################
#                            TX                                         #
#########################################################################

#%% # signal parameters
n_dims = 2
LASER_LINEWIDTH = 1*100e3 # [Hz]
DAC_SR = 16e9
ADC_MEM = 3e3
DAC_MEM = 1e3
EXPERIMENT = False
UPLOAD_SAMPLES = False
HOLD_SHOT = False
USE_PREDIST = False
SINC_CORRECTION = False
SNR = [60]*n_dims

#%% # contruct signal
sig_tx = comm.signal.Signal(n_dims=n_dims)
sig_tx.symbol_rate = 12.8e9

TX_UPSAMPLE_FACTOR = DAC_SR / sig_tx.symbol_rate[0]

#%% # generate bits
# fix for 4QAM right now
wanted_symbols = DAC_MEM/TX_UPSAMPLE_FACTOR
n_bits = 2*wanted_symbols
sig_tx.generate_bits(n_bits=2**15, seed=1)

#%% # set constellation (modulation format)
sig_tx.generate_constellation(format='QAM', order=4)

#%% # create symbols
sig_tx.mapper()

#%% # upsampling and pulseshaping
ROLL_OFF = 0.1
sig_tx.pulseshaper(upsampling=TX_UPSAMPLE_FACTOR, pulseshape='rrc', roll_off=ROLL_OFF)

#########################################################################
#                            CH                                         #
#########################################################################

#grab tx object
sig_ch = copy.deepcopy(sig_tx)

#concatenate fix for the moment
for dim in range(sig_ch.n_dims):
    sig_ch.samples[dim] = np.tile(sig_ch.samples[dim], 4)

# roll different sample vectors for timeshift in percent
if sig_ch.n_dims != 1:
    max_rnd_sample_diff_in_percent = 1
    rng = np.random.default_rng()
    diff_in_percent = rng.uniform(low=0,high=max_rnd_sample_diff_in_percent, size=((sig_ch.n_dims,)))
    lags = np.floor((ADC_MEM)/100*diff_in_percent)
    lags[0] = 0

    for dim in range(sig_ch.n_dims):
        sig_ch.samples[dim] = np.roll(sig_ch.samples[dim], int(lags[dim]))

#%% ## add amplitude noise
#for dim in range(sig_ch.n_dims):
#    sig_tx.samples[dim] = comm.channel.set_snr(sig_ch.samples[dim], snr_dB=SNR[dim], sps=sig_ch.sample_rate[dim]/sig_ch.symbol_rate[dim], seed=None)

##%% ## phase noise emulation
# TODO: call function with unique phase noise seeds per dimension, if different
# LOs are being used 
#for dim in range(sig_ch.n_dims):
#    sig_ch.samples[dim] = comm.channel.add_phase_noise(sig_ch.samples[dim] ,sig_ch.sample_rate[dim] , LASER_LINEWIDTH, seed=1)['samples']

#########################################################################
#                            RX                                         #
#########################################################################

#%% # contruct rx signal structure
sig_rx = copy.deepcopy(sig_tx)

## sampling error and ADC mem
sig_rx = comm.utils.add_sampling_error(sig_rx, ratio=1.0)
#sig_rx.samples[0] = sig_rx.samples[0][0:int(ADC_MEM)]
for dim in range(0,sig_rx.n_dims):
    sig_rx.samples[dim] = sig_rx.samples[dim][0:int(ADC_MEM)]

#sig_rx.plot_spectrum(tit='spectrum from scope')

#%% # From here: "standard" coherent complex baseband signal processing ############
#%% # resample to 2 sps
sig_rx = comm.utils.resample(sig_rx,target_sps=2)

#%% # normalize samples to mean magnitude of original constellation
#sig_rx = comm.utils.normalize_samples(sig_rx)

#### combining
# =============================================================================
# combining block goes here:

# 0. SNR estimation & shifting to samples[0] = Highest SNR, samples[1] = second highest SNR ... also needed for MRC comb!

# 1. Timedelay compensation
for i in range(1,sig_rx.n_dims):
    sig_rx.samples[0], sig_rx.samples[i], lag = comm.rx.comb_timedelay_compensation(sig_rx.samples[0], sig_rx.samples[i], method="crop", xcorr="abs")
    print("dim {}: lag: {}".format(i,lag))
    print(lags)

# 2. Phase comp. and combining
samples_rolling_sum = sig_rx.samples[0]
for i in range(1,sig_rx.n_dims):
    #phase compensation
    _, samples_second_sig_phase_align, _ = comm.rx.comb_phase_compensation(samples_rolling_sum, sig_rx.samples[i])
    #EGC combining
    samples_rolling_sum = sum(samples_rolling_sum, samples_second_sig_phase_align)

# 2.1 for continiously
for i in range(sig_rx.n_dims):
    sig_rx.samples[i] = samples_rolling_sum


# =============================================================================

sig_rx.plot_constellation(hist=True, tit='constellation before EQ')

adaptive_filter = False
#%% # either blind adaptive filter....
if adaptive_filter == True:    
    results = comm.rx.blind_adaptive_equalizer(sig_rx, n_taps=51, mu_cma=1e-4, 
                                               mu_rde=1e-5, mu_dde=0.5, decimate=True, 
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
     # plot last filter frequency response
    f = np.fft.fftshift(np.fft.fftfreq(h[0,:].size, d=1/sig_rx.sample_rate[0]))
    plt.plot(f, np.abs(np.fft.fftshift(np.fft.fft(h[-1,:]))))
    plt.title('last equalizer frequency response')
    plt.xlabel('frequency / Hz')
    plt.ylabel('amplitude a.u.')
    plt.show()               
        
    # cut away init symbols
    sps = int(sig_rx.sample_rate[0]/sig_rx.symbol_rate[0])
    cut = 5000
    sig_rx.samples = sig_rx.samples[0][int(cut)*sps:]

#%% # ... or matched filtering
else:
    # Rx matched filter
    sig_rx.raised_cosine_filter(roll_off=ROLL_OFF,root_raised=True) 
    
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
    
#%% # sampling (if necessary)
START_SAMPLE = 0
sps = sig_rx.sample_rate[0] / sig_rx.symbol_rate[0] # CHECK FOR INTEGER SPS!!!
sig_rx.samples = sig_rx.samples[0][START_SAMPLE::int(sps)]
sig_rx.plot_constellation(0, hist=True, tit='constellation after EQ')

#%% # CPE
viterbi = True
# ...either VV
if viterbi:
    cpe_results = comm.rx.carrier_phase_estimation_VV(sig_rx.samples[0], n_taps=31, 
                                                      filter_shape='wiener', mth_power=4, 
                                                      rho=.001)
    sig_rx.samples = cpe_results['rec_symbols']
    est_phase = cpe_results['phi_est'].real
# ...or BPS
else:
    cpe_results = comm.rx.carrier_phase_estimation_bps(sig_rx.samples[0], sig_rx.constellation[0], 
                                               n_taps=15, n_test_phases=45, const_symmetry=np.pi/2)
    sig_rx.samples = cpe_results['samples_corrected']
    est_phase = cpe_results['est_phase_noise']
    
plt.plot(est_phase)
plt.title('estimated phase noise')
plt.xlabel('time / symbols')
plt.ylabel('phase / rad')
plt.grid()
plt.show()

sig_rx.plot_constellation(hist=True, tit='constellation after CPE')

#%% # delay and phase ambiguity estimation and compensation
sig_rx = comm.rx.symbol_sequence_sync(sig_rx, dimension=-1)
    
#%% # calc EVM
evm = comm.utils.calc_evm(sig_rx, norm='max')
print("EVM: {:2.2%}".format(evm[0]))

#%% # estimate SNR
snr = comm.utils.estimate_SNR_evm(sig_rx, norm='rms', method='data_aided', opt=False)
if EXPERIMENT:
    print("est. SNR: {:.2f} dB".format(snr[0]))
else:
    print("real SNR: {:.2f} dB, est. SNR: {:.2f} dB".format(SNR[0], snr[0]))

#%% # decision and demapper
sig_rx.decision()
sig_rx.demapper()

#%% # BER counting
ber_res = comm.rx.count_errors(sig_rx.bits[0], sig_rx.samples[0])
print('BER = {}'.format(ber_res['ber']))