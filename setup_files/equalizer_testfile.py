import time, copy, datetime

import numpy as np
import scipy.signal as ssignal
import matplotlib.pyplot as plt
import scipy.interpolate as sinterp

import skcomm as skc

#%% Tx 
#%% # signal parameters
LASER_LINEWIDTH = 1*100e3 # [Hz]
DAC_SR = 16e9
SNR = 20
F_OFFSET = 100e6 # frequency offset

#%% # contruct signal
sig_tx = skc.signal.Signal(n_dims=1)
sig_tx.symbol_rate = 12.8e9

TX_UPSAMPLE_FACTOR = DAC_SR / sig_tx.symbol_rate[0]

#%% # generate bits
sig_tx.generate_bits(n_bits=2**15, seed=1)

#%% # set constellation (modulation format)
sig_tx.generate_constellation(format='QAM', order=16)

#%% # create symbols
sig_tx.mapper()

#%% # upsampling and pulseshaping
ROLL_OFF = 0.1
sig_tx.pulseshaper(upsampling=TX_UPSAMPLE_FACTOR, pulseshape='rrc', roll_off=[ROLL_OFF])

# format samples so that driver can handle them (range +-1)
maxVal = np.max(np.abs(np.concatenate((np.real(sig_tx.samples), np.imag(sig_tx.samples)))))
samples = np.asarray(sig_tx.samples) / maxVal
samples = np.concatenate((np.real(samples), np.imag(samples)))

# build ideal complex signal from Tx samples (no ampl. and phase noise)
samples = samples[0] + 1j*samples[1] 

sps = sig_tx.sample_rate[0] / sig_tx.symbol_rate[0]

# get samples from scope (repeat rx sequence)
ext = 40000*sps + 4000*sps
ratio_base = int(ext // samples.size)
ratio_rem = int(ext % samples.size        )
samples = np.concatenate((np.tile(samples, ratio_base), samples[:ratio_rem]), axis=0)


#%% ## add amplitude noise
samples = skc.channel.set_snr(samples, snr_dB=SNR, sps=sig_tx.sample_rate[0]/sig_tx.symbol_rate[0], seed=None)

##%% ## add frequency offset
samples = skc.channel.add_frequency_offset(samples,
                                            sample_rate=sig_tx.sample_rate[0], 
                                            f_offset=F_OFFSET)

##%% ## phase noise emulation
samples = skc.channel.add_phase_noise(samples ,sig_tx.sample_rate[0] , LASER_LINEWIDTH, seed=1)['samples']
sr = sig_tx.sample_rate[0]

#%% ## add artificial sample clock error
ratio = 1.0 # ratio of sampling frequency missmatch     
n_old = np.size(samples, axis=0)
t_old = np.arange(n_old) / sr
n_new = int(np.round(ratio * n_old))
t_new = np.linspace(start=t_old[0], stop=t_old[-1], num=n_new, endpoint=True)
sr_new = 1 / (t_new[1] - t_new[0])
# interpolate signal at different timing / sampling instants
f = sinterp.interp1d(t_old, samples, kind='cubic')
samples = f(t_new)    


#%% Rx 

#%% # contruct rx signal structure
sig_rx = copy.deepcopy(sig_tx)
sig_rx.samples = samples
sig_rx.sample_rate = sr

sig_rx.plot_spectrum(tit='spectrum from scope')

#%% # From here: "standard" coherent complex baseband signal processing ############
#%% # resample to 2 sps
sps_new = 2
sps = sig_rx.sample_rate[0]/sig_rx.symbol_rate[0]
new_length = int(sig_rx.samples[0].size/sps*sps_new)
sig_rx.samples = ssignal.resample(sig_rx.samples[0], new_length, window='boxcar')
sig_rx.sample_rate = sps_new*sig_rx.symbol_rate[0]

sig_rx.plot_spectrum(tit='spectrum after resampling')

#%% # estimate SNR
sig_range = np.asarray([-1, 1])*sig_rx.symbol_rate[0]/2*(1+ROLL_OFF) + F_OFFSET
noise_range = np.asarray([-1.1, -1.05, 1.05, 1.1]) * sig_rx.symbol_rate[0]/2 * (1+ROLL_OFF) + F_OFFSET
                
spec = np.abs(np.fft.fftshift(np.fft.fft(sig_rx.samples[0])))**2
freq = np.fft.fftshift(np.fft.fftfreq(sig_rx.samples[0].size, 1/sig_rx.sample_rate[0]))
snr = skc.utils.estimate_snr_spectrum(freq, spec, sig_range=sig_range, 
                                       noise_range=noise_range, order=1, 
                                       noise_bw=sig_rx.symbol_rate[0], 
                                       scaling='lin', plotting=True)

print('est. SNR (from spectrum): {:.1f} dB'.format(snr))


#%% # frequency offset estimation / correction
results_foe = skc.rx.frequency_offset_estimation(sig_rx.samples[0], 
                                                  sample_rate=sig_rx.sample_rate[0],
                                                  order=4)

sig_rx.samples[0] = results_foe['samples_corrected']
print('estimated frequency offset: {:.0f} MHz'.format(results_foe['estimated_fo']/1e6))


#%% # normalize samples to mean magnitude of original constellation
mag_const = np.mean(abs(sig_rx.constellation[0]))
mag_samples= np.mean(abs(sig_rx.samples[0]))
sig_rx.samples = sig_rx.samples[0] * mag_const / mag_samples

sig_rx.plot_constellation(hist=True, tit='constellation before EQ')

adaptive_filter = True
#%% # either blind adaptive filter....
if adaptive_filter == True:   
    t0 = datetime.datetime.now()
    results_p = skc.rx.blind_adaptive_equalizer(copy.deepcopy(sig_rx), n_taps=51, mu_cma=1e-5, 
                                               mu_rde=1e-5, mu_dde=1e-4, decimate=True, 
                                               return_info=True, stop_adapting=-1, 
                                               start_rde=5000*1, start_dde=10000*1,
                                               compiled=False)
    dt_p = datetime.datetime.now() - t0
    print(f'Python EQ time {dt_p.total_seconds():f} s')
    
    t0 = datetime.datetime.now()
    results_c = skc.rx.blind_adaptive_equalizer(copy.deepcopy(sig_rx), n_taps=51, mu_cma=1e-5, 
                                               mu_rde=1e-5, mu_dde=1e-4, decimate=True, 
                                               return_info=True, stop_adapting=-1, 
                                               start_rde=5000*1, start_dde=10000*1,
                                               compiled=True)
    dt_c = datetime.datetime.now() - t0
    print(f'Cython EQ time {dt_c.total_seconds():f} s')
    print(f'Cython time gain {dt_p.total_seconds()/dt_c.total_seconds():f}')
    
    if np.allclose(results_p['sig'].samples[0], results_c['sig'].samples[0],rtol=1e-5, atol=1e-8):
        print('EQ results EQUAL')
    else:
        print('EQ results NOT EQUAL')
    
    sig_rx = results_p['sig']
    sig_rx.plot_constellation(0, hist=True, tit='constellation after Python EQ')
    
    sig_rx = results_c['sig']
    sig_rx.plot_constellation(0, hist=True, tit='constellation after Cython EQ')    
    
    h = results_c['h'][0]
    eps = results_c['eps'][0]
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
    cut = 10000
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
# sig_rx.plot_constellation(0, hist=True, tit='constellation after EQ')

#%% # CPE
viterbi = True
# ...either VV
if viterbi:
    cpe_results = skc.rx.carrier_phase_estimation_VV(sig_rx.samples[0], n_taps=31, 
                                                      filter_shape='wiener', mth_power=4, 
                                                      rho=.001)
    sig_rx.samples = cpe_results['rec_symbols']
    est_phase = cpe_results['phi_est'].real
# ...or BPS
else:
    cpe_results = skc.rx.carrier_phase_estimation_bps(sig_rx.samples[0], sig_rx.constellation[0], 
                                               n_taps=15, n_test_phases=45, const_symmetry=np.pi/2)
    sig_rx.samples = cpe_results['samples_corrected']
    est_phase = cpe_results['est_phase_noise']
    
# plt.plot(est_phase)
# plt.title('estimated phase noise')
# plt.xlabel('time / symbols')
# plt.ylabel('phase / rad')
# plt.grid()
# plt.show()

# sig_rx.plot_constellation(hist=True, tit='constellation after CPE')

#%% # delay and phase ambiguity estimation and compensation
sig_rx = skc.rx.symbol_sequence_sync(sig_rx, dimension=-1)
    
#%% # calc EVM
evm = skc.utils.calc_evm(sig_rx, norm='max')
print("EVM: {:2.2%}".format(evm[0]))

#%% # estimate SNR
snr = skc.utils.estimate_SNR_evm(sig_rx, norm='rms', method='data_aided', opt=False)
print("real SNR: {:.2f} dB, est. SNR (from EVM): {:.2f} dB".format(SNR, snr[0]))

#%% # decision and demapper
sig_rx.decision()
sig_rx.demapper()

#%% # BER counting
ber_res = skc.rx.count_errors(sig_rx.bits[0], sig_rx.samples[0])
print('BER = {:.2e}'.format(ber_res['ber']))