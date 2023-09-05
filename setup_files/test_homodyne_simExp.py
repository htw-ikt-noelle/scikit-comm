import time, copy

import numpy as np
import scipy.signal as ssignal
import matplotlib.pyplot as plt
import scipy.interpolate as sinterp

import skcomm as skc

#%% Tx 
#%% # signal parameters
LASER_LINEWIDTH = 1*100e3 # [Hz]
DAC_SR = 16e9
EXPERIMENT = False
UPLOAD_SAMPLES = False
HOLD_SHOT = False
USE_PREDIST = False
SINC_CORRECTION = True
SNR = 20
F_OFFSET = 100e6 # frequency offset

#%% # contruct signal
sig_tx = skc.signal.Signal(n_dims=1)
sig_tx.symbol_rate = 12.8e9

TX_UPSAMPLE_FACTOR = DAC_SR / sig_tx.symbol_rate[0]

#%% # generate bits
sig_tx.generate_bits(n_bits=2**15, seed=1)

#%% # set constellation (modulation format)
sig_tx.generate_constellation(format='QAM', order=4)

#%% # create symbols
sig_tx.mapper()

#%% # upsampling and pulseshaping
ROLL_OFF = 0.1
sig_tx.pulseshaper(upsampling=TX_UPSAMPLE_FACTOR, pulseshape='rrc', roll_off=[ROLL_OFF])

#%% # sinc correction
if SINC_CORRECTION:
    sig_tx.samples[0] = skc.pre_distortion.dac_sinc_correction(sig_tx.samples[0],
                                                                f_max=1.0)

#%% # pre-equalization of AWG frequency response
if USE_PREDIST:
    #filtershape = np.load('setup_files/preDistFilter.npy')
    filtershape = np.transpose(np.array([[-16e9, -8e9, -7.8e9, 0, 7.8e9, 8e9, 16e9], [-100, -100, 2, 0 , 2,-100,-100]]))
    sig_tx.samples[0] = skc.filters.filter_arbitrary(sig_tx.samples[0], 
                                                      filtershape, 
                                                      sample_rate=sig_tx.symbol_rate[0]*TX_UPSAMPLE_FACTOR)

# TODO: equalization of cosine MZM transfer function

# format samples so that driver can handle them (range +-1)
maxVal = np.max(np.abs(np.concatenate((np.real(sig_tx.samples), np.imag(sig_tx.samples)))))
samples = np.asarray(sig_tx.samples) / maxVal
samples = np.concatenate((np.real(samples), np.imag(samples)))


#%% Link 
#%% # Experiment 
if EXPERIMENT:
    if UPLOAD_SAMPLES:                    
        # write samples to AWG        
        skc.instrument_control.write_samples_Tektronix_AWG70002B(samples, ip_address='192.168.1.21', 
                                                        sample_rate=[DAC_SR], amp_pp=[0.25, 0.25], 
                                                        channels=[1, 2], log_mode = False)
        time.sleep(2.0)
        
    if not HOLD_SHOT:
    	# get samples from scope
        sr, samples = skc.instrument_control.get_samples_Tektronix_MSO6B(channels=[1, 2], 
                                                                      ip_address='192.168.1.20',
                                                                      number_of_bytes = 1,
                                                                      log_mode = False)
    	
        tmp_shot = copy.deepcopy(samples)
    else:
        samples = copy.deepcopy(tmp_shot) # see hint at https://github.com/spyder-ide/spyder/issues/11558
        # or, dump to and load from file:
        # see https://stackoverflow.com/questions/4530611/saving-and-loading-objects-and-using-pickle/4531859
        # import pickle
        # with open('tmp_shot.obj', 'wb') as ts:
        #   pickle.dump(tmp_shot, ts)
        # with open('tmp_shot.obj', 'rb') as ts:
        #   tmp_shot = pickle.load(ts)
    
    # # I-Q imbalance correction (for suppressing the DD-term for high SSI-suppression)
    # IQ_imbalance = -0.45  # [dB] (a pos. value means channel I is smaller than channel Q after coh. frontend, and v.v.)
    # samples[0] = samples[0]*10**(0.5*IQ_imbalance/20)
    # samples[1] = samples[1]*10**(-0.5*IQ_imbalance/20)
    
    # # remove mean of signal
    # samples = samples - np.mean(samples)
    
    # construct complex signal
    samples = samples[0] + 1j * samples[1]
    samples = samples - np.mean(samples)

#%% # Simulation 
else: # Simulation
    # build ideal complex signal from Tx samples (no ampl. and phase noise)
    samples = samples[0] + 1j*samples[1] 
    
    sps = sig_tx.sample_rate[0] / sig_tx.symbol_rate[0]
    
    # get samples from scope (repeat rx sequence)
    ext = 40000*sps + 4000*sps
    ratio_base = int(ext // samples.size)
    ratio_rem = int(ext % samples.size        )
    samples = np.concatenate((np.tile(samples, ratio_base), samples[:ratio_rem]), axis=0)
    
    # # add artificial delay 
    # delay = 10*sps
    # samples = samples[delay:]
    
    # add phase ambiguity (phase rotation and conmplex conjugation)
    # w/o conjugate...
    # samples = samples * np.exp(1j*np.pi/3)
    # ...w/ conjugate
    # ATTENTION: if conj is applied before linear phase rotation, sign of the
    # additional phase is flipped and subsequently "misinterpreted" (but compensated
    # correctly) by ambiguity compensation
    # samples = np.conj(samples * np.exp(-1j*np.pi/3))
    
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

sig_rx.plot_spectrum(tit='spectrum from scope',fNum=1)

#%% # From here: "standard" coherent complex baseband signal processing ############
#%% # resample to 2 sps
sps_new = 2
sps = sig_rx.sample_rate[0]/sig_rx.symbol_rate[0]
new_length = int(sig_rx.samples[0].size/sps*sps_new)
sig_rx.samples = ssignal.resample(sig_rx.samples[0], new_length, window='boxcar')
sig_rx.sample_rate = sps_new*sig_rx.symbol_rate[0]

sig_rx.plot_spectrum(tit='spectrum after resampling',fNum=2)

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

sig_rx.plot_constellation(hist=True, tit='constellation before EQ', fNum=3)

adaptive_filter = True
#%% # either blind adaptive filter....
if adaptive_filter == True:    
    results = skc.rx.blind_adaptive_equalizer(sig_rx, n_taps=51, mu_cma=1e-4, 
                                               mu_rde=1e-5, mu_dde=0.5, decimate=True, 
                                               return_info=True, stop_adapting=-1, 
                                               start_rde=5000*0, start_dde=5000*0,
                                               compiled=False)
    
    sig_rx = results['sig']
    h = results['h'][0]
    eps = results['eps'][0]
    # plot error evolution
    f = plt.figure(4)
    f.clear()
    plt.plot(np.abs(eps))
    plt.title('evolution of equalizer error')
    plt.xlabel('time / symbols')
    plt.ylabel('error /a.u.')
    plt.show()       
    # plot evolution of filters frequency response
    f = plt.figure(5)
    f.clear()
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
    f = plt.figure(6)
    f.clear()
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
sig_rx.plot_constellation(0, hist=True, tit='constellation after EQ',fNum=7)

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
                                               n_taps=15, n_test_phases=45, const_symmetry=np.pi/2,
                                               compiled=True)
    sig_rx.samples = cpe_results['samples_corrected']
    est_phase = cpe_results['est_phase_noise']

f = plt.figure(8)
f.clear()
plt.plot(est_phase)
plt.title('estimated phase noise')
plt.xlabel('time / symbols')
plt.ylabel('phase / rad')
plt.grid()
plt.show()

sig_rx.plot_constellation(hist=True, tit='constellation after CPE',fNum=9)

#%% # delay and phase ambiguity estimation and compensation
sig_rx = skc.rx.symbol_sequence_sync(sig_rx, dimension=-1)
    
#%% # calc EVM
evm = skc.utils.calc_evm(sig_rx, norm='max')
print("EVM: {:2.2%}".format(evm[0]))

#%% # estimate SNR
snr = skc.utils.estimate_SNR_evm(sig_rx, norm='rms', method='data_aided', opt=False)
if EXPERIMENT:
    print("est. SNR: {:.2f} dB".format(snr[0]))
else:
    print("real SNR: {:.2f} dB, est. SNR (from EVM): {:.2f} dB".format(SNR, snr[0]))

#%% # decision and demapper
sig_rx.decision()
sig_rx.demapper()

#%% # BER counting
ber_res = skc.rx.count_errors(sig_rx.bits[0], sig_rx.samples[0])
print('BER = {:.2e}'.format(ber_res['ber']))