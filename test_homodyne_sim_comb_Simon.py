import copy
import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '..'))
import comm as comm
from scipy.special import erfc

#########################################################################
#                            special functions                          #
#########################################################################

#%% helper functions
def gen_SIMO_samples(sig, max_phase_offset_in_rad=np.pi/3, max_timedelay_in_percent=10, n_apertures=2, repeat=5, cut_to=3, seed=None, subsample_shift=True):
    """
    NOTE: This function is just a temporary way to generate signals (close to the experimental setup as possible) on signal object level. 

    generate some "realistic" SIMO signals 
    """
    len_orig = len(sig.samples[0])
    sig.samples[0] = np.tile(sig.samples[0],reps=repeat)
    len_vectors = len(sig.samples[0])
    
    rng = np.random.default_rng(seed=seed)
    time_delay = rng.uniform(0,max_timedelay_in_percent,size=n_apertures) #[%]
    time_delay_in_samples = np.around(time_delay*(len_orig/100),0) #[samples]
    time_delay_in_subsamples = time_delay*(len_orig/100)-time_delay_in_samples #[samples]

    phase_offset = rng.uniform(0,max_phase_offset_in_rad,size=n_apertures)
            
    time_delay_in_samples[0] = 0 
    time_delay_in_subsamples[0] = 0
    phase_offset[0] = 0
    
    if subsample_shift != True:
        time_delay_in_subsamples[:] = 0

    samples_to_proceed = np.zeros((n_apertures, int(len_vectors+max(time_delay_in_samples))), dtype=complex)
    for n in range(n_apertures):
        samples_to_proceed[n,0:len_vectors] = sig.samples[0]
        samples_to_proceed[n,:] = samples_to_proceed[n,:]*np.exp(1j*phase_offset[n])
        samples_to_proceed[n,:] = np.roll(samples_to_proceed[n,:], shift=int(time_delay_in_samples[n]))
        samples_to_proceed[n,:] = comm.filters.time_shift(samples_to_proceed[n,:], sample_rate=sig.sample_rate[0],tau=time_delay_in_subsamples[n]/sig.symbol_rate[0])

    #crop signals to cut_to
    return_samples = np.zeros((n_apertures, len_orig*cut_to), dtype=complex)
    for n in range(n_apertures):
        return_samples[n,:] = samples_to_proceed[n,(int(max(time_delay_in_samples))):(int(len_orig*cut_to+max(time_delay_in_samples)))]
        sig.samples[n] = return_samples[n,:] 

    return_dict = {
        "time_delay_in_samples": time_delay_in_samples,
        "time_delay_in_subsamples": time_delay_in_subsamples,
        "phase_offset": np.rad2deg(phase_offset)
    }

    return sig, return_dict, return_samples

def plot_signal_timebase(sig1, tit=""):
    """
    NOTE: This function is just a temporary function to plot complex signals comparable in the time domain.
    """
    t = np.arange(0, (len(sig1.samples[0])/sig1.sample_rate[0]),1/sig1.sample_rate[0])
    plt.figure()
    plt.title(tit+str(", constellation"))
    plt.scatter(sig1.samples[0].real, sig1.samples[0].imag, label="Samples 1")
    plt.scatter(sig1.samples[1].real, sig1.samples[1].imag, label="Samples 2")
    plt.legend()
    plt.grid()
    plt.show()

    plt.figure()
    plt.title(tit+str(", timebase"))
    plt.plot(t,sig1.samples[0].real, label="samples 1 real")
    plt.plot(t,sig1.samples[0].imag, label="samples 1 imag",ls="--", marker="o")
    plt.plot(t,sig1.samples[1].real, label="samples 2 real")
    plt.plot(t,sig1.samples[1].imag, label="samples 2 imag",ls="--", marker="*")
    plt.xlim(0,20/sig1.sample_rate[0])
    plt.ylabel("Amplitude")
    plt.xlabel("time [s]")
    plt.legend()
    plt.show()

#%% MC loop

#### parameter setup
MC = 10
SNR_vec = np.arange(0,11)
distribution = 'rayleigh'

#### signal parameters
n_dims = 2
amount_of_symbols = 2**12
mod_format = "QAM"
mod_order = 4
ROLL_OFF = 0.1 #rrc
plotting = False

comb_method = "MRC"
adaptive_filter = False
LASER_LINEWIDTH = 0*100e3

#### init output array
BER_vec = np.full_like(SNR_vec,0.,dtype='float')
SNR_comb_vec = np.full_like(SNR_vec,0.,dtype='float')
mean_set_SNR_vec = np.full_like(SNR_vec,0.,dtype='float')

#%%% loop over SNR val
for SNR_idx, SNR_val in enumerate(SNR_vec):
    
    # init output arrays
    BER_tmp = np.zeros((MC,),dtype='float')
    SNR_comb_tmp = np.zeros((MC,),dtype='float')
    set_SNR_vec = np.zeros((MC,),dtype='float')
    
    # gen amplitude coefficients
    mean_snr_lin = 10**(SNR_val/10)
    rng = np.random.default_rng(seed=None)
    if distribution == 'rayleigh':
        rayleigh_amplitude = rng.rayleigh(scale=np.sqrt(mean_snr_lin/2),size=(MC,n_dims))
    else:
        raise ValueError('Only rayleigh distribution is implemented atm.')
        
#%%% loop over MC runs
    for mc_idx in range(MC):
        #########################################################################
        #                            Settings                                   #
        #########################################################################
        

        
        #########################################################################
        #                            TX                                         #
        #########################################################################
        
        # construct signal
        sig_tx = comm.signal.Signal(n_dims=n_dims)
        sig_tx.symbol_rate = 12.8e9
        
        #fix TX upsample for now to 2, to avoid signal degradation in resample with interpolation
        #TX_UPSAMPLE_FACTOR = DAC_SR / sig_tx.symbol_rate[0]
        TX_UPSAMPLE_FACTOR = 2
        
        # generate bits
        # fix for 4QAM right now
        wanted_symbols = amount_of_symbols
        n_bits = int((np.log2(mod_order)*wanted_symbols)-((np.log2(mod_order)*wanted_symbols)%np.log2(mod_order)))
        sig_tx.generate_bits(n_bits=n_bits, seed=1)
        
        # set constellation (modulation format)
        sig_tx.generate_constellation(format=mod_format, order=mod_order)
        
        # create symbols
        sig_tx.mapper()
        
        # upsampling and pulseshaping
        sig_tx.pulseshaper(upsampling=TX_UPSAMPLE_FACTOR, pulseshape='rrc', roll_off=ROLL_OFF)
        
        #########################################################################
        #                            CH                                         #
        #########################################################################
        
        #grab tx object
        sig_ch = copy.deepcopy(sig_tx)
        
        sig_ch, return_dict, _ = gen_SIMO_samples(sig_ch, max_phase_offset_in_rad=np.pi, max_timedelay_in_percent=10, n_apertures=n_dims, repeat=5, cut_to=3)
        
        # # add amplitude noise
        # for dim in range(sig_ch.n_dims):
        #    sig_ch.samples[dim] = comm.channel.set_snr(sig_ch.samples[dim], snr_dB=SNR[dim], sps=sig_ch.sample_rate[dim]/sig_ch.symbol_rate[dim], seed=None)
        
        # AWGN
        sps = np.array(sig_ch.sample_rate)/np.array(sig_ch.symbol_rate)
        # normalize each dimension to have a mean power of 1
        for j in range(sig_ch.n_dims):
            sig_ch.samples[j] = sig_ch.samples[j]/np.sqrt(np.mean(np.abs(sig_ch.samples[j])**2))
        # scale samples with rayleigh_amplitude, add AWGN with mean=0,var=1, and power distributed equally among real and imaginary part
        for j in range(sig_ch.n_dims):
            #sig.samples[j] = sig.samples[j] * rayleigh_amplitude[i,j] + np.sqrt(0.5)*(rng.standard_normal(size=(sig.samples[j].size,)) + 1j*rng.standard_normal(size=(sig.samples[j].size,)))
            n = (rng.standard_normal(size=(sig_ch.samples[j].size,)) + 1j*rng.standard_normal(size=(sig_ch.samples[j].size,)))
            sig_ch.samples[j] = sig_ch.samples[j] * rayleigh_amplitude[mc_idx,j]/np.sqrt(sps[j]) + np.sqrt(0.5)*n
        
        ## phase noise emulation
        # TODO: call function with unique phase noise seeds per dimension, if different
        # LOs are being used 
        for dim in range(sig_ch.n_dims):
           sig_ch.samples[dim] = comm.channel.add_phase_noise(sig_ch.samples[dim] ,sig_ch.sample_rate[dim] , LASER_LINEWIDTH, seed=1)['samples']
        
        #########################################################################
        #                            RX                                         #
        #########################################################################
        
        # contruct rx signal structure
        sig_rx = copy.deepcopy(sig_ch)
        
        # From here: "standard" coherent complex baseband signal processing ############
        # resample to 2 sps
        sig_rx = comm.utils.resample(sig_rx,target_sps=2)
        
        #### combining
        # =============================================================================
        # combining block goes here:
        
        # 0. SNR estimation & shifting to samples[0] = Highest SNR, samples[1] = second highest SNR ... also needed for MRC comb!
        # 0.1 estimate snr of diff signals with spectrum method 
        # TODO: Try to wweak spectrums method to be more precise
        sig_range = np.array([-sig_rx.symbol_rate[0]/2-(sig_rx.symbol_rate[0]/2*ROLL_OFF),sig_rx.symbol_rate[0]/2+(sig_rx.symbol_rate[0]/2*ROLL_OFF)])
        noise_range = np.array([-sig_rx.symbol_rate[0]/2-(sig_rx.symbol_rate[0]/2*ROLL_OFF)-1e9,-sig_rx.symbol_rate[0]/2-(sig_rx.symbol_rate[0]/2*ROLL_OFF),sig_rx.symbol_rate[0]/2+(sig_rx.symbol_rate[0]/2*ROLL_OFF),sig_rx.symbol_rate[0]/2+1e9+(sig_rx.symbol_rate[0]/2*ROLL_OFF)])
        
        snr_lin_list = np.zeros(n_dims)
        for i in range(0,sig_rx.n_dims):
            x_in = np.fft.fftshift(np.fft.fftfreq(len(sig_rx.samples[i]),1/sig_rx.sample_rate[0]))
            y_in = np.abs(np.fft.fftshift(np.fft.fft(sig_rx.samples[i])))**2 
            snr_lin_list[i] = 10**(comm.utils.estimate_snr_spectrum(x_in,y_in,sig_range=sig_range, noise_range= noise_range, 
                                                              order=1,noise_bw=sig_rx.symbol_rate[0],plotting=False)/10)
        
        # 0.2 shift apertures, first one should be the best in terms of SNR
        order_list = np.zeros(n_dims)
        snr_lin_list_temp = copy.deepcopy(snr_lin_list)
        for i in range(0,len(snr_lin_list)):
            order_list[snr_lin_list_temp.argmax()] = i
            snr_lin_list_temp[snr_lin_list_temp.argmax()] = -100
        
        sig_rx.samples = [sig_rx.samples[int(i)] for i in order_list]
        
        # 1. Timedelay compensation (subsample)
        for i in range(0,sig_rx.n_dims):
            results = comm.rx.sampling_phase_adjustment(sig_rx.samples[i], sample_rate=sig_rx.sample_rate[0], symbol_rate=sig_rx.symbol_rate[0], shift_dir='both')
            sig_rx.samples[i] = results['samples_out']
        
        # 1.2 Timedelay compensation (sample)
        sig_rx, lag_list = comm.rx.comb_timedelay_compensation(sig_rx, method="crop", xcorr="abs")
        
        # cut sample dimensions to same size, since there seems to be a discrepancy when
        # going beyond n_dims = 2
        min_len = np.min(np.array([sig_rx.samples[dim].size for dim in range(n_dims)]))
        for dim in range(n_dims):
            sig_rx.samples[dim] = sig_rx.samples[dim][:min_len] 
        
        # 2. Phase comp. and combining
        samples_rolling_sum = np.zeros(len(sig_rx.samples[0]), dtype=complex)
        for i in range(0,sig_rx.n_dims):
            #phase compensation
            _, samples_second_sig_phase_align, est_phase = comm.rx.comb_phase_compensation(samples_rolling_sum, sig_rx.samples[i])
            sig_rx.samples[i] = samples_second_sig_phase_align
            if comb_method == "EGC":
                #EGC combining
                samples_rolling_sum += sig_rx.samples[i]
            elif comb_method == "MRC":
                samples_rolling_sum += sig_rx.samples[i]*snr_lin_list[i]
            else: 
                raise Exception("No/wrong comb_method specified.")
                
        # wrinting back combined signal
        sig_rx.samples[0] = samples_rolling_sum
        
        # =============================================================================
        
        # normalize samples to mean magnitude of original constellation
        mag_const = np.mean(abs(sig_rx.constellation[0]))
        mag_samples = np.mean(abs(sig_rx.samples[0]))
        sig_rx.samples[0] = sig_rx.samples[0] * mag_const / mag_samples
        
        
        # either blind adaptive filter....
        if adaptive_filter == True:    
            results = comm.rx.blind_adaptive_equalizer(sig_rx, n_taps=51, mu_cma=1e-4, 
                                                       mu_rde=1e-5, mu_dde=0.5, decimate=True, 
                                                       return_info=True, stop_adapting=-1, 
                                                       start_rde=5000*0, start_dde=5000*0)
            
            sig_rx = results['sig']
            h = results['h'][0]
            eps = results['eps'][0]
            if plotting:
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
        
        # ... or matched filtering
        else:
            # Rx matched filter
            sig_rx.raised_cosine_filter(roll_off=ROLL_OFF,root_raised=True) 
            if plotting:
                sig_rx.plot_constellation(0, hist=True, tit='constellation after matched filter')
        
            # crop samples here, if necessary
            sps = int(sig_rx.sample_rate[0] / sig_rx.symbol_rate[0])
            crop = 10*sps
            if crop != 0:
                sig_rx.samples = sig_rx.samples[0][crop:-crop]
            else:
                sig_rx.samples = sig_rx.samples[0]
        
        # To one sample per symbol
        START_SAMPLE = 0
        sps = sig_rx.sample_rate[0] / sig_rx.symbol_rate[0] # CHECK FOR INTEGER SPS!!!
        sig_rx.samples = sig_rx.samples[0][START_SAMPLE::int(sps)]
        if plotting:
            sig_rx.plot_constellation(0, hist=True, tit='constellation after EQ')
        
        # CPE
        viterbi = True
        # ...either VV
        if viterbi:
            cpe_results = comm.rx.carrier_phase_estimation_VV(sig_rx.samples[0], n_taps=101, 
                                                              filter_shape='wiener', mth_power=4, 
                                                              rho=.00001)
            sig_rx.samples = cpe_results['rec_symbols']
            est_phase = cpe_results['phi_est'].real
        # ...or BPS
        else:
            cpe_results = comm.rx.carrier_phase_estimation_bps(sig_rx.samples[0], sig_rx.constellation[0], 
                                                       n_taps=15, n_test_phases=45, const_symmetry=np.pi/2)
            sig_rx.samples = cpe_results['samples_corrected']
            est_phase = cpe_results['est_phase_noise']
        
        if plotting:
            plt.plot(est_phase)
            plt.title('estimated phase noise')
            plt.xlabel('time / symbols')
            plt.ylabel('phase / rad')
            plt.grid()
            plt.show()
        
            sig_rx.plot_constellation(hist=True, tit='constellation after CPE')
        
        # delay and phase ambiguity estimation and compensation
        sig_rx = comm.rx.symbol_sequence_sync(sig_rx, dimension=-1)
            
        # calc EVM
        evm = comm.utils.calc_evm(sig_rx, norm='max')
        if plotting:
            print("EVM: {:2.2%}".format(evm[0]))
        
        # estimate SNR
        # snr = comm.utils.estimate_SNR_evm(sig_rx, norm='rms', method='data_aided', opt=False)
        snr = np.zeros((n_dims,))
        for dim in range(n_dims):
            snr[dim] = 10*np.log10(comm.utils.estimate_SNR_m2m4(sig_rx.samples[dim], sig_rx.constellation[dim]))
        if plotting:
            print("set SNR: {:.2f} dB @ {} apertures, est. SNR: {:.2f} dB, comb with {}".format(SNR[0], n_dims, snr[0], comb_method))
        
        # decision and demapper
        sig_rx.decision()
        sig_rx.demapper()
        
        # BER counting
        ber_res = comm.rx.count_errors(sig_rx.bits[0], sig_rx.samples[0])
        if plotting:
            print('BER = {}'.format(ber_res['ber']))
            
#%%% write to output arrays
        
        BER_tmp[mc_idx] = ber_res['ber']
        SNR_comb_tmp[mc_idx] = snr[0]
        set_SNR_vec[mc_idx] = np.mean(rayleigh_amplitude[mc_idx,:]**2)
        
    BER_vec[SNR_idx] = np.mean(BER_tmp)
    SNR_comb_vec[SNR_idx] = 10*np.log10(np.mean(10**(SNR_comb_tmp/10)))
    mean_set_SNR_vec[SNR_idx] = 10*np.log10(np.mean(set_SNR_vec))

#%% plotting

plt.figure(1)
plt.plot(SNR_vec,mean_set_SNR_vec)
# EbN0 = 10*np.log10((10**(SNR_vec/10))/np.log2(mod_order))
    
# plt.figure(1)
# plt.semilogy(SNR_vec,BER_vec,color='r')
# plt.semilogy(10*np.log10((10**(SNR_vec/10))*np.log2(mod_order)),0.5*erfc(np.sqrt(10**(SNR_vec/10))),color='salmon')
# plt.legend(('Eb/N0 sim w/ combining','Eb/N0 theory w/ 1 aperture'))
# plt.xlabel('Eb/N0 [dB]')
# plt.xticks(SNR_vec[::2])
# plt.ylabel('BER')
# plt.title('BER over Eb/N0 with MRC @ {} apertures'.format(n_dims))
# plt.grid()
# plt.ylim([1e-12,1])
# plt.show()

# plt.figure(2)
# plt.plot(SNR_vec,SNR_comb_vec,color='r')
# plt.plot(SNR_vec,SNR_vec,color='salmon')
# plt.legend(('Eb/N0 per aperture vs. Eb/N0 post-combining','Eb/N0 single aperture'))
# plt.xlabel('Eb/N0 per aperture [dB]')
# plt.xticks(SNR_vec[::2])
# plt.ylabel('Eb/N0 post-combining [dB]')
# plt.title('Eb/N0 per aperture vs. Eb/N0 with MRC @ {} apertures'.format(n_dims))
# plt.grid()
# plt.show()