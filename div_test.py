import comm 
import numpy as np
import matplotlib.pyplot as plt
import copy

#### diversity gain test
# number of apertures/antennas
n_dims = np.arange(2,11)
# number of simulation runs
MC_runs = 10
# Rayleigh distribution parameters
rayleigh_mean_dB = 20 # in dB
# rayleigh_scale = 2 # must be positive
# signal parameters
roll_off = 0.01
symb_rate = 5e9
samp_rate = 15e9
n_bits = 2**6
mod_format = 'QAM'
mod_order = 4
n_samples = int(n_bits/np.log2(mod_order)*(samp_rate/symb_rate)) # number of samples per MC run

# combining function
def combining(sig_div,comb_method='MRC',est_method='spectrum',snr_true=None):
    """
    Performs Diversity Combining of the rows of an incoming signal matrix,
    where each row represents the signal captured by an antenna of a SIMO 
    system. The SNR is calculated per row and the individual rows are then
    scaled according to the chosen combination method and added together.
    
    Different SNR estimation methods can be selected by altering the est_method
    parameter. 

    Parameters
    ----------
    sig_div : signal object
        n-dimensional signal object with list of sample arrays in the 'samples'
        attribute.
    comb_method : str, optional
        Combining method. MRC and EGC are available. The default is 'MRC'.
    est_method : str, optional
        SNR estimation method. The default is 'spectrum'.
    snr_true : ndarray
        If the true SNR values (in dB) are known and should be used instead of 
        the estimates, they can be passed as an array here.

    Returns
    -------
    sig_combined : signal object
        SIgnal after combining. The sample attribute now has the combined sample
        array in every dimension.

    """
    # deep copy of signal
    sig = copy.deepcopy(sig_div)
    # error checks
    if len(sig.samples) < 2:
        raise TypeError('Signal must have at least two sample arrays in list for diversity combining to be performed.')
    if est_method != 'spectrum':
        raise ValueError('No other SNR estimation methods besides spectrum are implemented yet.')
    
    # init vector of SNR estimates
    snr_vec = np.zeros((len(sig.samples),),dtype='float')
    # SNR estimation
    for i in range(len(sig.samples)):
        df = sig.sample_rate[i]/sig.samples[i].size
        x = np.linspace(-(sig.sample_rate[i]+df)/2,(sig.sample_rate[i]-df)/2,sig.samples[i].size,)
        y = np.abs(np.fft.fftshift(np.fft.fft(sig.samples[i])))**2
        # TODO: adjust signal range according to roll off factor
        snr_vec[i] = comm.utils.estimate_snr_spectrum(x,y,sig_range=np.array([-sig.symbol_rate[i]/2,sig.symbol_rate[i]/2]),
                                                      noise_range=np.array([-sig.symbol_rate[i]/2-3e9,-sig.symbol_rate[i]/2,sig.symbol_rate[i]/2,sig.symbol_rate[i]/2+3e9]),
                                                      order=1,noise_bw=sig.symbol_rate[i],plotting=False)
        # print estimated SNR
        # print('True vs. estimated SNR for channel {}: {:.2f} vs. {:.2f} dB.'.format(i,snr[i],snr_vec[i]))
        
    # replace estimated SNRs with true SNRs, eliminating possible estimation
    # error, if desired
    if snr_true:
        snr_vec = np.asarray(snr_true)
    # scaling           
    if comb_method == 'MRC':
        for i in range(len(sig.samples)):
            sig.samples[i] = sig.samples[i] * (10**(snr_vec[i]/20)) 
    elif comb_method == 'EGC':
        pass
    elif comb_method == 'SDC':
        mask = np.where(snr_vec == np.max(snr_vec),1,0)
        for i in range(len(sig.samples)):
            sig.samples[i] = sig.samples[i] * mask[i]
    else:
        raise ValueError("Combining method not implemented.")
    # combination
    # TODO: samples attribute is currently being overwritten - is this practical
    # for our simulation? Should we add a new attribute?
    sig.samples = np.sum(sig.samples,axis=0)
    return sig

#### MC simulation
mean_MRC_AVG_gain = np.full_like(n_dims,0,dtype='float')
mean_EGC_AVG_gain = np.full_like(n_dims,0,dtype='float')
# arrays for dumping combined sample arrays from each MC run into
MRC_MC_array = np.zeros((int(MC_runs*n_samples),),dtype='complex')
EGC_MC_array = np.zeros((int(MC_runs*n_samples),),dtype='complex')
# arrays for mean SNRs after combining
MRC_mean_snr = np.zeros((len(n_dims),),dtype='float')
EGC_mean_snr = np.zeros((len(n_dims),),dtype='float')
# mean true SNR across all MC runs
mean_true_snr = np.zeros((len(n_dims),),dtype='float')

# n_apertures loop
for i in n_dims:
    # MC loop
    MRC_AVG_gain = np.zeros((MC_runs,),dtype='float')
    EGC_AVG_gain = np.zeros((MC_runs,),dtype='float')
    mean_true_snr_MC = np.zeros((MC_runs,),dtype='float') # mean SNR across all antennas per run

    for j in range(MC_runs):
        # init signal
        sig = comm.signal.Signal(n_dims=int(i))
        sig.symbol_rate = symb_rate
        sig.sample_rate = samp_rate
        bit_seed = 1 
        df = sig.sample_rate[0]/n_samples
        # init RNG and rayleigh distribution
        rng = np.random.default_rng(seed=None) 
        # linearize rayleigh mean
        rayleigh_mean_lin = 10**(rayleigh_mean_dB/10)
        # calc rayleigh distribution scale from desired mean value
        rayleigh_scale = rayleigh_mean_lin/np.sqrt(np.pi/2)
        # set linear SNRs with rayleigh distribution
        snr_lin = rng.rayleigh(scale=rayleigh_scale,size=(i,))
        # convert linear SNRs to log
        snr_dB = (10*np.log10(snr_lin)).tolist()
        # init seeds for SNRs and bit vector
        # snr_seeds = np.random.randint(0,10000000,size=(i,)).tolist()
        
        #### TX
        sig.generate_bits(n_bits=n_bits,seed=bit_seed)
        sig.generate_constellation(format=mod_format,order=mod_order)
        sig.mapper()
        sig.pulseshaper(upsampling=int(sig.sample_rate[0]/sig.symbol_rate[0]),pulseshape='rrc',roll_off=roll_off)
        #### CH
        sig.set_snr(snr_dB=snr_dB,seed=None)
        #### RX
        # SDC combining
        sig_comb_SDC = combining(sig,comb_method='SDC',snr_true=snr_dB)
        x_comb_SDC = np.linspace(-(sig_comb_SDC.sample_rate[0]+df)/2,(sig_comb_SDC.sample_rate[0]-df)/2,sig_comb_SDC.samples[0].size)
        y_comb_SDC = np.abs(np.fft.fftshift(np.fft.fft(sig_comb_SDC.samples[0])))**2
        snr_comb_SDC = comm.utils.estimate_snr_spectrum(x_comb_SDC,y_comb_SDC,sig_range=np.array([-sig_comb_SDC.symbol_rate[0]/2,sig_comb_SDC.symbol_rate[0]/2]),
                                                      noise_range=np.array([-sig_comb_SDC.symbol_rate[0]/2-1e9,-sig_comb_SDC.symbol_rate[0]/2,sig_comb_SDC.symbol_rate[0]/2,sig_comb_SDC.symbol_rate[0]/2+1e9]),
                                                      order=1,noise_bw=sig_comb_SDC.symbol_rate[0],plotting=False)
        
        # MRC combining
        sig_comb_MRC = combining(sig,comb_method='MRC',snr_true=snr_dB)
        x_comb_MRC = np.linspace(-(sig_comb_MRC.sample_rate[0]+df)/2,(sig_comb_MRC.sample_rate[0]-df)/2,sig_comb_MRC.samples[0].size)
        y_comb_MRC = np.abs(np.fft.fftshift(np.fft.fft(sig_comb_MRC.samples[0])))**2
        snr_comb_MRC = comm.utils.estimate_snr_spectrum(x_comb_MRC,y_comb_MRC,sig_range=np.array([-sig_comb_MRC.symbol_rate[0]/2,sig_comb_MRC.symbol_rate[0]/2]),
                                                      noise_range=np.array([-sig_comb_MRC.symbol_rate[0]/2-1e9,-sig_comb_MRC.symbol_rate[0]/2,sig_comb_MRC.symbol_rate[0]/2,sig_comb_MRC.symbol_rate[0]/2+1e9]),
                                                      order=1,noise_bw=sig_comb_MRC.symbol_rate[0],plotting=False)
        
        # EGC combining
        sig_comb_EGC = combining(sig,comb_method='EGC',snr_true=snr_dB)
        x_comb_EGC = np.linspace(-(sig_comb_EGC.sample_rate[0]+df)/2,(sig_comb_EGC.sample_rate[0]-df)/2,sig_comb_EGC.samples[0].size)
        y_comb_EGC = np.abs(np.fft.fftshift(np.fft.fft(sig_comb_EGC.samples[0])))**2
        snr_comb_EGC = comm.utils.estimate_snr_spectrum(x_comb_EGC,y_comb_EGC,sig_range=np.array([-sig_comb_EGC.symbol_rate[0]/2,sig_comb_EGC.symbol_rate[0]/2]),
                                                      noise_range=np.array([-sig_comb_EGC.symbol_rate[0]/2-1e9,-sig_comb_EGC.symbol_rate[0]/2,sig_comb_EGC.symbol_rate[0]/2,sig_comb_EGC.symbol_rate[0]/2+1e9]),
                                                      order=1,noise_bw=sig_comb_EGC.symbol_rate[0],plotting=False)
        
        # matched filter
        sig_comb_SDC.samples[0] = comm.filters.raised_cosine_filter(sig_comb_SDC.samples[0],
                                                                    sig_comb_SDC.sample_rate[0],
                                                                    sig_comb_SDC.symbol_rate[0],
                                                                    roll_off,root_raised=True)
        sig_comb_MRC.samples[0] = comm.filters.raised_cosine_filter(sig_comb_MRC.samples[0],
                                                                    sig_comb_MRC.sample_rate[0],
                                                                    sig_comb_MRC.symbol_rate[0],
                                                                    roll_off,root_raised=True)
        sig_comb_EGC.samples[0] = comm.filters.raised_cosine_filter(sig_comb_EGC.samples[0],
                                                                    sig_comb_EGC.sample_rate[0],
                                                                    sig_comb_EGC.symbol_rate[0],
                                                                    roll_off,root_raised=True)        
        # # downsampling to 1sps
        # sig_comb_SDC.samples[0] = sig_comb_SDC.samples[0][::int(sig_comb_SDC.sample_rate[0]/sig_comb_SDC.symbol_rate[0])]
        # sig_comb_MRC.samples[0] = sig_comb_MRC.samples[0][::int(sig_comb_MRC.sample_rate[0]/sig_comb_MRC.symbol_rate[0])]
        # sig_comb_EGC.samples[0] = sig_comb_EGC.samples[0][::int(sig_comb_EGC.sample_rate[0]/sig_comb_EGC.symbol_rate[0])]
        # # symbol sequence sync
        # sig_comb_SDC = comm.rx.symbol_sequence_sync(sig_comb_SDC, dimension=-1)
        # sig_comb_MRC = comm.rx.symbol_sequence_sync(sig_comb_MRC, dimension=-1)
        # sig_comb_EGC = comm.rx.symbol_sequence_sync(sig_comb_EGC, dimension=-1)
        
        # write combined sample array to MC array
        MRC_MC_array[int(j*n_samples):int((j+1)*n_samples)] = sig_comb_MRC.samples[0]
        EGC_MC_array[int(j*n_samples):int((j+1)*n_samples)] = sig_comb_EGC.samples[0]
        # calc mean SNR across all antennas per run
        mean_true_snr_MC[j] = np.mean(snr_lin)
        
        # SNR gain over single antenna with average SNR
        # MRC_AVG_gain[j] = 10*np.log10((10**(snr_comb_MRC/10)) / (np.mean(10**(np.array(snr)/10))))
        # EGC_AVG_gain[j] = 10*np.log10((10**(snr_comb_EGC/10)) / (np.mean(10**(np.array(snr)/10))))
    
    # calc mean of true SNRs across all MC runs
    mean_true_snr[i-2] = 10*np.log10(np.mean(mean_true_snr_MC))
    # estimate mean SNR of MC array
    est_x_axis = np.linspace(-(samp_rate+(samp_rate/n_samples))/2,(samp_rate-(samp_rate/n_samples))/2,MRC_MC_array.size)
    MRC_y_axis = np.abs(np.fft.fftshift(np.fft.fft(MRC_MC_array)))**2
    EGC_y_axis = np.abs(np.fft.fftshift(np.fft.fft(EGC_MC_array)))**2
    
    MRC_mean_snr[i-2] = comm.utils.estimate_snr_spectrum(est_x_axis,MRC_y_axis,sig_range=np.array([-symb_rate/2,symb_rate/2]),
                                                        noise_range=np.array([-symb_rate/2-3e9,-symb_rate/2,symb_rate/2,symb_rate/2+3e9]),
                                                        order=1,noise_bw=symb_rate,scaling='lin',plotting=True)
    EGC_mean_snr[i-2] = comm.utils.estimate_snr_spectrum(est_x_axis,EGC_y_axis,sig_range=np.array([-symb_rate/2,symb_rate/2]),
                                                        noise_range=np.array([-symb_rate/2-3e9,-symb_rate/2,symb_rate/2,symb_rate/2+3e9]),
                                                        order=1,noise_bw=symb_rate,scaling='lin',plotting=False) 
    
    # downsample
    MRC_MC_array_tmp = MRC_MC_array[::int(samp_rate/symb_rate)]
    EGC_MC_array_tmp = EGC_MC_array[::int(samp_rate/symb_rate)]
    # MRC_mean_snr[i-2] = comm.utils.estimate_SNR_evm()
    
    # calc mean SNR gain
    mean_MRC_AVG_gain[i-2] = 10*np.log10((10**(MRC_mean_snr[i-2]/10)) / (10**(mean_true_snr[i-2]/10)))
    mean_EGC_AVG_gain[i-2] = 10*np.log10((10**(EGC_mean_snr[i-2]/10)) / (10**(mean_true_snr[i-2]/10)))

# close all figures
plt.close('all')

plt.figure(1)
plt.plot(n_dims,mean_MRC_AVG_gain,color='r')
plt.plot(n_dims,mean_EGC_AVG_gain,color='b')
plt.plot(n_dims,10*np.log10(n_dims),color='salmon')
plt.plot(n_dims,10*np.log10(n_dims*np.pi/4),color='steelblue')
plt.grid()
plt.xticks(ticks=n_dims)
plt.title('mean MRC/EGC SNR gain over single (average) antenna over {} runs'.format(MC_runs))
plt.xlabel('Number of antennas')
plt.ylabel('SNR gain [dB]')
plt.legend(("MRC gain over AVG",'EGC gain over AVG','MRC theory','EGC theory'))
plt.show()