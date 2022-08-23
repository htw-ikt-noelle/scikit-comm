import comm 
import numpy as np
import matplotlib.pyplot as plt
import copy

#### diversity gain test
# number of apertures/antennas
n_dims = np.arange(2,21)
# number of simulation runs
MC_runs = 50

# combining function
def combining(sig_div,comb_method='MRC',est_method='spectrum'):
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

    Returns
    -------
    sig_combined : 1d array
        MRC combined signal.

    """
    # deep copy of signal
    sig = copy.deepcopy(sig_div)
    # error checks
    if len(sig.samples) < 2:
        raise TypeError('Signal must have at least two sample arrays in list for diversity combining to be performed.')
    # if len(sig.samples.shape) > 2:
    #     raise TypeError('Samples attribute of signal object may not have more than 2 dimensions.')
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
        # TODO: Polynomial fit fails for some reason, carefully check dimensions of y!
        snr_vec[i] = comm.utils.estimate_snr_spectrum(x,y,sig_range=np.array([-sig.symbol_rate[i]/2,sig.symbol_rate[i]/2]),
                                                      noise_range=np.array([-sig.symbol_rate[i]/2-3e9,-sig.symbol_rate[i]/2,sig.symbol_rate[i]/2,sig.symbol_rate[i]/2+3e9]),
                                                      order=1,noise_bw=sig.symbol_rate[i],plotting=False)
        # print estimated SNR
        # print('True vs. estimated SNR for channel {}: {:.2f} vs. {:.2f} dB.'.format(i,snr[i],snr_vec[i]))
        
    # scaling    
    if comb_method == 'MRC':
        for i in range(len(sig.samples)):
            sig.samples[i] = sig.samples[i] * (10**(snr_vec[i]/10))
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
mean_MRC_SDC_gain = np.full_like(n_dims,0,dtype='float')
mean_EGC_SDC_gain = np.full_like(n_dims,0,dtype='float')
mean_MRC_AVG_gain = np.full_like(n_dims,0,dtype='float')
mean_EGC_AVG_gain = np.full_like(n_dims,0,dtype='float')
# n_apertures loop
for i in n_dims:
    # MC loop
    MRC_SDC_gain = np.zeros((MC_runs,),dtype='float')
    EGC_SDC_gain = np.zeros((MC_runs,),dtype='float')
    MRC_AVG_gain = np.zeros((MC_runs,),dtype='float')
    EGC_AVG_gain = np.zeros((MC_runs,),dtype='float')
    for j in range(MC_runs):
        sig = comm.signal.Signal(n_dims=int(i))
        sig.symbol_rate = 5e9
        sig.sample_rate = 15e9
        snr = np.random.randint(0,20,size=(i,)).tolist()
        snr_seeds = np.random.randint(0,1000,size=(i,)).tolist()
        bit_seeds = np.tile(np.random.randint(0,1000,size=(1,)),i).tolist()
        n_bits = 2**14
        df = sig.sample_rate[0]/n_bits
        #### TX
        sig.generate_bits(n_bits=n_bits,seed=bit_seeds)
        sig.generate_constellation(format='QAM',order=4)
        sig.mapper()
        sig.pulseshaper(upsampling=sig.sample_rate[0]/sig.symbol_rate[0],pulseshape='rrc',roll_off=.0)
        #### CH
        sig.set_snr(snr_dB=snr,seed=snr_seeds)
        #### RX
        # SDC combining
        sig_comb_SDC = combining(sig,comb_method='SDC')
        x_comb_SDC = np.linspace(-(sig_comb_SDC.sample_rate[0]+df)/2,(sig_comb_SDC.sample_rate[0]-df)/2,sig_comb_SDC.samples[0].size)
        y_comb_SDC = np.abs(np.fft.fftshift(np.fft.fft(sig_comb_SDC.samples[0])))**2
        snr_comb_SDC = comm.utils.estimate_snr_spectrum(x_comb_SDC,y_comb_SDC,sig_range=np.array([-sig_comb_SDC.symbol_rate[0]/2,sig_comb_SDC.symbol_rate[0]/2]),
                                                      noise_range=np.array([-sig_comb_SDC.symbol_rate[0]/2-1e9,-sig_comb_SDC.symbol_rate[0]/2,sig_comb_SDC.symbol_rate[0]/2,sig_comb_SDC.symbol_rate[0]/2+1e9]),
                                                      order=1,noise_bw=sig_comb_SDC.symbol_rate[0],plotting=False)
        
        # MRC combining
        sig_comb_MRC = combining(sig,comb_method='MRC')
        x_comb_MRC = np.linspace(-(sig_comb_MRC.sample_rate[0]+df)/2,(sig_comb_MRC.sample_rate[0]-df)/2,sig_comb_MRC.samples[0].size)
        y_comb_MRC = np.abs(np.fft.fftshift(np.fft.fft(sig_comb_MRC.samples[0])))**2
        snr_comb_MRC = comm.utils.estimate_snr_spectrum(x_comb_MRC,y_comb_MRC,sig_range=np.array([-sig_comb_MRC.symbol_rate[0]/2,sig_comb_MRC.symbol_rate[0]/2]),
                                                      noise_range=np.array([-sig_comb_MRC.symbol_rate[0]/2-1e9,-sig_comb_MRC.symbol_rate[0]/2,sig_comb_MRC.symbol_rate[0]/2,sig_comb_MRC.symbol_rate[0]/2+1e9]),
                                                      order=1,noise_bw=sig_comb_MRC.symbol_rate[0],plotting=False)
        
        # EGC combining
        sig_comb_EGC = combining(sig,comb_method='EGC')
        x_comb_EGC = np.linspace(-(sig_comb_EGC.sample_rate[0]+df)/2,(sig_comb_EGC.sample_rate[0]-df)/2,sig_comb_EGC.samples[0].size)
        y_comb_EGC = np.abs(np.fft.fftshift(np.fft.fft(sig_comb_EGC.samples[0])))**2
        snr_comb_EGC = comm.utils.estimate_snr_spectrum(x_comb_EGC,y_comb_EGC,sig_range=np.array([-sig_comb_EGC.symbol_rate[0]/2,sig_comb_EGC.symbol_rate[0]/2]),
                                                      noise_range=np.array([-sig_comb_EGC.symbol_rate[0]/2-1e9,-sig_comb_EGC.symbol_rate[0]/2,sig_comb_EGC.symbol_rate[0]/2,sig_comb_EGC.symbol_rate[0]/2+1e9]),
                                                      order=1,noise_bw=sig_comb_EGC.symbol_rate[0],plotting=False)
        
        # SNR gain over SDC
        MRC_SDC_gain[j] = snr_comb_MRC - snr_comb_SDC
        EGC_SDC_gain[j] = snr_comb_EGC - snr_comb_SDC
        # SNR gain over single antenna with average SNR
        MRC_AVG_gain[j] = snr_comb_MRC - np.mean(np.array(snr))
        EGC_AVG_gain[j] = snr_comb_EGC - np.mean(np.array(snr))
        
    mean_MRC_SDC_gain[i-2] = np.mean(MRC_SDC_gain)
    mean_EGC_SDC_gain[i-2] = np.mean(EGC_SDC_gain)
    
    mean_MRC_AVG_gain[i-2] = np.mean(MRC_AVG_gain)
    mean_EGC_AVG_gain[i-2] = np.mean(EGC_AVG_gain)
# print('Average SNR over all {} channels: {:.2f} dB.'.format(n_dims,np.mean(np.array(snr))))
# print('Estimated SNR after EGC combining: {:.2f} dB.'.format(snr_comb_EGC))
# print('Estimated SNR after MRC combining: {:.2f} dB.'.format(snr_comb_MRC))
# print('MRC SNR gain over EGC = {:.2f} dB.'.format(snr_comb_MRC-snr_comb_EGC))
# print('MRC gain over Selection Combining = {:.2f} dB.'.format(snr_comb_MRC-np.max(np.array(snr))))

plt.figure(1)
plt.plot(n_dims,mean_MRC_SDC_gain,color='r')
plt.plot(n_dims,mean_EGC_SDC_gain,color='b')
plt.grid()
plt.xticks(ticks=n_dims)
plt.title('mean MRC/EGC SNR gain over SDC over {} runs'.format(MC_runs))
plt.xlabel('Number of antennas')
plt.ylabel('SNR gain [dB]')
plt.legend(("MRC gain over SDC",'EGC gain over SDC'))
plt.show()

plt.figure(2)
plt.plot(n_dims,mean_MRC_AVG_gain,color='r')
plt.plot(n_dims,mean_EGC_AVG_gain,color='b')
plt.grid()
plt.xticks(ticks=n_dims)
plt.title('mean MRC/EGC SNR gain over single antenna over {} runs'.format(MC_runs))
plt.xlabel('Number of antennas')
plt.ylabel('SNR gain [dB]')
plt.legend(("MRC gain over AVG",'EGC gain over AVG'))
plt.show()