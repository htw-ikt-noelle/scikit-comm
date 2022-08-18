import comm 
import numpy as np
import matplotlib.pyplot as plt

#### diversity gain test
sig = comm.signal.Signal(n_dims=2)
sig.symbol_rate = 5e9
sig.sample_rate = 15e9
x_pol_snr = np.random.randint(0,20)
y_pol_snr = np.random.randint(0,20)
n_bits = 2**16
df = sig.sample_rate[0]/n_bits
#### TX
sig.generate_bits(n_bits=n_bits,seed=None)
sig.generate_constellation(format='QAM',order=4)
sig.mapper()
sig.pulseshaper(upsampling=sig.sample_rate[0]/sig.symbol_rate[0],pulseshape='rrc',roll_off=.2)
#### CH
sig.set_snr(snr_dB=[x_pol_snr,y_pol_snr],seed=None)
# x_pol = comm.channel.set_snr(sig.samples[0],snr_dB=x_pol_snr,seed=1)
# y_pol = comm.channel.set_snr(sig.samples[1],snr_dB=y_pol_snr,seed=1)
#### RX
# pre-combining SNR estimation
x_snr_est = comm.utils.estimate_snr_spectrum(np.linspace(-sig.sample_rate[0]/2+df/2,sig.sample_rate[0]/2-df/2,sig.samples[0].size), 
                                             np.abs(np.fft.fftshift(np.fft.fft(sig.samples[0])))**2, 
                                             sig_range=np.array([-3e9,3e9]), 
                                             noise_range=np.array([-4e9,-3e9,3e9,4e9]),
                                             order=1,
                                             noise_bw=sig.sample_rate[0],
                                             plotting=False)

y_snr_est = comm.utils.estimate_snr_spectrum(np.linspace(-sig.sample_rate[0]/2+df/2,sig.sample_rate[0]/2-df/2,sig.samples[0].size), 
                                             np.abs(np.fft.fftshift(np.fft.fft(sig.samples[1])))**2, 
                                             sig_range=np.array([-3e9,3e9]), 
                                             noise_range=np.array([-4e9,-3e9,3e9,4e9]),
                                             order=1,
                                             noise_bw=sig.sample_rate[0],
                                             plotting=False)
# EGC combining
sig_EGC = sig.samples[0] + sig.samples[1]
EGC_snr_est = comm.utils.estimate_snr_spectrum(np.linspace(-sig.sample_rate[0]/2+df/2,sig.sample_rate[0]/2-df/2,sig.samples[0].size), 
                                             np.abs(np.fft.fftshift(np.fft.fft(sig_EGC)))**2, 
                                             sig_range=np.array([-3e9,3e9]), 
                                             noise_range=np.array([-4e9,-3e9,3e9,4e9]),
                                             order=1,
                                             noise_bw=sig.sample_rate[0],
                                             plotting=False)
# MRC combining
sig_MRC = (10**(x_snr_est/10))*sig.samples[0] + (10**(y_snr_est/10))*sig.samples[1]
MRC_snr_est = comm.utils.estimate_snr_spectrum(np.linspace(-sig.sample_rate[0]/2+df/2,sig.sample_rate[0]/2-df/2,sig.samples[0].size), 
                                             np.abs(np.fft.fftshift(np.fft.fft(sig_MRC)))**2, 
                                             sig_range=np.array([-3e9,3e9]), 
                                             noise_range=np.array([-4e9,-3e9,3e9,4e9]),
                                             order=1,
                                             noise_bw=sig.sample_rate[0],
                                             plotting=False)

print('Sum of individual channel SNRs = {:.2f} dB.'.format(x_snr_est+y_snr_est))
print('MRC SNR = {:.2f} dB.'.format(MRC_snr_est))
print('MRC SNR gain = {:.2f} dB.'.format(MRC_snr_est-EGC_snr_est))

def combining(sig,comb_method='MRC',est_method='spectrum'):
    """
    Performs Diversity Combining of the rows of an incoming signal matrix,
    where each row represents the signal captured by an antenna of a SIMO 
    system. The SNR is calculated per row and the individual rows are then
    scaled according to the chosen combination method and added together.
    
    Different SNR estimation methods can be selected by altering the est_method
    parameter. 

    Parameters
    ----------
    sig : signal object
        n-dimensional signal object with list of sample arrays in the 'samples'
        attribute.
    comb_method : str, optional
        Combination method. The default is 'MRC'.
    est_method : str, optional
        SNR estimation method. The default is 'spectrum'.

    Returns
    -------
    sig_combined : 1d array
        MRC combined signal.

    """
    # error checks
    if sig.samples.shape[0] < 2:
        raise TypeError('Signal must have at least two sample arrays in list for diversity combining to be performed.')
    if len(sig.samples.shape) > 2:
        raise TypeError('Samples attribute of signal object may not have more than 2 dimensions.')
    if est_method != 'MRC':
        raise ValueError('No other SNR estimation methods besides spectrum are implemented yet.')
    
    # vector of SNR estimates
    snr_vec = np.zeros((sig.samples.shape[0],),dtype='float')
    # SNR estimation
    for i in range(sig.samples.shape[0]):
        df = sig.sample_rate[i]/sig.samples[i].size
        x = np.linspace(-(sig.sample_rate[i]+df)/2,(sig.sample_rate[i]-df)/2,sig.samples[i].size)
        y = np.abs(np.fft.fftshift(np.fft.fft(sig.samples[i])))**2
        # TODO: adjust signal range according to roll off factor
        snr_vec[i] = comm.utils.estimate_snr_spectrum(x,y,sig_range=np.array([-sig.sample_rate[i]/2,sig.sample_rate[i]/2]),
                                                      noise_range=np.array([-sig.sample_rate[i]/2-1e9,-sig.sample_rate[i]/2,sig.sample_rate[i]/2,sig.sample_rate[i]/2+1e9]),
                                                      order=1,plotting=False)
        
    # scaling
    for i in range(sig.samples.shape[0]):
        sig.samples[i] = sig.samples[i] * (10**(snr_vec[i]/10))
    # combination
    # TODO: samples attribute is currently being overwritten - is this practical
    # for our simulation? Should we add a new attribute?
    sig.samples = np.sum(sig.samples,axis=0)
    return sig