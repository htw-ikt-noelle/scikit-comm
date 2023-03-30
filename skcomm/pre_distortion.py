import numpy as np
import scipy.signal as ssignal
import matplotlib.pyplot as plt

from . import filters


def generate_wn_probesignal(n_samples=2**17, f_max=1.0):
    """
    Generate complex white noise samples.
    
    The amplitudes of the real and imaginary 
    parts are each uniformly distributed between -1.0 and 1.0.
    
    

    Parameters
    ----------
    n_samples : integer, optional
        Number of noise samples to generate. The default is 2**17.
    f_max : float, optional
        cut off frequency, 0.0 < f_max <= 1.0, where 1.0 specifies the Nyquist 
        frequency (half the sampling frequency). The default is 1.0.
    

    Returns
    -------
    samples : 1D numpy array, complex
        Random complex noise samples where the real and imaginary parts are each
        unifomrly distributed bewtween -1.0 and 1.0.

    """
    # generate complex white noise samples
    samples = np.random.uniform(low=-1, high=1, size=n_samples*2).view(np.complex128)
    
    if f_max != 1.0:        
        # pre-filtering input noise
        samples = filters.ideal_lp(samples, fc=f_max)['samples_out']
        
    return samples


def estimate_tf_welch(samples_in, sample_rate_in, samples_out, sample_rate_out, f_max, nperseg=64, visualize=0):
    """
    Estimate magnitude transfer function (in linear scale) using Welch method 
    and calculate the inverted transfer function up to f_max.
    
    For calculation of the transfer function, the samples_out are resampled to the 
    input sample rate, if the samples rates do not match.
    
    See documentation of scipy.signal.welch for mor information.
    

    Parameters
    ----------
    samples_in : 1D numpy array, real or complex
        input samples to the DUT.
    sample_rate_in : float
        samples rate of samples_in in Hz.
    samples_out : 1D numpy array, real or complex
        output samples from DUT.
    sample_rate_out : float
        samples rate of samples_out in Hz..
    nperseg : int, optional
        block size of Welch method. The default is 64.    
    f_max : float
        frequency up to which the inversion of the transfer function is calculeted.
        The inverted magnitude transfer function is set to 0 outside this frequency
        range.
    visualize : int, optional
        should debug plots (spectra) be generted? Value also specifies the figure
        number. No plots for visualize=0. The default is 0.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    results : dict containing following keys
        tf : 1D numpy array, real or complex
            estimated linear magnitude transfer function.
        tf_inv : 1D numpy array, real or complex
            inverted linear magnitude transfer function up to frequency f_max.
        freq : 1D numpy array, real
            frequency axis of tf / tf_inv

    """    
    if sample_rate_in != sample_rate_out:
        #  do we need an AA filter?
        if sample_rate_out > sample_rate_in:
            f_c = sample_rate_in/sample_rate_out
            samples_out = filters.ideal_lp(samples_out, fc=f_c)['samples_out']            
            
        # resample output to input samplerate
        # check that this is really an integer, otherwise the samplerate is asynchronous with the data afterwards!!!
        len_dsp = sample_rate_in / sample_rate_out * np.size(samples_out)
        if len_dsp % 1:
            raise ValueError('DSP samplerate results in asynchronous sampling of the data symbols')
            # resampling to input samplerate
        samples_out = ssignal.resample(samples_out, num=int(len_dsp), window=None)   
        sample_rate_out = sample_rate_in
    
        
    
    spectrogram_in = ssignal.welch(samples_in, fs=sample_rate_in, window='hann', 
                                 nperseg=nperseg, noverlap=None, nfft=None, detrend=False, 
                                 return_onesided=False, scaling='spectrum', axis=- 1, average='mean')
    
    # shifting the zero-frequency component to the center of the spectrum
    freq = np.fft.fftshift(spectrogram_in[0]) 
    
    # magnitude spectrum of input signal
    mag_in = np.fft.fftshift(np.sqrt(spectrogram_in[1]))
    
    # prevent zeros in input magnitude spectrum
    mag_in[mag_in < 1e-100] = 1e-100
    
    
    if visualize != 0:
        # mag. spectrum in dB
        mag_in_dB = 20*np.log10(mag_in) 
        plt.figure(visualize)
        plt.plot(freq,mag_in_dB,'og-')
        plt.xlabel('Frequency (Hz)'); plt.ylabel('magnitude (dB)')           
        
        
    spectrogram_out = ssignal.welch(samples_out, fs=sample_rate_out, window='hann', 
                                  nperseg=nperseg, noverlap=None, nfft=None, detrend=False, 
                                  return_onesided=False, scaling='spectrum', axis=- 1, average='mean')

    # spec_out_re = ssignal.welch(np.real(samples_out), fs=sample_rate_out, window='hann', 
    #                               nperseg=nperseg, noverlap=None, nfft=None, detrend=False, 
    #                               return_onesided=False, scaling='spectrum', axis=- 1, average='mean')
    
    # spec_out_im = ssignal.welch(np.imag(samples_out), fs=sample_rate_out, window='hann', 
    #                               nperseg=nperseg, noverlap=None, nfft=None, detrend=False, 
    #                               return_onesided=False, scaling='spectrum', axis=- 1, average='mean')

    
    mag_out = np.fft.fftshift(np.sqrt(spectrogram_out[1])) # mag. spectrum
    
    # psd_out_re_dB = 20*np.log10(np.fft.fftshift(spec_out_re[1])) # to dB
    # psd_out_im_dB = 20*np.log10(np.fft.fftshift(spec_out_im[1])) # to dB
    
    if visualize != 0:
        mag_out_dB = 20*np.log10(mag_out) # to dB
        plt.figure(visualize);
        plt.plot(freq,mag_out_dB,'b.-')
        # plt.plot(freq,psd_out_re_dB,'b.-')
        # plt.plot(freq,psd_out_im_dB,'r.-')       
                
    # Magnitude transfer Function
    tf = mag_out / mag_in
    
    if visualize != 0:
        tf_dB = 20*np.log10(tf) # to dB
        plt.figure(visualize);
        plt.plot(freq,tf_dB,'k.-')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('magnitude (dB)')
        plt.legend(('input','output', 'transfer function'))                
        plt.grid(True)
        plt.show()
        
    
    # prevent zeros in tf
    tf[tf < 1e-100] = 1e-100
    # invert tf
    tf_inv = 1/tf
    # crop to usable min/max freq. range 
    tf_inv[np.abs(freq) > f_max] = 1e-100 
    
    # generate results dict
    results = dict()
    results['tf'] = tf
    results['tf_inv'] = tf_inv
    results['freq'] = freq
    
    return results


def dac_sinc_correction(samples, f_max=1.0):
    """
    Compensates for the Sinc rolloff of a zero order hold digital-to-analogue converter.
    
    Parameters
    ----------
    samples : 1D numpy array, real or complex
        input signal.
    f_max : float, optional
        frequency up to which the rolloff is compensated for. 0.0 < f_max <= 1.0, 
        where 1.0 specifies the Nyquist frequency (half the sampling frequency).
        The default is 1.0.
        

    Returns
    -------
    samples_out : 1D numpy array, real or complex
        output signal.

    """
    
    # frequency axis
    f = np.fft.fftshift(np.fft.fftfreq(samples.size))
    # transfer function
    Hf = 1 / np.sinc(f / 1)
    # all frequencies above f_max will not be touched
    if f_max:
        Hf[np.abs(f) > f_max/2] = 1
    # filter / pre-distort
    samples_out = filters.filter_samples(samples, Hf, domain='freq')
    
    return samples_out