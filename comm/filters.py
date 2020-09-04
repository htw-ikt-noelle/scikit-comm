import numpy as np
import scipy.signal as signal
import math as math
import matplotlib.pyplot as plt


def filter_samples(samples, h, domain='freq'):
    """ Filter signal with given impulse response.
    
    Filter is either implemented in either in the 
    
        -time domain (convolution), filter data along one-dimension with an FIR filter, see 
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lfilter.html
    
        -frequency domain (multiplication of spectra equivalent to a cyclic convolution), in this case the given impulse
        response h is zero-padded or truncated to fit the length of the input singal.
    
    """
    
    # check, if input is real    
    isreal = np.alltrue(np.isreal(samples))
    
    if domain == 'time':
        samples_out = signal.lfilter(h, np.ones(1), samples)
    elif domain == 'freq':
        samples_out = np.fft.ifft(np.fft.fft(samples) * np.fft.fft(h, n=samples.size))
    else:
        raise ValueError('filter_samples: domain must either be "time" or "freq" ...')    
        
    if isreal:
        samples_out = np.real(samples_out)
    
    return samples_out
    


def raised_cosine_filter(samples, sample_rate=1.0, symbol_rate=1.0, roll_off=0.0, 
                         length=-1, root_raised=False, domain='freq'):
    """
    Filter a given signal with a (root) raised cosine filter.
    
     -time domain (convolution), filter data along one-dimension with an FIR filter, see 
     https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lfilter.html
        
     -frequency domain (multiplication of spectra equivalent to a CYCLIC convolution), in this case the given impulse
     response h is zero-padded or truncated to fit the length of the input singal.
     
     CAUTION: This filter generates a group delay which is equal to ceil(size(h)/2).

    Parameters
    ----------
    samples : 1D numpy array, real or complex
        input signal.
    sample_rate : float, optional
        sample rate of input signal in Hz. The default is 1.0.
    symbol_rate : float, optional
        symbol rate of the input signal in Bd. The default is 1.0.
    roll_off : float, optional
        roll off factor of the filter. The default is 0.0.
    length : int, optional
        length of the filter impulse response, -1 equals to the length of the input singal. The default is -1.
    root_raised : bool, optional
        is the filter a root-raised cosine filter (True) or a raised-cosine filter (False). The default is False.
    domain : string, optional
        implementation of the filter either in 'time' or in 'freq' domain. The default is 'freq'.

    Returns
    -------
    samples_out : 1D numpy array, real or complex
        filtered output signal.

    """
    
    
    # set parameters
    sps = sample_rate / symbol_rate
    if length == -1:
        N = samples.size
    else:
        N = length # length of impulse response in number of samples
    t_filter = np.arange(-np.ceil(N/2)+1, np.floor(N/2)+1)
    T = sps
    
    # generate impulse response
    if root_raised:
        # root-raised cosine filter
        with np.errstate(divide='ignore',invalid='ignore'):# avoid raising a divide by zero / NaN warning
            h = (np.sin(np.pi * t_filter / T * (1-roll_off)) + 4 * roll_off * t_filter / T * np.cos(np.pi * t_filter / T * (1 + roll_off))) / (np.pi * t_filter / T * (1 - (4 * roll_off * t_filter / T)**2))
        h[t_filter==0] = (1 - roll_off + 4 * roll_off / np.pi)
        if roll_off != 0.0:
            h[np.abs(t_filter)==T/4/roll_off] = roll_off / np.sqrt(2) * ((1 + 2 / np.pi) * np.sin(np.pi / 4 / roll_off) + (1 - 2 / np.pi) * np.cos(np.pi / 4 / roll_off))
    else:
        # raised cosine filter
        with np.errstate(divide='ignore',invalid='ignore'):# avoid raising a divide by zero / NaN warning
            h = ((np.sin(np.pi*t_filter/T)) / (np.pi*t_filter/T)) * ((np.cos(roll_off*np.pi*t_filter/T)) / (1-(2*roll_off*t_filter/T)**2))
        h[t_filter==0] = 1
        if roll_off != 0.0:
            h[np.abs(t_filter) == (T/(2*roll_off))] = np.sin(np.pi/2/roll_off) / (np.pi/2/roll_off) * np.pi / 4
    
    # actual filtering
    samples_out = filter_samples(samples, h, domain)
    
    return samples_out




def moving_average(samples, average, domain='freq'):
    """ Filter a given signal with a moving average filter.
    
    In the time domain, a causal moving average filter is used, while the last 
    average-1 samples of the convolution are truncated.
    
    In the frequency domain a cyclic convoultion is used.
    
    """
    
    # generate impulse response
    h = np.ones(average)/average
    
    # actual filtering
    samples_out = filter_samples(samples, h, domain)    
    
    return samples_out


def windowed_sinc(samples, fc=0.5, order=111, window=None, domain='freq'):
    """
    Filter a given signal windowed Si-funtion as impulse response.

    Parameters
    ----------
    samples : 1D numpy array, real or complex
        input signal..
    fc : float, optional
        cut off frequency, 0.0 <= fc <= 1.0, where 1.0 specifies the Nyquist frequency (half the sampling frequency). The default is 0.5.
    order : int, optional
        length of the inpulse response, which has to be odd. If order=-1, the length is chosen to be the length of the input signal. The default is 111.
    window : string, optional
        type of the window funtion which is multiplied with the sinc-impulse response before filtering.
        'none' --> Si impulse response
        'Hamming'
        'Blackmann-Harris'. The default is None.
    domain : string, optional
        implementation of the filter either in 'time' or in 'freq' domain. The default is 'freq'.
        
    CAUTION: This filter generates a group delay which is equal to ceil(order/2).

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    samples_out : 1D numpy array, real or complex
        filtered output signal.

    """    
    
    # calc paramters
    if (order % 2) == 0:
        raise ValueError('windowed_sinc: order has to be odd...')        
    
    if order == -1:
        order = np.size(samples)
        
    n = np.arange(-np.floor(order/2), np.ceil(order/2), 1)
    
    # generate Si impulse response
    with np.errstate(divide='ignore',invalid='ignore'):# avoid raising a divide by zero / NaN warning
        h = np.sin(n * np.pi * fc) / (n * np.pi)
    h[n==0] = fc
    h /= np.max(h)
    
    # window in time domain
    if window == None:
        pass
    elif window == 'Hamming':
        h *= signal.windows.hamming(order)
    elif window == 'Blackman-Harris':
        h *= signal.windows.blackmanharris(order)
        
    # # debug plots...
    # plt.figure(1)
    # f = np.fft.fftshift(np.fft.fftfreq(order))
    # plt.plot(f, np.abs(np.fft.fftshift(np.fft.fft(h))))
    # plt.show()
    
    # actual filtering
    samples_out = filter_samples(samples, h, domain)
    
    return samples_out


def ideal_lp(samples, fc):
    """
    Filter a given signal with an ideal lowpass filter.
    
    This filter is only implemented in the frequency domain (cyclic convolution)
    and has NO group delay.
    
    Parameters
    ----------
    samples : 1D numpy array, real or complex
        input signal.
    fc : float, optional
        cut off frequency, 0.0 <= fc <= 1.0, where 1.0 specifies the Nyquist frequency (half the sampling frequency). The default is 0.5.
    
    Returns
    -------
    results : dict containing following keys
        samples_out : 1D numpy array, real or complex
            filtered output signal.
        real_fc : float
            actual applied cut off frequency (matching the frequency grid of the FFT)
    """ 
    # generate ideal frequency response
    f = np.fft.fftfreq(samples.size, d=0.5)
    H = np.zeros_like(samples)
    H[np.abs(f) <= fc] = 1
    
    real_fc = np.max(np.abs(f[np.abs(f) <= fc]))
    
    # calc impulse response with ifft    
    h = np.real(np.fft.ifft(H))
    
    # actual filtering
    samples_out = filter_samples(samples, h, domain='freq')
    
    # generate results dict
    results = dict()
    results['samples_out'] = samples_out
    results['real_fc'] = real_fc
    
    return results

def time_shift(samples, sample_rate=1.0, tau=0.0):
    """
    Add a cyclic time shift to the input signal.
    
    A positve time shift tau delays (right shifts) the signal, while a negative 
    time shift advances (left shifts) it. For time shifts equal to an integer sampling duration, 
    the signal is simply rolled.

    Parameters
    ----------
    samples :  1D numpy array, real or complex
        input signal.
    sample_rate : float, optional
        sample rate of input signal in Hz. The default is 1.0.
    tau : float, optional
        time shift in s. The default is 0.0.

    Returns
    -------
    samples_out : 1D numpy array, real or complex
        cyclic time shifted input signal.
    """    
       
    # integer sample shift is simple roll operation
    if tau%(1/sample_rate) == 0.0:
        shift = int(tau*sample_rate)
        samples_out = np.roll(samples, shift)    
    # fractional sample shift
    else:    
        # check, if input is real    
        isreal = np.alltrue(np.isreal(samples))
        # frequency vector
        w = np.fft.fftfreq(np.size(samples, axis=0), d=1/sample_rate) * 2 * np.pi    
        
        samples_out = np.fft.ifft(np.fft.fft(samples) * np.exp(-1j * w * tau))
        
        if isreal:
            samples_out = np.real(samples_out)
            
    # # for debugging purpose
    # plt.plot(np.abs(samples[:100]), 'C0')
    # plt.plot(np.abs(samples_out[:100]), 'C1')
    
    return samples_out
    
    