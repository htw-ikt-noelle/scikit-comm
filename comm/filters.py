import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt


def filter_samples(samples, h, domain='freq'):
    """ Filter signal with fiven impulse response.
    
    Filter is either implemented in time domain (convolution) or in frequency
    domain (multiplication).
    
    """
    
    # check, if input is real    
    isreal = np.alltrue(np.isreal(samples))
    
    if domain == 'time':
        samples_out = signal.convolve(samples, h, mode='same')
    elif domain == 'freq':
        samples_out = np.fft.ifft(np.fft.fft(samples) * np.fft.fft(h, n=samples.size))
    else:
        raise ValueError('filter_samples: domain must either be "time" or "freq" ...')    
        
    if isreal:
        samples_out = np.real(samples_out)
    
    return samples_out
    


def raised_cosine_filter(samples, sample_rate=1.0, symbol_rate=1.0, roll_off=0.0, 
                         length=-1, root_raised=False, domain='freq'):
    """ Filter a given signal with a (root) raised cosine filter.
    
    time domain
    bla bla bla
    
    """
    
    # set parameters
    sps = sample_rate / symbol_rate
    if length == -1:
        N = samples.size
    else:
        N = length * sps # length of impulse response in number of symbols
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
    
    time domain
    bla bla bla
    
    """
    
    # generate impulse response
    h = np.ones(average)/average
    
    # actual filtering
    samples_out = filter_samples(samples, h, domain)    
    
    return samples_out


def windowed_sinc(samples, fc=0.5, order=111, window=None, domain='freq'):
    """ Filter a given signal windowed Si-funtion as impulse response.
    
    time domain
    0<fc<1, where 1 specifies the Nyquist frequency (half the sampling frequency)
    window can be 
    'none' --> Si impulse response
    'Hamming'
    'Blackmann-Harris'
    order has to be odd    
    
    """
    
    # calc paramters
    if (order % 2) == 0:
        raise ValueError('windowed_sinc: order has to be odd...')        
    
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
        
    # debug plots...
    # f = np.fft.fftshift(np.fft.fftfreq(order))
    # plt.plot(f, np.abs(np.fft.fftshift(np.fft.fft(h))))
    
    # actual filtering
    samples_out = filter_samples(samples, h, domain)
    
    return samples_out


def ideal_lp(samples, fc):
    """ Filter a given signal with an ideal lowpass filter.
    
    frequency domain
    0<fc<1, where 1 specifies the Nyquist frequency (half the sampling frequency)
    
    
    """
    
    # generate ideal frequency response
    f = np.fft.fftfreq(samples.size)
    H = np.zeros_like(samples)
    H[np.abs(f) <= fc/2] = 1
    
    # calc impulse response with ifft    
    h = np.real(np.fft.ifft(H))
    
    # actual filtering
    samples_out = filter_samples(samples, h, domain='freq') 
    
    return samples_out