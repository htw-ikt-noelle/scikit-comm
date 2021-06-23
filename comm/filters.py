import numpy as np
import scipy.signal as signal
import math as math
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def filter_samples(samples, filter, domain='freq'):
    """ Filter the input signal.
    
    
    Filter is either implemented in either in the 
    
        -time domain filter (convolution). CAUTION: size of signal vector changes!
    
        -frequency domain (multiplication of input signal spectrum with transfer 
                           function...equivalent to a cyclic convolution).
        
        
    Parameters
    ----------
    samples : 1D numpy array, real or complex
        input signal.
    filter : 1D numpy array, real or complex
        Either inpulse response of the filter (when domain='time') or transfer 
        function (double-sided, starting from negative frequencies) of the 
        filter (when domain='freq').    
    domain : string, optional
        implementation of the filter either in 'time' or in 'freq' domain. The default is 'freq'.

    Returns
    -------
    samples_out : 1D numpy array, real or complex
        filtered input signal.
    
    """
    
    if domain == 'time':
        samples_out = signal.convolve(samples, filter)
    elif domain == 'freq':
        if samples.shape != filter.shape:
            raise TypeError('shape of samples and filter must be equal')
        samples_out = np.fft.ifft(np.fft.ifftshift(np.fft.ifftshift(np.fft.fft(samples)) * filter))
    else:
        raise ValueError('filter_samples: domain must either be "time" or "freq" ...')  
    
    return samples_out


    


def raised_cosine_filter(samples, sample_rate=1.0, symbol_rate=1.0, roll_off=0.0, 
                         length=-1, root_raised=False, domain='freq'):
    """
    Filter a given signal with a (root) raised cosine filter.
    
     -time domain (convolution): data is convolved with an (root) raised cosine 
     impulse response (h). CAUTION: This filter generates a group delay which is 
     equal to ceil(size(h)/2). Further, the number of samples of the output signal
     differs from the number of samples of the input signal.
        
     -frequency domain (multiplication of spectra equivalent to a CYCLIC 
    convolution), the frequency response of an (root) raised cosine filter is
    multiplied with the spectrum of the signal. This filter is acausal and does
    not produce any group delay.
     
     

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
        length of the filter impulse response, -1 equals to the length of the 
        input singal. The default is -1. This parameter is only reasonably in 
        case of time domain filtering (domain='time').
    root_raised : bool, optional
        is the filter a root-raised cosine filter (True) or a raised-cosine 
        filter (False). The default is False.
    domain : string, optional
        implementation of the filter either in 'time' or in 'freq' domain. 
        The default is 'freq'.

    Returns
    -------
    samples_out : 1D numpy array, real or complex
        filtered output signal.

    """
    
    if samples.ndim != 1:
        raise ValueError('signal vector has to be a 1D array...')
    
    if (roll_off < 0.0) or (roll_off > 1.0):
        raise ValueError('roll_off needs to be between 0.0 and 1.0')
    
    # set parameters
    sps = sample_rate / symbol_rate
    
    # time domain implementation
    if domain == 'time':
    
        if length > samples.size:
            raise ValueError('impulse response should be shorter or equal than signal')
                
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
    
    # frequency domain implementation
    elif domain == 'freq':
        
        f = np.fft.ifftshift(np.fft.fftfreq(samples.size, 1/sample_rate))
        H = np.zeros(f.size)
        T = 1/symbol_rate
        
        # set constant part to 1 (see transfer function in Wikipedia 
        # https://en.wikipedia.org/wiki/Raised-cosine_filter)        
        H[(np.abs(f) <= (1-roll_off)/2/T)]  = 1.0
        
        # define transition part
        if roll_off > 0.0:
            transition = (np.abs(f) <= (1+roll_off)/2/T) & (np.abs(f) > (1-roll_off)/2/T)
            H[transition] = 0.5 * (1 + np.cos((np.pi*T/roll_off) * (np.abs(f[transition])-((1-roll_off)/2/T))))
        
        # root-raised cosine?
        if root_raised:
            H = np.sqrt(H)
        
        # actual filtering
        samples_out = filter_samples(samples, H, domain)
        
 
    
    return samples_out




def moving_average(samples, average=4, domain='freq'):
    """ Filter a given signal with a moving average filter.
    
    In the time domain implementation, a causal impulse response is used, which
    generates a filter group delay of floor(average/2) samples. CAUTION: The 
    number of samples of the output signal differs from the number of samples 
    of the input signal.
    
    In the frequency domain implementation (cyclic convolution), an acausal 
    filter is used which does not generate any group delay.
    
     Parameters
    ----------
    samples : 1D numpy array, real or complex
        input signal.
    average : int, optional
        number of samples to average. The default is 4.    
    domain : string, optional
        implementation of the filter either in 'time' or in 'freq' domain. 
        The default is 'freq'.

    Returns
    -------
    samples_out : 1D numpy array, real or complex
        filtered output signal.
    
    """
    
    if samples.ndim != 1:
        raise ValueError('signal vector has to be a 1D array...')
    
    if average == 1:
        return samples
    
    # time domain implementation
    if domain == 'time':
        # generate impulse response
        h = np.ones(average)/average
        # actual filtering
        samples_out = filter_samples(samples, h, domain)
    elif domain == 'freq':
        # generate causal impulse response
        h = np.zeros(samples.shape)
        h[0:np.int(np.ceil(average/2))] = 1/average
        h[-np.int(np.floor(average/2))::] = 1/average
        # calc frequency response
        H = np.fft.ifftshift(np.fft.fft(h))
        samples_out = filter_samples(samples, H, domain)    
    
    return samples_out


def windowed_sinc(samples, fc=0.5, order=111, window=None):
    """
    Filter a given signal with a windowed Si-funtion as impulse response.

    Parameters
    ----------
    samples : 1D numpy array, real or complex
        input signal.
    fc : float, optional
        cut off frequency, 0.0 <= fc <= 1.0, where 1.0 specifies the Nyquist 
        frequency (half the sampling frequency). The default is 0.5.
    order : int, optional
        length of the inpulse response, which has to be odd. If order=-1, the
        length is chosen to be the length of the input signal. The default is 111.
    window : string, optional
        type of the window funtion which is multiplied with the sinc-impulse response before filtering.
        'none' --> Si impulse response
        'Hamming'
        'Blackmann-Harris'. The default is None.
            
    CAUTION: This filter generates a group delay which is equal to ceil(order/2).
    Further, the number of samples of the output signal differs from the 
    number of samples of the input signal (see scipy convolution).

    Returns
    -------
    samples_out : 1D numpy array, real or complex
        filtered output signal of lenth samples.size+order-1.

    """  
    
    if samples.ndim != 1:
        raise ValueError('signal vector has to be a 1D array...')
    
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
    else:
        raise ValueError('window has to be either None, "Hamming" or "Blackman-Harris"')
        
    # # debug plots...
    # plt.figure(1)
    # f = np.fft.fftshift(np.fft.fftfreq(order))
    # plt.plot(f, np.abs(np.fft.fftshift(np.fft.fft(h))))
    # plt.show()
    
    # actual filtering
    samples_out = filter_samples(samples, h, domain='time')
    
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
        cut off frequency, 0.0 < fc <= 1.0, where 1.0 specifies the Nyquist 
        frequency (half the sampling frequency). The default is 0.5.
    
    Returns
    -------
    results : dict containing following keys
        samples_out : 1D numpy array, real or complex
            filtered output signal.
        real_fc : float
            actual applied cut off frequency (matching the frequency grid of the FFT)
    """ 
    # generate ideal frequency response
    f = np.fft.ifftshift(np.fft.fftfreq(samples.size, d=0.5))
    H = np.zeros_like(samples)
    H[np.abs(f) <= fc] = 1
    
    real_fc = np.max(np.abs(f[np.abs(f) <= fc]))    
    
    # actual filtering
    samples_out = filter_samples(samples, H, domain='freq')
    
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


def filter_arbitrary(samples,FILTER,sample_rate=1.0):
    """
    Filters a given signal with an arbitrary filter defined in a table.
    This filter is implemented in frequency domain (cyclic convolution).
    The filter can be defined as single-sided frequency axis (real filter response is assumed),
    or with double-sided frequency axis (complex filter response is assumed).

    Parameters
    ----------
    samples : 1D numpy array, real or complex input signal.
        shape = (N,)
    sample_rate : float
        Sample rate of the signal in Hz. The default is 1.0.
    FILTER : 2D numpy array
        1st col.: frequency in Hz
        2nd col.: magnitude in dB
        3rd col.: phase in deg; optional (zero phase assumed if 3rd col. is omitted)

    Returns
    -------
    samples_out : 1D numpy array, real or complex
        filtered output signal.

    """
    # convert Nx1 matrix to vector (see https://stackoverflow.com/questions/39549331/reshape-numpy-n-vector-to-n-1-vector?rq=1)
    samples = samples.reshape(-1,)
    
    # check if FILTER has two or three columns
    
    if isinstance(FILTER,np.ndarray):
    
        if FILTER[0].shape == (3,):
                FILTER = FILTER
        elif FILTER[0].shape == (2,):
            #FILTER = np.append(FILTER, np.zeros((FILTER.shape[0], 1), dtype=FILTER.dtype), axis=1)
            FILTER = np.column_stack((FILTER, np.zeros((FILTER.shape[0], 1), float)))  # add the third column with phase in it                     
        else:
            raise ValueError('FILTER must be a NX2 or NX3 numpy array')
    else :
        TypeError('FILTER must be a Numpy ND-array')
           
    ###### Sort the data in which the frequency is in ascending order
    u,udx = np.unique(FILTER[:,0],return_index=True) # unique ascending frequencies
    FILTER_sorted = FILTER[udx,:] # adds sorted and unique 1st row into the array table
    
    f_Hz = FILTER_sorted[:,0] # freq axis from the table
    
    # fill the FILTER-sorted with mag_lin and phase_rad in the 2nd adn 3rd column
    #FILTER_sorted[:,1] = 10**(FILTER_sorted[:,1]/20)      # Magnitude dB into linear
    FILTER_sorted[:,2] = np.pi/180 * FILTER_sorted[:,2]   # phase angle from degree to radian
    
    
    if all(f_Hz>=0):# check if all frequencies are positive: condition for real filter
        filter_is_real = True
        FILTER_flip= np.flip(FILTER_sorted,0).copy() # flipping  all the arrays columnwise vertically;
        FILTER_flip[:,0]= -FILTER_flip[:,0]          # set all frequencies to negative and adds to FILTER_flip
        FILTER_flip[:,2]= -FILTER_flip[:,2]           # set all phase to negative and adds to FILTER_flip
        
        H = np.concatenate(((FILTER_flip,FILTER_sorted)))  # H = final freq. dom. filter array; concatenate FILTER_flip and FILTER_sorted, generates table of Nx3 order
        
        
        # if f_Hz[0]==0:
        i,idx = np.unique(H[:,0],return_index=True) # dicard double zero  frequency entries if any
        H = H[idx,:] 
            
    else: # doublesided(complex) filter definition
        filter_is_real = False # check correctly for real filter (complex conj. filter definition)
        H = FILTER_sorted
    
    #phase unwrap 
    H[:,2] = np.unwrap(H[:,2])
    
    
    ### Interpolator for magnitude and phase 
    # f_mag = interp1d(H[:,0],H[:,1],bounds_error=False,fill_value='extrapolate', kind='linear') # magnitude interpolation
    # f_ph = interp1d(H[:,0],H[:,2],bounds_error=False,fill_value='extrapolate', kind='linear') # phase interpolation
    f_mag = interp1d(H[:,0],H[:,1],bounds_error=False,fill_value=(H[0,1],H[-1,1]), kind='linear') # magnitude interpolation
    f_ph  = interp1d(H[:,0],H[:,2],bounds_error=False,fill_value=(H[0,2],H[-1,2]), kind='linear') # phase interpolation
    
    # Frequency from input signal,FFT frequencies
    f_sig = np.fft.fftshift(np.fft.fftfreq(samples.size,1/sample_rate)) 
    
    #interpolation to input signal frequency axis
    f_mag_ip = 10**(f_mag(f_sig)/20)   # back to linear scale
    f_ph_ip = f_ph(f_sig) # with phase
    
    # FFT of the input samples
    X_f = np.fft.fft(samples) 
    
    #Calculates Transfer function H(f) and Inverse FFTSHIFT(!!!Important)
    # H_f = (H[:,1])*np.exp(1j*H[:,2])  # H*e^j*phi
    H_interp = np.fft.ifftshift(f_mag_ip * np.exp(1j*f_ph_ip)) 
        
    ### Multiplying interpolated Transfer function and FFT of the input signal
    Y_f = X_f * H_interp
    
    ### IFFT : back to time domain
    samples_out = np.fft.ifft(Y_f)
    
    ### check for real output samples
    if all(np.isreal(samples)) and filter_is_real:
        samples_out = np.real(samples_out)
    
    return samples_out

    
    