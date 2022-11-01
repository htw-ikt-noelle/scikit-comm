import numpy as np


def set_snr(samples, snr_dB=10, sps=1.0, seed=None):
    """
    Add noise to an array according to a given SNR (in dB).
    
    CAUTION: this function assumes the input signal to be noise free!
    
    If input signal is of type "real" only real noise (with noise power according to 
    SNR) is generated, if signal is of complex type also complex noise (with noise
    power according to SNR/2 in each quadrature) is added.

    Parameters
    ----------
    samples : numpy array, real or complex
        input signal.
    snr_dB : float, optional
        The desired SNR per symbol in dB. The default is 10.
    sps : float, optional
        samples per symbol of the input signal. The default is 1.0
    seed : int, optional
        random seed of the generated noise samples. The default is None.

   
    Returns
    -------
    samples_out : numpy array, real or complex
        output singal with desired SNR (input signal plus random noise samples).

    """
    
    if samples.ndim > 1:
        raise ValueError('number of dimensions of samples should be <= 1')        
        
    snr = 10**(snr_dB/10)
    
    power_samples = np.mean(np.abs(samples)**2, axis=-1)
    power_noise = (power_samples / snr * sps)
    
    rng = np.random.default_rng(seed=seed)
    
    # check, if samples are of complex type to decide if noise should also be complex
    if samples.dtype == complex:
        noise = np.sqrt(power_noise/2) * (rng.standard_normal(size=samples.shape) + 
                                          1j * rng.standard_normal(size=samples.shape))        
    else:
        noise = np.sqrt(power_noise) * rng.standard_normal(size=samples.shape)
        
    samples_out = samples + noise
    
    return samples_out

def add_phase_noise(samples, s_rate=1.0, linewidth=1.0, seed=None):
    """
    Add laser phase noise to complex signal in 1D ndarray 'samples'.
    
    See https://github.com/htw-ikt-noelle/OptischeKommunikationssysteme/blob/master/LaserPhaseNoise.ipynb
    
    TODO: expand to two polatizations (or higher dimension signals)!!!

    Parameters
    ----------
    samples : numpy array, complex
        complex signal.
    s_rate : float, optional
        sample rate of the incoming singal. The default is 1.0.
    linewidth : float, optional
        3 dB bandwidth of the generated phase noise in Hz. The default is 1.0.
    seed : int, optional
        seed of the random number generator. The default is None.

    Returns
    -------
    results : dict containing following keys
        samples : numpy array, complex
            complex singal including phase noise.
        phaseAcc : numpy array, real
            phase noise vector in rad.
        varPN : float
            variance of generated phase noise in rad**2.
    """          
    # helper calculations
    dt = 1/s_rate   # [s] sample interval for discrete phase noise model
    varPN = 2*np.pi*linewidth*dt; #[rad²] phase variance increase after time-step dt;   proportional to linewidth and observation time dt [Barry1990]/eq.(112)
    # phase noise (Wiener) processes
    np.random.seed(seed=seed)
    phaseInc = np.sqrt(varPN)*np.random.normal(loc=0,scale=1,size=np.size(samples,0)); # [rad] ensemble of Gaussian i.i.d. phase increments with variance varPN
    phaseAcc = np.cumsum(phaseInc,0); # [rad] accumulated phase = random walks
    phaseAcc = phaseAcc - phaseAcc[0]    # [rad] rotate (shift) all phase processes back to initial zero phase

    samples = samples * np.exp(1j*phaseAcc); 
    
    results = dict()
    results['samples'] = samples
    results['phaseAcc'] = phaseAcc
    results['varPN'] = varPN
    
    return results

def add_frequency_offset(samples, sample_rate=1.0, f_d = 100e6):
    """
    Add frequency offset to complex signal in 1D ndarray 'samples'.
    
    Parameters
    ----------
    samples : numpy array, complex
        complex signal.
    sample_rate : float, optional
        sample rate of the incoming singal. The default is 1.0.
    f_d : float, optional
        frequency deviation / frequency offset in Hz. The default is 100 MHz.

    Returns
    -------
    results : dict containing following keys
        samples : numpy array, complex
            complex singal including frequency offset.
    """  
    
    #Creating time axis
    t = np.arange(0, np.size(samples)) / sample_rate
    
    #Adding frequency offset   
    samples = samples *  np.exp(1j*2*np.pi*f_d*t)
            
    results = dict()
    results['samples'] = samples
    return results