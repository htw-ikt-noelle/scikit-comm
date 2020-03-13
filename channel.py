import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
from . import utils

def set_snr(samples, snr_dB=10, sps=1, seed=None):
    """ Set SNR for given sampples.
	
	
	
    """
    if samples.ndim > 2:
        raise ValueError('number of dimensions of samples should be <= 2')
        
    if samples.ndim == 1:
        # promote to 2D array for processing
        samples = samples[np.newaxis, :]   
        
    snr = 10**(snr_dB/10)
    
    power_samples = np.mean(np.abs(samples)**2, axis=-1)
    power_noise = (power_samples / snr * sps).reshape(-1,1)
    
    rng = np.random.default_rng(seed=seed)
    
    # check for real or complex samples
    if np.all(np.isreal(samples)):
        noise = np.sqrt(power_noise) * rng.standard_normal(size=samples.shape)
    else:
        noise = np.sqrt(power_noise/2) * (rng.standard_normal(size=samples.shape) + 
                                          1j * rng.standard_normal(size=samples.shape))
        
    samples_out = samples + noise
    
    return samples_out