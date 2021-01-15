import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
from . import utils
from . import filters



def generate_bits(n_bits=2**15, type='random', seed=None):
    """
    Generate an array of size (n_bits,) binary values.

    Parameters
    ----------
    n_bits : TYPE, optional
        DESCRIPTION. The default is 2**15.
    type : TYPE, optional
        DESCRIPTION. The default is 'random'.
    seed : TYPE, optional
        DESCRIPTION. The default is None.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    bits : TYPE
        DESCRIPTION.

    """
    
    if type == 'random':
        rng = np.random.default_rng(seed=seed)
        bits = rng.integers(0, high=2, size=n_bits, dtype=bool)
    else:
        raise ValueError('type not implemented yet...')
    
    return bits



def mapper(bits, constellation):
    """ Map bits to a given constellation alphabet.
	
	Bits are grouped into blocks of log2(constellation.size) and converted to
	decimals. These decimals are used to index the particular constellation 
	value in the constellation array.
	
    """
    if constellation.ndim > 1:
        raise ValueError('multiple, different constellations not allowed yet...')    
    
    if bits.ndim > 1:
        raise ValueError('number of dimensions of bits should be 1')   
        
    m = int(np.log2(constellation.size))
    
    if bits.shape[0] % m:
        raise ValueError('number of bits mus be an integer multiple of m')
    
    decimals = np.full((int(bits.shape[0]/m),), np.nan)
    
    if m == 1:
        decimals = bits        
    else:
        decimals = utils.bits_to_dec(bits, m)
    
    symbols = constellation[decimals.astype(int)]
    
    return symbols


def pulseshaper(samples, upsampling=2, pulseshape='rc', roll_off=0.2):
        
    if samples.ndim > 1:
        raise ValueError('number of dimensions of samples should be 1...')   
        
    if upsampling%1:        
        raise ValueError('upsampling factor has to be an integer...')   
        # suggestion for the future (non integer upsampling: e.g. 2.34):
        #   1) upsample to the next smaller integer: 2
        #   2) filtering / pulseshaping
        #   3) resample (scipy.singal.resample (FFT)) to the desired upsampling 2.34
        
    if upsampling == 1:        
        return samples
    
    # upsampling (insert zeros between sampling points)
    # changed implementation necessary due to bug in scipy from version 1.5.0
    # samples_up = signal.upfirdn(np.asarray([1]), samples, up=upsampling, down=1)
    tmp = np.zeros((samples.size, upsampling-1))
    samples_up = np.c_[samples, tmp]
    samples_up = np.reshape(samples_up,-1)
    
    # actual pulseshaping filter
    if pulseshape == 'rc':
        samples_out = filters.raised_cosine_filter(samples_up, 
                                                   sample_rate=upsampling, 
                                                   roll_off=roll_off,
                                                   domain='freq')
    elif pulseshape == 'rrc':
        samples_out = filters.raised_cosine_filter(samples_up, 
                                                   sample_rate=upsampling, 
                                                   roll_off=roll_off, 
                                                   root_raised=True,
                                                   domain='freq')
    elif pulseshape == 'rect':
        samples_out = filters.moving_average(samples_up, upsampling, 
                                             domain='time')
    elif pulseshape == 'None':
        samples_out = samples_up
    else:
        raise ValueError('puseshape can only be either rc, rrc, None or rect...')   
        
            
    
    
    
    return samples_out