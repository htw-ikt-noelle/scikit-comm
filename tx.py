import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
from . import utils



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
    
    if bits.ndim > 2:
        raise ValueError('number of dimensions of bits should be <= 2')
        
    if bits.ndim == 1:
        # promote to 2D array for loop processing
        bits = bits[np.newaxis, :]
        
    m = int(np.log2(constellation.size))
    
    if bits.shape[1] % m:
        raise ValueError('number of bits mus be an integer multiple of m')
    
    decimals = np.full((bits.shape[0], np.int(bits.shape[1]/m)), np.nan)
    
    if m == 1:
        decimals = bits        
    else:
        for idx, row in enumerate(bits):
            decimals[idx] = utils.bits_to_dec(row, m)
    
    symbols = constellation[decimals.astype(int)]
    
    return symbols.squeeze()