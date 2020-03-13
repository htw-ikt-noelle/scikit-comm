import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
from . import utils

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