import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt


def generate_constellation(format='QAM', order=4):
    """ Generate constellation vector for given modulation format.
	
    """
    
    if ((np.log2(order) % 1) != 0) | (order == 1):
        raise ValueError('gen_constellation: order must be a power of two...')
        
    if format == 'QAM':
        if order == 4:
            constellation = np.array([-1-1j, -1+1j, 1-1j, 1+1j])
        else:
            raise ValueError('gen_constellation: not implemented yet...')
    elif format == 'PAM':
        # TODO: not gray coded yet!!!
        constellation = np.arange(0, order, 1)
        constellation = constellation - np.mean(constellation)
    elif format == 'PSK':
        raise ValueError('gen_constellation: not implemented yet...')
    else:
        raise ValueError('gen_constellation: unknown modulation format...')
    
    return constellation


def bits_to_dec(bits, m):
    """ Convert 1D array of bits into 1D array of decimals, using a resolution of m bits.
    
    """
    bits = np.asarray(bits)
    
    if bits.ndim > 1:
        raise ValueError('dimension of bits should be <=1...')    
        
    if bits.size % m:
        raise ValueError('amount of bits not an integer multiple of m...')
    
    bits_reshape = bits.reshape((-1, m))
    bit_val = np.reshape(2**(np.arange(m-1,-1,-1)),(1,-1))
    decimals = np.sum(bits_reshape * bit_val, axis=-1).astype(int)
    
    return decimals



def dec_to_bits(decimals, m):
    """ Convert 1D array of decimals into 1D array of bits, using a resolution of m bits.
    
    """
    decimals = np.asarray(decimals)
    
    if decimals.ndim > 1:
        raise ValueError('dimension of bits should be <=1...')        
    
    bits = np.full((decimals.size,m), np.nan)
    tmp = decimals
    
    for bit in np.arange(m-1,-1,-1):        
        bits[:, bit] = tmp % 2
        tmp = tmp // 2
                        
    bits = bits.reshape(decimals.size*m).astype(int)
    
    return bits

def create_time_axis(sample_rate=1.0, n_samples=1000):
    """
    Generate a time axis array.

    Parameters
    ----------
    sample_rate : float, optional
        Sample rate of the time axis array. The default is 1.0
    n_samples : int, optional
        Length of the time axis in samples. The default is 1000.

    Returns
    -------
    t : 1D array
        Time axis array of length n_samples sampled at equidistant points with
        time difference 1/sample_rate.

    """
    t = np.arange(0, n_samples) / sample_rate
    return t

    