import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
from . import utils


def demapper(samples, constellation):
    """ Demap decided samples to bits using a given constellation alphabet.
	
	samples are compared to a given constellation constellation alphabet
    array and the position of the corresponding constellation (integer) is 
    converted to the corresponding bit value.
    
    TODO: change function so that samples is not allowed to have ndim > 1!!!
	
	"""
    
    samples = np.asarray(samples)
    constellation = np.asarray(constellation)
    
    if constellation.ndim > 1:
        raise ValueError('multiple, different constellations not allowed yet...')
    
    if samples.ndim > 2:
        raise ValueError('number of dimensions of samples should be <= 2')
        
    if samples.ndim == 1:
        # promote to 2D array for processing
        samples = samples[np.newaxis, :]     
    
    decimals = np.full_like(samples.real, np.nan)
    n_bits = np.int(np.log2(constellation.size))  
    bits = np.full((samples.shape[0], samples.shape[1]*n_bits), np.nan)
    
    for idx_row, row in enumerate(samples):
        for idx_const, cost_point in enumerate(constellation):
            decimals[idx_row, row == cost_point] = idx_const
        bits[idx_row] = utils.dec_to_bits(decimals[idx_row], n_bits)    

                    
    # ## TODO: CHECK!!! for higher order constellations!!!
    # bits = np.reshape(bits, (-1,), order='c').astype(int)
    return bits.squeeze()
        
    

def decision(samples, constellation):
    """ Decide samples samples to a given constellation alphabet.
	
	Find for every samples sample the closest constellation point in a
    constellations array and return this value.
    
    TODO: change function so that samples is not allowed to have ndim > 1!!!
	
	"""    
    if constellation.ndim > 1:
        raise ValueError('multiple, different constellations not allowed yet...')    
    
    if samples.ndim > 2:
        raise ValueError('number of dimensions of samples should be <= 2')
        
    if samples.ndim == 1:
        # promote to 2D array for processing
        samples = samples[np.newaxis, :]     
    
    # normalize samples to mean magnitude of original constellation
    mag_const = np.mean(abs(constellation))    
    mag_samples = np.mean(abs(samples), axis=-1).reshape(-1,1)
    samples_norm = samples * mag_const / mag_samples
    
    # shape to 2D array and repeat in order to match size of samples
    const = np.tile(constellation.reshape(-1,1), (1, samples_norm.shape[1]))
    
    dec_symbols = np.full_like(samples_norm, np.nan)
    for row_idx, row in enumerate(samples_norm):
        const_idx = np.argmin(np.abs(row-const), axis=0)
        dec_symbols[row_idx] = constellation[const_idx]
        
    return dec_symbols.squeeze()



def count_errors(bits_tx, bits_rx):
    """ Count bit errors and return the bit error rate.
	
	"""
    
    if (bits_rx.ndim > 2) | (bits_tx.ndim > 2):
        raise ValueError('number of dimensions of bits should be <= 2')
    
    err_idx = np.not_equal(bits_tx, bits_rx)
    ber = np.sum(err_idx, axis=-1) / bits_tx.shape[-1]
    
    return ber, err_idx