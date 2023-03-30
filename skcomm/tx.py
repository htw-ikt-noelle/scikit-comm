import numpy as np
import scipy.signal as ssignal

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


def pulseshaper(samples, upsampling=2.0, pulseshape='rc', roll_off=0.2):
    """
    Upsample and pulseshape a given sample sequence.
    
    The provided samples are upsampled by the factor upsampling.
    
    This is done by inserting ceil(upsampling)-1 zeros between each sample followed 
    by applying a pulseshaping filter. (Root) raised cosine and rectangular filter
    impulse respnses are available. pulseshape = None does not apply any filter 
    to the upsampled sequence.
    
    After the pulseshaping a resampling (downsampling) is performed in case of a
    fractional upsampling factor.    

    Parameters
    ----------
    samples : 1D numpy array, real or complex
        input signal.
    upsampling : float, optional
        upsampling factor. The default is 2.0
    pulseshape : sting, optional
        pulseshaping filter, can either be 'rc', 'rrc', 'rect' or 'None', 
        meaning raised cosine filter, root raised cosine filter, erctangular 
        filter or no filter, respectively. The default is 'rc'.
    roll_off : float, optional
        rolloff factor in case of (root) raised consine filter. The default is 0.2.

    
    Returns
    -------
    samples_out : 1D numpy array, real or complex
        upsampled and pulseshaped signal samples.

    """ 
    if samples.ndim > 1:
        raise ValueError('number of dimensions of samples should be 1...') 
        
    # integer upsampling factor before pulseshaping
    upsampling_int = int(np.ceil(upsampling))
        
    if upsampling%1:        
        # remaining fractional factor for downsampling after pulseshaping
        resampling_rem = upsampling / upsampling_int
    else:
        # no resampling necessary
        resampling_rem = None      
       
    if upsampling == 1:
        # shortcut
        return samples
    
    # upsampling (insert zeros between sampling points)
    # changed implementation necessary due to bug in scipy from version 1.5.0
    # samples_up = ssignal.upfirdn(np.asarray([1]), samples, up=upsampling, down=1)
    tmp = np.zeros((samples.size, upsampling_int-1))
    samples_up = np.c_[samples, tmp]
    samples_up = np.reshape(samples_up,-1)
    
    # check if symbols are real
    if np.isrealobj(samples_up):
        real = True
    else:
        real = False
    
    # actual pulseshaping filter
    if pulseshape == 'rc':
        samples_out = filters.raised_cosine_filter(samples_up, 
                                                   sample_rate=upsampling_int, 
                                                   roll_off=roll_off,
                                                   domain='freq') * upsampling_int
    elif pulseshape == 'rrc':
        samples_out = filters.raised_cosine_filter(samples_up, 
                                                   sample_rate=upsampling_int, 
                                                   roll_off=roll_off, 
                                                   root_raised=True,
                                                   domain='freq') * upsampling_int
    elif pulseshape == 'rect':
        samples_out = filters.moving_average(samples_up, upsampling_int, 
                                             domain='freq') * upsampling_int
    elif pulseshape == 'None':
        samples_out = samples_up
    else:
        raise ValueError('puseshape can only be either rc, rrc, None or rect...') 
        
    # if symbols are real, the pulseshaped samples should be real as well
    if real:
        samples_out = np.real(samples_out)
        
    if resampling_rem:
        # check for an integer number of samples after resampling
        if (samples_out.size * resampling_rem) % 1:            
            raise ValueError('Length of resampled vector must be an integer. Modify vector length or upsampling factor.')
        # resampling
        samples_out = ssignal.resample(samples_out, int(resampling_rem*len(samples_out)))
    
    return samples_out