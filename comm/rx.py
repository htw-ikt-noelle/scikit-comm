import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
from . import utils
from . import filters
from . import visualizer


def demapper(samples, constellation):
    """
    Demap decided samples to bits using a given constellation alphabet.

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



def sampling_phase_adjustment(samples, sample_rate=1.0, symbol_rate=2.0):
    """
    Estimate the sampling phase offset and compensate for it.
    
    The sampling phase offset is estimated by finding the phase of the oszillation
    with the frequency of the symbol rate in the signal abs(samples)**2. This offset
    is compensated for by a temporal cyclic shift of the input signal.
    
    To contain this frequency component, the signal has to be sampled at least 
    with a rate of three times the symbol rate. If the input signal is sampled
    with lower frequency, the signal is temporally upsampled.
    
    Literature: see ???

    Parameters
    ----------
    samples : 1D numpy array, real or complex
        input signal.        
    sample_rate : float, optional
        sample rate of input signal in Hz. The default is 1.0.
    symbol_rate : float, optional
        symbol rate of input signal in Hz. The default is 2.0.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    results : dict containing following keys
        samples_out : 1D numpy array, real or complex
            cyclic shifted output signal.
        est_shift : float
            estimated (and inversely applied) temporal shift.

    """
    
    # do all dsp on a copy of the samples
    samples_tmp = samples
    # sample rate of dsp (must be at least 3 time symbol rate)
    sr_dsp = sample_rate
    
    # if signal has less than three samples per symbol --> oszillation with
    # symbol rate not present in spectrum of abs(signal)
    if sample_rate < (3 * symbol_rate):
        # upsample to 3 samples per symol
        sr_dsp = symbol_rate * 3
        # watch out, that this is really an integer
        len_dsp = sr_dsp / sample_rate * np.size(samples, axis=0)
        if len_dsp % 1:
            raise ValueError('DSP samplerate results in asynchronous sampling of the data symbols')
        samples_tmp = signal.resample(samples_tmp, num=int(len_dsp), window=None)    
    
    # calc length of vector so that spectrum exactly includes the symbol rate
    tmp = np.floor(symbol_rate * np.size(samples_tmp, axis=0) / sr_dsp)
    n_sam = int(tmp / symbol_rate * sr_dsp)
    
    # cut vector to size and take amplitude square
    samples_tmp = np.abs(samples_tmp[:n_sam])**2
    
    # calc phase of frequency component with frequency equal to the symbol rate
    t_tmp = np.arange(n_sam) / sr_dsp
    est_phase = np.angle(np.sum(samples_tmp * np.exp(1j * 2 * np.pi * symbol_rate * t_tmp)))
    est_shift = est_phase / 2 / np.pi / symbol_rate
    
    # # for debugging purpose
    # visualizer.plot_eye(samples, sample_rate=sample_rate, bit_rate=symbol_rate)
    
    # compensate for found sample phase offset
    samples_out = filters.time_shift(samples, sample_rate, -est_shift)
    
    # # for debugging purpose
    # visualizer.plot_eye(samples_out, sample_rate=sample_rate, bit_rate=symbol_rate)
    
    # generate results dict
    results = dict()
    results['samples_out'] = samples_out
    results['est_shift'] = est_shift
    
    return results
    
    
    
    
    
    
    
    
    


