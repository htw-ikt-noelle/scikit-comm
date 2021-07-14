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



def sampling_phase_adjustment(samples, sample_rate=1.0, symbol_rate=2.0, shift_dir='both'):
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
    shift_dir : string, optional
        defines the bahaviour of the compensation of the found temporal shift. Can be 
        either 'delay', 'advance' or 'both'. 
        'delay' / 'advance' will only compensate for the time shift by exclusively 
        shifting the signal to the right / left. This means that the time shift 
        can be in the range between [-1/symbol_rate, 0] or [0, 1/symbol_rate], respectively. 
        The option 'both' will shift the signal either to the left or to the 
        right, depending on the minimum amount of time shift. The time shift will 
        therefore be in the range between [-0.5/symbol_rate, 0.5/symbol_rate].

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
    
    # ensure to only advance the signal (shift to the left)
    if (shift_dir == 'advance') and (est_phase < 0.0):
        est_phase += 2 * np.pi        
    
    # ensture to only delay the singal (shift to the right)
    if (shift_dir == 'delay') and (est_phase > 0.0):
        est_phase -= 2 * np.pi        
    
    est_shift = est_phase / 2 / np.pi / symbol_rate
    
    # # for debugging purpose
    # sps = int(sample_rate / symbol_rate)
    # visualizer.plot_eye(samples, sample_rate=sample_rate, bit_rate=symbol_rate)
    # plt.plot(np.abs(samples[:10*sps]), 'C0')
    
    # compensate for found sample phase offset
    samples_out = filters.time_shift(samples, sample_rate, -est_shift)
    
    # # for debugging purpose
    # visualizer.plot_eye(samples_out, sample_rate=sample_rate, bit_rate=symbol_rate)
    # plt.plot(np.abs(samples[:10*sps]), 'C1')
    
    # generate results dict
    results = dict()
    results['samples_out'] = samples_out
    results['est_shift'] = est_shift
    
    return results

def sampling_clock_adjustment(samples, sample_rate=1.0, symbol_rate=2.0, block_size=500):
    """
    Estimate sampling frequency offset and correct for it.
    
    This function estimates the sampling frequency offset between nominal (given) 
    and actual (sample rate of samples) sampling frequency and corrects for the
    found value.
    This is done by splitting the singal into blocks of block_size SYMBOLS and 
    estimating the sampling time offset for each of the blocks (see comm.rx.sampling_phase_adjustment).
    Then, the sampling time offest is corrected for by cyclic time shift of the individual blocks.
    Longer block sizes enhance the accuracy of the time offset estimatin, but cause the sampling
    clock offset to be larger within one block.
    CAUTION: method will fail for large sampling clock (frequency) offsets and if the sampling frequency
    offset leads to a timing error larger than +- 1/symbol_rate over the whole samples.
    

    Parameters
    ----------
    samples : 1D numpy array, real or complex
        input signal.        
    sample_rate : float, optional
        sample rate of input signal in Hz. The default is 1.0.
    symbol_rate : float, optional
        symbol rate of input signal in Hz. The default is 2.0.
    block_size : int, optional
        signal is split into blocks of block_size SYMBOLS. The default is 500.

    Returns
    -------
    results : dict containing following keys
        samples_out : 1D numpy array, real or complex
            output signal.
        est_shift : float or 1D numpy array of floats
            estimated (and inversely applied) temporal shift per block.

    """
    
    # if only one block is questioned -> do only time shift
    if block_size == -1:
        results = sampling_phase_adjustment(samples, sample_rate=sample_rate,
                                                         symbol_rate=symbol_rate, 
                                                         shift_dir='both')
        return results
    # otherwise, do time shift for every block
    else:
        n_samples = np.size(samples, axis=0)    
        sps = np.int(sample_rate / symbol_rate)
        s_per_block = block_size * sps    
        
        # cut signal, to an integer number of blocks and create array
        n_samples_new = int(np.floor(n_samples / s_per_block) * s_per_block)
        samples_array = np.reshape(samples[:n_samples_new], (-1, s_per_block))
        
        tmp = list()
        # run sampling_phase_adjustment once per block
        for idx, block in enumerate(samples_array):
            tmp.append(sampling_phase_adjustment(block, sample_rate=sample_rate,
                                                             symbol_rate=symbol_rate, 
                                                             shift_dir='both'))
        
    # generate output dict containing samples and estimated time shifts per block
    results = dict()
    results['samples_out'] = np.asarray([block['samples_out'] for block in tmp]).reshape(-1)
    results['est_shift'] = np.asarray([block['est_shift'] for block in tmp])
    return results



    ########### MULTIPLE, DIFFERENT METHODS TO COMPENSATE SAMPLING CLOCK OFFSETS (have to be tested) ######################
    """
    These could be used instead of the above implemented method (within the "else" branch)
    
    METHOD 1: "sample dropping", prabably better for real time implementation, but not better in first tests
    
        n_blocks = np.size(samples_array, axis=0)        
        correction  = 0
        tmp = list()
        # run sampling_phase_adjustment multiple times
        for block in np.arange(n_blocks):
            # shift estimation, "dry run"
            samples_tmp = samples[block*s_per_block + correction:block*s_per_block + s_per_block + correction]
            tmp = sampling_phase_adjustment(samples_tmp, sample_rate=sample_rate, symbol_rate=symbol_rate, shift_dir='advance')
            # "sample dropping"
            correction += int(tmp['est_shift']//(1/sample_rate))
            # actual phase adjustemnt
            samples_block = samples[block*s_per_block + correction:block*s_per_block + s_per_block + correction]
            tmp.append(sampling_phase_adjustment(samples_block, sample_rate=sample_rate, symbol_rate=symbol_rate, shift_dir='advance'))
            
            
    MATHOD 2: estimate slope of sampling clock offset over blocks and do resampling (only works in case of almost constant sampling frequency missmatch)
        
        tmp = list()
        # run sampling_phase_adjustment multiple times, once for every block
        for idx, block in enumerate(samples_array):
            tmp.append(sampling_phase_adjustment(block, sample_rate=sample_rate, symbol_rate=symbol_rate, shift_dir='both'))
            # print(results[-1]['est_shift'])
        # generate shifts as ndarrays    
        shifts = np.asarray([block['est_shift'] for block in tmp])
        # do unwrap of shifts (has to be converted in rad???), and shifts have to be fftshifted, because of current filter implementation...can be changed in the future???
        shifts = np.unwrap(np.fft.fftshift(shifts) / (0.5/symbol_rate) * 2 * np.pi) * (0.5/symbol_rate) / 2 / np.pi
        # estimate mean sample time offset per block
        sample_time_offset = np.mean(np.diff(shifts))
        
        t_block = s_per_block / sample_rate
        # ratio between nominal and actual sampling frequency
        ratio = t_block / (t_block + sample_time_offset)
        
        t_old = np.arange(n_samples) / sample_rate
        # interpolate signal at different timing / sampling instants, but keep start and end time the same
        n_samples_new = int(np.round(ratio * n_samples))
        t_new = np.linspace(start=t_old[0], stop=t_old[-1], num=n_samples_new, endpoint=True)        
        f = sinterp.interp1d(t_old, samples, kind='cubic')
        samples = f(t_new)
        
        # do clock phase adjumstemnt once after resampling        
        results = sampling_phase_adjustment(samples, sample_rate=sample_rate, symbol_rate=symbol_rate, shift_dir='both')
    """
    #####################################################################################
    

def carrier_phase_estimation_VV(symbols, n_taps=21, filter_shape='wiener', mth_power=4, rho=0.2):
    """
    Viterbi-Viterbi carrier phase estimation and recovery.
    
    This function estimates the phase noise of the carrier by using the Viterbi-
    Viterbi method [1]. Either a rectangular or a Wiener filter shape can be used.    
    
    
    
    [1] A. Viterbi, "Nonlinear estimation of PSK-modulated carrier phase with 
    application to burst digital transmission," in IEEE Transactions on 
    Information Theory, vol. 29, no. 4, pp. 543-551, July 1983, doi: 10.1109/TIT.1983.1056713.
    
    [2] Ezra Ip, Alan Pak Tao Lau, Daniel J. F. Barros, and Joseph M. Kahn, 
    "Coherent detection in optical fiber systems," Opt. Express 16, 753-791 (2008)
    
    [3] E. Ip and J. M. Kahn, "Feedforward Carrier Recovery for Coherent 
    Optical Communications," in Journal of Lightwave Technology, vol. 25, 
    no. 9, pp. 2675-2692, Sept. 2007, doi: 10.1109/JLT.2007.902118.

    Parameters
    ----------
    symbols :  1D numpy array, real or complex
        input symbols.  
    n_taps : int, optional
        Number of symbols to average over. The default is 21.
    filter_shape : string, optional
        Specifies the filter shape: either 'rect' or 'wiener'. The default is 'wiener'.
    mth_power : int, optional
        Specifies the power to which the constellation is raised in order to 
        remove the modulation (needs to be the number of equidistant phase 
        states of the PSK modulation). The default is 4.
    rho : float, optional
        Tuning factor for 'wiener' filter shape, rho>0; rho is the ratio between
        the magnitude of the phase noise variance sigma²_phi and the additive
        noise variance sigma² (for more informations see [2],[3]). The default is 0.2.

    Returns
    -------
     results : dict containing following keys
        rec_symbols : 1D numpy array, real or complex
            recovered symbols.
        est_shift : float or 1D numpy array of floats
            estimated phase noise.

    """
    
    if filter_shape == 'rect':
        # filter the unwraped phase
        phi_est = filters.moving_average(np.unwrap(mth_power * np.angle(symbols)) / mth_power, n_taps, domain='freq')
       
        
    elif filter_shape == 'wiener':        
        a = 1 + rho / 2 - np.sqrt( ( 1 + rho / 2)**2 - 1) # alpha         
        h_wiener = a * rho / (1 - a**2) * a**np.arange(n_taps // 2 + 1) # postive half
        h_wiener = np.concatenate((np.flip(h_wiener[1:]), h_wiener)) # make symmetric
        h_wiener = h_wiener / np.sum(h_wiener) # normalize to unit sum (make unbiased estimator)        
       
        H_wiener = np.fft.ifftshift(np.fft.fft(h_wiener, n=symbols.size))
        # filter the unwraped phase        
        phi_est = filters.filter_samples(np.unwrap(mth_power*np.angle(symbols))/mth_power, H_wiener, domain='freq')
        # undo group delay of Wiener filter
        phi_est = np.roll(phi_est,  -(n_taps//2))
    
    # for QPSK: shift recoverd constellation by pi/4
    if mth_power == 4:
        phase_correction = np.pi/4
    else:
        phase_correction = 0
    
    # actual phase recovery
    rec_symbols = symbols * np.exp(-1j*(phi_est + phase_correction))
    # crop start and end (necessary???)
    rec_symbols = rec_symbols[1*n_taps+1:-n_taps*1] 
    
    
    # generate output dict containing recoverd symbols and estimated phase noise
    results = dict()
    results['rec_symbols'] = rec_symbols
    results['phi_est'] = phi_est
    return results


def calc_evm(symbols, constellation, norm='max'):
    """
    Calculate the error vector magnitude (EVM).
    
    The EVM [1] is calculated for given received modulation symbols considering
    the given ideal constellation points.
    
    Therefore, the received symbols are normalized to the same power as the ideal
    constellation points before the received symbols are decided to these ideal 
    constellation points. 
    NOTE: the error vector is calculated between the received symbols and these
    DECIDED constellation points and not between the received symbols and the 
    ACTUALLY ("really") sent constellations. This method will therefore lead to 
    an optimistic EVM in case of low SNR (and many wrong symbol decisions e.g.
    high BER).  
    
    The EVM Normalization Reference [2] can be specivied as constellation maximum
    'max' or as reference RMS 'rms'.
    
    [1] https://rfmw.em.keysight.com/wireless/helpfiles/89600b/webhelp/subsystems/digdemod/Content/digdemod_symtblerrdata_evm.htm
    [2] https://rfmw.em.keysight.com/wireless/helpfiles/89600b/webhelp/subsystems/digdemod/Content/dlg_digdemod_comp_evmnormref.htm

    Parameters
    ----------
    symbols : 1D numpy array, real or complex
        input symbols. 
    constellation : 1D numpy array, real or complex
        ideal (sent) constellation points. 
    norm : string, optional
        Specifies the EVM Normalization Reference [2] and can either be 
        constellation maximum 'max' or reference RMS 'rms'. The default is 'max'.

    Returns
    -------
    evm : float
        calculated EVM value as ratio (to convert to percent, the ratio has to 
        be multiplied by 100).

    """
    if norm == 'max':
        evm_norm_ref = np.max(np.abs(constellation))
    elif norm == 'rms':
        evm_norm_ref = np.sqrt(np.mean(np.abs(constellation)**2))
            
    
    # normalize received constellation symbols to ideal constellation
    symbols_norm = symbols * np.sqrt(np.mean(np.abs(constellation)**2) / np.mean(np.abs(symbols)**2))
    
    # decide symbols
    symbols_dec = decision(symbols_norm, constellation)
    
    # calc evm
    error = symbols_norm - symbols_dec
    evm = np.sqrt(np.mean(np.abs(error)**2)) / evm_norm_ref
    
    return evm
    


    
    
    
    
    
    
    
    


