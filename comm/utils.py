import numpy as np
import scipy.signal as signal
import scipy.special as sspecial
import matplotlib.pyplot as plt
import math


def generate_constellation(format='QAM', order=4):
    """ Generate constellation vector for given modulation format.
	
    """
    
    if ((np.log2(order) % 1) != 0) | (order == 1):
        raise ValueError('gen_constellation: order must be a power of two...')
    if type(order) != int:
        raise TypeError('Order parameter must be passed as integer.')
        
    if format == 'QAM':
        # derive number of bits encoded in one symbol from QAM order
        n = int(np.log2(order))
        ### build rectangular matrix of gray-coded bit words
        # generate separate Gray codes for I and Q using XOR (if n%2 == 1, n+1 bits
        # will be encoded in the I- and n bits in the Q-branch, respectively, if 
        # n%2 == 0, n bits will be encoded in both the I and Q branch)
        width_bit = int(-(-n // 2))
        height_bit = int(n // 2)
        gray_I_dec = np.arange(2**width_bit) ^ ((np.arange(2**width_bit))>>1)
        gray_Q_dec = np.arange(2**height_bit) ^ ((np.arange(2**height_bit))>>1)
        # generate indices for I and Q values to build matrix to project 
        # constellation points onto later
        x_I = np.arange(int(2**(-(-n // 2))))
        y_Q = np.arange(int(2**(n // 2)))
        # combine into meshgrid
        xx,yy = np.meshgrid(x_I,y_Q)
        # build matrix of decimal values whose binary representations have
        # a Hamming distance of 1 in both vertical and horizontal direction by
        # shifting bits of I-branch Gray code left by floor(n/2) and adding bits 
        # of Q-branch Gray code
        bits = (gray_I_dec[xx]<<int(n // 2)) + gray_Q_dec[yy]
        
        # for odd-bit QAM constellations:
        if int(np.log2(order))%2:
            ### build symmetrical matrix of complex symbols
            # calc dimension for symmetrical constellation matrix
            cross_dim = int((2**(-(- n // 2)) + 2**(n // 2)) / 2)
            # generate evenly spaced values with euclidean distance of 2
            cross_values = np.linspace(-cross_dim + 1, cross_dim - 1, cross_dim)
            # generate meshgrid
            cross_I, cross_Q = np.meshgrid(cross_values,cross_values)
            # build complex symbols
            symbols = cross_I + 1j*cross_Q
            # cut away corners
            cut = int(cross_dim / 6)
            symbols[:cut,:cut] = 0
            symbols[-cut:,:cut] = 0
            symbols[:cut,-cut:] = 0
            symbols[-cut:,-cut:] = 0
            # copy matrix for assigning gray-coded decimal values to
            bits_symm = np.full_like(symbols,0,dtype='int')
            # write 'middle block' of rectangular constellation into symmetrical
            # matrix
            bits_symm[cut:-cut,:] = bits[:,cut:-cut]
            # manipulate the 8 'end blocks' of rectangular constellation and 
            # write them into new positions in the symmetrical matrix
            # top left block
            bits_symm[:cut,cut:2*cut] = np.flipud(bits[:cut,:cut])
            # upper middle left block
            bits_symm[:cut,2*cut:3*cut] = np.fliplr(bits[cut:2*cut,:cut])
            # lower middle left block
            bits_symm[-cut:,2*cut:3*cut] = np.fliplr(bits[-(2*cut):-cut,:cut])
            # bottom left block
            bits_symm[-cut:,cut:2*cut] = np.flipud(bits[-cut:,:cut])
            # top right block
            bits_symm[:cut,-(2*cut):-cut] = np.flipud(bits[:cut,-cut:])
            # upper middle right block
            bits_symm[:cut,-(3*cut):-(2*cut)] = np.fliplr(bits[cut:2*cut,-cut:])
            # lower middle right block
            bits_symm[-cut:,-(3*cut):-(2*cut)] = np.fliplr(bits[-(2*cut):-cut,-cut:])
            # bottom right block
            bits_symm[-cut:,-(2*cut):-cut] = np.flipud(bits[-cut:,-cut:])
            
            ### manipulate and reshape symmetrical matrix into array of connstellation points
            # flatten matrices out and delete entries at indices where
            # cross_symbols == 0 (the corners that were cut away)
            bits_symm = np.delete(bits_symm.flatten(),np.argwhere(symbols.flatten()==0))
            symbols = np.delete(symbols.flatten(),np.argwhere(symbols.flatten()==0))
            # write into bits for naming convention
            bits = bits_symm
            
        # for even-bit QAM
        else:
            # generate evenly space values for I and Q and build matrix of complex
            # symbols
            values_I = np.linspace(-(np.sqrt(order))+1,np.sqrt(order)-1,int(np.sqrt(order)))
            values_Q = np.linspace(-(np.sqrt(order))+1,np.sqrt(order)-1,int(np.sqrt(order)))
            II,QQ = np.meshgrid(values_I,values_Q)
            symbols = (II + 1j*QQ).flatten()
            
        
        # convert gray-coded sequence to binary for control and labelling purposes
        # change dtype if more than 8 bits are encoded per symbol
        bits_bin = np.full_like(bits.flatten(),0,dtype='<U8')
        for i in range(len(bits.flatten())):
            bits_bin[i] =  np.binary_repr(bits.flatten()[i], width=int(n))
        ### sort and output values as numpy arrays
        # initialize lists for return values
        constellation = []
        bits_tmp = []
        # iterate over flattened symbols and bits matrices and append complex
        # constellation points and binary number labels to respective lists
        for i in range(order):
            constellation.append(symbols[np.argwhere(bits.flatten() == i)][0][0])
            bits_tmp.append(bits_bin.flatten()[np.argwhere(bits.flatten() == i)][0][0])
        # convert into arrays
        bits = np.asarray(bits_tmp)
        constellation = np.asarray(constellation)
    ### PAM    
    elif format == 'PAM':
        # TODO: not gray coded yet!!!
        constellation = np.arange(0, order, 1)
        constellation = constellation - np.mean(constellation)
    ### ASK    
    elif format == 'ASK':
        raise ValueError('ASK not implemented yet')
    ### PSK    
    elif format == 'PSK':
        n = int(np.log2(order))
        # generate Gray code
        gray = np.arange(2**n) ^ (np.arange(2**n)>>1)
        gray_bin = np.full_like(gray,0,dtype='<U8')
        for i in range(len(gray)):
            gray_bin[i] = np.binary_repr(gray[i],width=int(n))
        # build constellation points
        symbols = np.asarray([np.exp(1j*2*np.pi*i/order) for i in range(len(gray))]) 
        # reorder symbols and label vector
        constellation = []
        bits = []
        for i in range(order):
            constellation.append(symbols.flatten()[np.argwhere(gray.flatten()==i)][0][0])
            bits.append(gray_bin.flatten()[np.argwhere(gray.flatten()==i)][0][0])
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
                        
    bits = bits.reshape(decimals.size*m).astype(bool)
    
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

def ber_awgn(value=np.arange(20), modulation_format='QAM', modulation_order=4, type='SNR', symbol_rate=32e9, d_lambda=0.1e-9, ref_wl=1550e-9, PDM=False):
    """
    Calculate theoritical BER performances in case of additive white Gaussian 
    noise (AWGN) for various modulation formats.
    
    Amount of noise can either be specified as signal-to-noise ratio (SNR), 
    signal-to-noise ratio per bit (SNRB) or as
    optical signal-to-noise ratio (OSNR), respectively.
    NOTE: The parameter value is interpreted as logarithmic scale (dB).
       
    References:
        
    [1] Essiambre et al., JLT, vol 28, no. 4, 2010, "Capacity Limits of Optical Fiber Networks"
    [2] Xiong, 2006, "Digital Modulation Techniques", second edition, Artech House
    
   

    Parameters
    ----------
    value : 1D array, float
        range which specifies the amount of noise for which the BER perfomance 
        is calculated. 
        NOTE: The value is interpreted as logarithmic value (dB). The default
        is np.arange(20).
    modulation_format : string, optional
        modulation format for which the BER performance is calculated. Can be 
        'QAM'. The default is 'QAM'.
    modulation_order : int, optional
        modulation order (number of bits per symbol) for which the BER performance 
        is calculated. Has to be a power of 2. The default is 4.
    type : string, optional
        specifies the type how the parameter value is interpreted. Can either 
        be 'SNR', SNRB or 'OSNR'. The default is 'SNR'.
    symbol_rate : float, optional
        symbol rate of the signal. Is used to convert from OSNR to SNR and 
        vice versa. Only affects the result in case of tpye='OSNR'. 
        The default is 32e9.
    d_lambda : float, optional
        bandwidth (in m) which is used to calculate the OSNR. Only affects 
        the result in case of tpye='OSNR'. The default is 0.1e-9.
    ref_wl : float, optional
        center wavelength of the optical signal. Only affects the result in 
        case of tpye='OSNR'. The default is 1550e-9.
    PDM : bool, optional
        is the optical signal a polarization-division multiplexed (PDM) 
        signal? Only affects the result in case of tpye='OSNR'. The default is False.

    
    Returns
    -------
    ber : 1D array, float
        theoretical BER performance for the specified amount of AWGN.

    """    
    # reference bandwidth to convert OSNR to SNR and vice versa, see [1], eq. (34)
    # 12.5GHz corresponds roughly to 0.1nm at 1550nm wavelength
    # TODO: add utils df_dlam
    # speed of light
    c0 = 299792458
    b_ref = d_lambda*c0/ref_wl**2
    
    # PDM scaling factor, see [1], eq (34)
    if PDM:
        p = 2
    else:
        p = 1
    
    # convert given value to snr_lin in order to calculate BER performance
    if type == 'OSNR':
        osnr_dB = value
        osnr_lin = 10**(value/10)
        # see [1] equation (34)
        snr_lin = osnr_lin * 2 * b_ref / p / symbol_rate
    elif type == 'SNR':
        snr_dB = value
        snr_lin = 10**(value/10)
    elif type == 'SNRB':
        snrb_dB = value
        snrb_lin = 10**(snrb_dB/10)
        # see [1], equation (31)
        snr_lin = snrb_lin * math.log2(modulation_order)
    else:
        raise ValueError('type should be "OSNR", "SNR" or "SNRB".')  
        
    ber = []
    
    if modulation_format == 'QAM':
        #  calc ber performance for M-QAM according to [2], p. 465, nonnumbered eq.
        for snr in snr_lin:
            s = 0.0
            for i in range(1,int(math.sqrt(modulation_order)/2)+1):
                # calc inner sum and use correspondence Q(x)=0.5*erfc(x/sqrt(2))
                s += 0.5 * sspecial.erfc((2*i-1) * math.sqrt(3*snr/(modulation_order-1))/math.sqrt(2)) 
            tmp = 4 / math.log2(modulation_order)*(1-1/math.sqrt(modulation_order)) * s
            ber.append(tmp)
    else:
        raise ValueError('modulation format not implemented yet.')  
    
    return np.asarray(ber)
        
    