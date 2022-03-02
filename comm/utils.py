import numpy as np
import scipy.signal as signal
import scipy.special as sspecial
import matplotlib.pyplot as plt
import math
import time
import json
import pickle
import csv
from numpy.polynomial import Polynomial


def generate_constellation(format='QAM', order=4):
    """
    Generate array of Gray-coded constellation points for a given modulation 
    format of a given order.
	
    For QAM formats, a pure Gray-coded square QAM constellation is generated 
    if order is an even power of 2 (even-bit), and Pseudo-Gray-coded symmetrical
    QAM constellation based on the proposed implementation in [1] is generated 
    if order is an odd power of 2 (odd-bit) and order is 32 or greater. For the
    special case of order = 8 (3-bit), an 8-Star-QAM is implemented.
    For PSK and PAM formats, pure Gray-Coded constellations are always generated, 
    as long as order is a power of 2.
    
    Constellation points in QAM and PAM modulation schemes are currently built
    with a Euclidean distance of 2 between two neighbouring points, as per the
    convention in [1]. This might require normalization of the constellation
    array.
    
    [1] J. G. Smith, "Odd-Bit Amplitude Quadrature-Shift Keying", IEEE Transactions
    on Communications, pp. 385-389, 1975
    
    Parameters
    ----------
    format : string, optional
        Modulation format to be used. Options are 'QAM', 'PAM', and 'PSK'.
        The default is 'QAM'.
    order : int, optional
        Order of the modulation scheme to be generated. Only powers of 2 are
        valid. The default is 4.

    Raises
    ------
    ValueError : 
        If order is not a power of 2 or if a QAM with order 2 is generated 
        (for order = 2, PAM/PSK should be used).
    TypeError : 
        If order is not passed as integer.

    Returns
    -------
    constellation : ndarray of complex
        Array of complex constellation points to map bit words to.

    """


    
    if ((np.log2(order) % 1) != 0) | (order == 1):
        raise ValueError('gen_constellation: order must be a power of two...')
    if type(order) != int:
        raise TypeError('gen_constellation: order must be passed as integer...')
    
    #### QAM
    if format == 'QAM':
        # check for reasonable order for QAM (for order == 2, PSK or ASK should
        # be used, for order == 8, Star QAM is implemented non-algorithmically)
        if order == 2:
            raise ValueError('gen_constellation: for order == 2, use PSK or ASK instead of QAM...')
        if order == 8:
            constellation = np.asarray([-1-1j, 1-1j, -1+1j, 1+1j, -1j*(1+np.sqrt(3)), 1+np.sqrt(3,), -1-np.sqrt(3), 1j*(1+np.sqrt(3))])
            return constellation
        
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
            
        ### sort and output values as numpy arrays
        # convert gray-coded sequence to binary for control and labelling purposes
        # change dtype if more than 8 bits are encoded per symbol
        bits_bin = np.full_like(bits.flatten(),0,dtype='<U8')
        for i in range(len(bits.flatten())):
            bits_bin[i] =  np.binary_repr(bits.flatten()[i], width=int(n))
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
    #### PAM    
    elif format == 'PAM':
        # https://electronics.stackexchange.com/questions/158754/what-is-the-difference-between-pam-and-ask
        n = int(np.log2(order))
        gray = np.arange(2**n) ^ (np.arange(2**n)>>1)
        symbols = np.linspace(-(2**n)+1,(2**n)-1,order)
        constellation = []
        for i in range(order):
            constellation.append(symbols[np.argwhere(gray==i)][0][0])
        constellation = np.asarray(constellation)
    #### PSK    
    elif format == 'PSK':
        # hardcoded BPSK for better accuracy of the constellation points
        if order == 2:
            constellation = np.asarray([-1+0j,1+0j])
        # other PSK orders
        else:
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
            constellation = np.asarray(constellation)
            bits = np.asarray(bits)
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


def osnr(power_vector = [], wavelength_vector = [], interpolation_points = [], integration_area = [], resolution_bandwidth = 0.1, polynom_order = 3, plotting = False):

    """ 
    Description
    -----------
    Function to calculate the OSNR from OSA (Optical spectrum analyzer) trace data via interpolation method. The function will interpolate the spectral noise shape 
    in a given spectral area, which can be definded by the user. From this data the noise power is estimated. Than the function will calculate the signal power and 
    will afterwards calculate the OSNR. 


    Paramters
    ---------
        power_vector: numpy array
            Vector with the power values of the OSA trace. 
            Must be dBm.
            Must be same length as wavelength_vector.

        wavelength_vector: numpy array
            Vector with the wavelength values of the OSA trace.
            Must be nm.
            Must be same length as power vector.

        interpolation_points: numpy array of length 4 [a,b,c,d]
            This array specifies the areas for creating the polynomial. This requires 4 points. The array elements a and b indicate the left area of ​​
            the signal spectrum and the elements c and d the right area.

            If the passed wavelength value is not present in the wavelength vector, the passed values ​​are rounded to the nearest existing value.

        integration_area: numpy array of length 2 [integration_start, integration_stop]
            These two points determine the bandwidth in which the noise and signal power are determined.

            If the passed wavelength value is not present in the wavelength vector, the passed values ​​are rounded to the nearest existing value.

        resolution_bandwidth: float
            Insert here the used resolution bandwidth (rbw) of the OSA.

        polynom_order: int
            Insert here the polynomial order for the noise interpolation.

        plotting: boolean, optional (default = False)
            If true, the spectrum is plotted with the interpolation area, integration area and interpolated noise shape. 
            To show the plot, plt.show() must be called in the main script. 

    Returns
    -------
        OSNR_01nm:
            The calculated OSNR normalized to a noise bandwidth of 0.1nm.

        OSNR_val: 
            The calculated OSNR of the integration area.


    Examples
    --------
        >>> import comm as comm
        >>> import numpy as np

        # Set area for polynom creation (Values were randomly selected for this example)
        >>> a = 1552.025
        >>> b = 1552.325
        >>> c = 1552.725
        >>> d = 1553.025

        # Set integration area (Values were randomly selected for this example)
        >>> integration_start = 1552.375
        >>> integration_stop = 1552.675

        # Set polynomial order
        >>> poly_ord = 2

        # Get optical spectrum data from OSA or another arbitary source
        >>> OSA_trace_dict = comm.instrument_control.get_samples_HP_71450B_OSA()
        >>> power = OSA_trace_dict['A']['Trace_data']
        >>> wavelength = OSA_trace_dict['A']['WL_Vector']
        >>> resolution_bw = OSA_trace_dict['A']['Resolution_BW']*1e9

        # Calculate OSNR with plot
        >>> [ONSR_0.1nm,OSNR] = comm.osnr.osnr(power_vector = power,
                            wavelength_vector = wavelength,
                            interpolation_points = np.array([a,b,c,d]),
                            integration_area = np.array([integration_start,integration_stop]),
                            resolution_bandwidth = resolution_bw,
                            polynom_order=poly_ord,
                            plotting = True)

    """

    # =============================================================================
    #  Check inputs of correctnes
    # ============================================================================= 

    try:
        if not (isinstance(power_vector, np.ndarray) and isinstance(wavelength_vector, np.ndarray)
           and isinstance(interpolation_points, np.ndarray) and isinstance(integration_area, np.ndarray)):
            raise TypeError('power_vector, wavelength_vector, interpolation_points or integration are are not of type np.array')

        if not (isinstance(resolution_bandwidth, float)):
            raise TypeError("resolution_bandwidth must be float")

        if not (isinstance(polynom_order, int)):
            raise TypeError("polynom_order must be int")

        if not (isinstance(plotting, bool)):
            raise TypeError("plotting must be bool")

        if not (power_vector.size == wavelength_vector.size):
            raise ValueError("power_vector and wavelength_vector must be same size")

        if not (interpolation_points.size == 4):
            raise ValueError("interpolation_points needs 4 elements") 

        if not (integration_area.size == 2):
            raise ValueError("integration_area needs 2 elements") 

        if not (interpolation_points[0] < interpolation_points[1] < interpolation_points[2] < interpolation_points[3]):
            raise ValueError("Values of interpolation_points must meet the following conditions: a < b < c < d")

        if not (integration_area[0] < integration_area[1]):
            raise ValueError("Values of integration_area must meet the following conditions: integration_start < integration_stop")

    except Exception as e:
        print("Error: {0}".format(e))
        exit()
    

    # =============================================================================
    #  Calculations
    # ============================================================================= 


    # Correct the input interpolation points to the nearest wavelength in the wavelength vector
    # - Find the position of these values
    closest_interpolation_wavelength = np.array([])
    closest_interpolation_wavelength_index = np.array([])
    for idx,inter_point in enumerate(interpolation_points):
        difference_vector = np.abs(wavelength_vector - inter_point)
        closest_interpolation_wavelength_index = np.int16(np.append(closest_interpolation_wavelength_index,difference_vector.argmin()))
        closest_interpolation_wavelength = np.append(closest_interpolation_wavelength,wavelength_vector[int(closest_interpolation_wavelength_index[idx])])

    # Correct the input integration area to the nearest wavelength in the wavelength vector
    # - Find the position of these values
    closest_integration_wavelength = np.array([])
    closest_integration_wavelength_index = np.array([])
    for idx,integration_point in enumerate(integration_area):
        difference_vector = np.abs(wavelength_vector - integration_point)
        closest_integration_wavelength_index = np.int16(np.append(closest_integration_wavelength_index,difference_vector.argmin()))
        closest_integration_wavelength = np.append(closest_integration_wavelength,wavelength_vector[int(closest_integration_wavelength_index[idx])])

    # Create the interpolated noise shape
    # Getting the wavelengths between lambda 0 (a) and lambda 1 (b)
    wavelengths_lambda_0_1 = wavelength_vector[closest_interpolation_wavelength_index[0]:closest_interpolation_wavelength_index[1]+1]

    # Getting the wavelengths between lambda 2 (c) and lambda 3 (d)
    wavelengths_lambda_2_3 = wavelength_vector[closest_interpolation_wavelength_index[2]:closest_interpolation_wavelength_index[3]+1]

    # Combine the both wavelengths vectors into one vector.
    sample_point_wavelengths_vector = np.append(wavelengths_lambda_0_1,wavelengths_lambda_2_3)

    # Getting the power between lambda 0 and lambda 1
    power_lambda_0_1 = power_vector[closest_interpolation_wavelength_index[0]:closest_interpolation_wavelength_index[1]+1]

    # Getting the power between lambda 2 and lambda 3
    power_lambda_2_3 = power_vector[closest_interpolation_wavelength_index[2]:closest_interpolation_wavelength_index[3]+1]

    # Combine the both power vectors into one vector.
    sample_point_power_vector = np.append(power_lambda_0_1,power_lambda_2_3)

    # Creation of the polynom
    # - The Polynomial.fit method will give back an scaled version of the the coefficients back. 
    # - To get unscaled values, the convert() method is needed.
    polynom_coeffs = Polynomial.fit(sample_point_wavelengths_vector,sample_point_power_vector,polynom_order).convert().coef

    # Creation of polynomial
    poly = Polynomial(polynom_coeffs)
    interpol_noise_powers_complete = poly(wavelength_vector)

    # For OSNR calculation needed span
    interpol_noise_powers = interpol_noise_powers_complete[closest_integration_wavelength_index[0]:closest_integration_wavelength_index[1]+1]

    # To calculate the power the power values must be converted from db to linear.
    delta_lambda = np.diff(wavelength_vector)
    delta_lambda = np.append(delta_lambda,delta_lambda[-1])
    power_vector_lin = 10**np.divide(power_vector,10) * delta_lambda / resolution_bandwidth
    interpol_noise_powers_lin = 10**np.divide(interpol_noise_powers,10) * delta_lambda[closest_integration_wavelength_index[0]:closest_integration_wavelength_index[1]+1] / resolution_bandwidth
    
    # Calculation noise power
    pseudo_noise_power = np.trapz(interpol_noise_powers_lin,wavelength_vector[closest_integration_wavelength_index[0]:closest_integration_wavelength_index[1]+1])

    # Calculation signal plus noise power
    pseudo_signal_noise_power = np.trapz(power_vector_lin[closest_integration_wavelength_index[0]:closest_integration_wavelength_index[1]+1],
                                        wavelength_vector[closest_integration_wavelength_index[0]:closest_integration_wavelength_index[1]+1])
    # Calculation signalpower
    pseudo_signal_power = pseudo_signal_noise_power - pseudo_noise_power

    # Calculation OSNR
    OSNR_val = 10*np.log10(pseudo_signal_power/pseudo_noise_power)
    bandwidth = closest_integration_wavelength[1] - closest_integration_wavelength[0]
    OSNR_01nm = 10*np.log10(pseudo_signal_power / (pseudo_noise_power * 0.1 / bandwidth))

    if plotting == True:
        plt.figure()
        plt.plot(wavelength_vector,power_vector,'-',
                integration_area,[power_vector[closest_integration_wavelength_index[0]],power_vector[closest_integration_wavelength_index[1]]],'ro',
                np.append(wavelengths_lambda_0_1,wavelengths_lambda_2_3),np.append(power_lambda_0_1,power_lambda_2_3),'g.',
                wavelength_vector,interpol_noise_powers_complete,'-',
                )

        plt.gca().legend(('Optical power from OSA','Integration borders','Area for polyfit','Interpolated noise power' ))
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('Power density [dBm/{0}nm]'.format(resolution_bandwidth))        
        plt.ylim(np.min(power_vector)-10,np.max(power_vector)+10)
        plt.grid()

    return OSNR_01nm,OSNR_val
    
        
def estimate_snr(sig,block_size=-1,bias_comp=False):
    """
    Estimates the SNR per symbol of a signal with AWGN depending on its modulation format.
    Different algorithms are used depending on whether the modulation format is
    BPSK, QPSK, or QAM of order 16 or upwards. The function assumes that the 'samples'
    attribute contains symbols, i.e. that the signal has been downsamples to 1 sample
    per symbol.

    References:
    [1]: Ijaz, A., Awoseyila, A.B., Evans, B.G.: "Improved SNR estimation for BPSK and
    QPSK signals", Electronics Letters Vol. 45 No. 16, 2009
        
    [2]: Qun, X., Jian, Z.: "Improved SNR Estimation Algorithm, International Conference
    on Computer Systems", Electronics, and Control (ICCSEC), 2017
    
    [3]: Xu, H., Li, Z., Zheng, H.: "A non-data-aided SNR Estimation Algorithm for QAM
    Signals", IEEE, 2004

    Parameters
    ----------
    sig : signal-class object
        The signal object on which to operate. The object's modulation format must
        be known.
    block_size : int
        The number of symbols to average SNR over. Greater block size means a 
        more accurate estimate. Default value is -1, which treats the entire 
        symbol vector as one block.
    bias_comp: bool
        Flag determining whether the bias of the employed algorithm is compensated
        according to the data obtained in a Monte Carlo trial in a setup file. 
        Default value is False.

    Returns
    -------
    snr_estimate : numpy array
        Array with the same shape as sig.samples containing the estimated SNR values.

    """
    # TODO: edit generate_constellation function so that the "modulation info"
    # attribute is more descriptive
    
    # TODO: run Monte Carlo trial in setup file and determine approximate bias
    # of each algorithm over a realistic SNR range, then compensate it if flag is set
    
    # TODO: add appropriate type and value error checks
    if type(block_size) != int:
        raise TypeError("Block size must be a positive integer or -1!")
    
    if block_size > 1 & block_size != -1:
        raise ValueError("Block size must be a positive integer or -1!")
    
    # nested function that performs the estimation
    def inner(samples,mod_info):
        
        # BPSK case
        if mod_info == '2-PSK':
            # ref.: [1], Eqs. 6, 17, valid for SNR values between 0 dB and 20 dB
            # calc 2nd order moment
            m2_hat = np.mean(np.abs(samples)**2)
            # calc novel signal power estimate
            s_hat = ((np.sum(np.abs(np.real(samples)))**2) + (np.sum(np.abs(np.imag(samples)))**2)) / (len(samples)**2)
            # calc SNR estimate
            snr_estimate = s_hat / (m2_hat - s_hat)
            
        # QPSK case
        elif mod_info in ['4-PSK','4-QAM']:
            # ref.: [2], Eqs. 11-15, valid for SNR values between -10 dB and 30 dB
            # calc mean
            symb_mean = np.mean(np.abs(samples))
            # calc variance
            symb_var = np.mean((np.abs(samples)-symb_mean)**2)
            # calc SNR estimate
            snr_estimate = 10*np.log10((np.abs(symb_mean)**2)/(2*symb_var))
            # for SNR values below 10 dB, modify the estimate
            if snr_estimate < 10:
                snr_estimate = np.sqrt((snr_estimate-2.5)*39.2)-7
            
        # QAM case (order between 16 and 256, either square or symmetrical QAM)
        elif mod_info in ['16-QAM','32-QAM','64-QAM','128-QAM','256-QAM']:
            # ref.: [3] Eqs. 11-12, valid for SNR values between -5 dB and 20 dB
            
            # TODO: fix large estimation deviation for 16, 32, 64 QAM (128 and 256 QAM are good)
            # dictionary with coefficient values for every available QAM order (ref.: [3] Eq. 11)
            coeff_dict = {'qam_order' : np.asarray([16,32,64,128,256]),
                          'coeff' : np.asarray([[0.36060357620798,-1.28034019700542,1.81839487545758,-1.29109105371317,0.45823466671427,-0.06503489292716],
                                                [0.53715056289170,-1.85210885301961,2.55864235321160,-1.76993734603623,0.61298700208470,-0.08502242157078],
                                                [1.81625572448046,-6.24952901412163,8.60050533607873,-5.91706608901663,2.03511551491328,-0.27993710478023],
                                                [0.64033054858630,-2.17678215614423,2.95932583860006,-2.01114864174439,0.68323069211818,-0.09282225372024],
                                                [0.33595278506244,-1.15419807009244,1.58563212231193,-1.08880229086714,0.37369521988006,-0.05128588224013]])}
            
            # pull coefficient vector that corresponds with signal's QAM order from dictionary and scale according to [Eq. 11]
            coeff_vec = coeff_dict['coeff'][np.where(['16-QAM','32-QAM','64-QAM','128-QAM','256-QAM']==mod_info)][0]
            # scale coeff values
            if mod_info == '256-QAM':
                coeff_vec = coeff_vec*1e7
            else:
                coeff_vec = coeff_vec*1e6
                
            # calc z_hat (ref.: [3] Eq. 12)
            r_kI = np.real(samples)
            r_kQ = np.imag(samples)
            z_hat = (np.mean(r_kI**2)+np.mean(r_kQ**2)) / (np.mean(np.abs(r_kI))**2+np.mean(np.abs(r_kQ))**2)

            # calc SNR estimate
            snr_estimate = 10*np.log10(coeff_vec[5]*(z_hat**5)+coeff_vec[4]*(z_hat**4)+coeff_vec[3]*(z_hat**3)+coeff_vec[2]*(z_hat**2)+coeff_vec[1]*(z_hat)+coeff_vec[0])
            
        # all other modulation formats
        else:
            raise ValueError('Modulation format is not supported!')
            
        return snr_estimate
    
    # init array for snr estimates
    snr_estimate = np.full_like(sig.samples,0,dtype='float')
    
    # loop over signal dimensions
    for dim in range(sig.ndims):
        # split samples into blocks
        if block_size != -1:
            # identify number of whole blocks
            n_blocks = sig.samples[dim].len // block_size
            # identify remainder after splitting operation
            block_rem = sig.samples[dim].len % block_size
            
            # call inner function per block
            for i in range(n_blocks):
                snr_tmp = inner(sig.samples[dim][i*block_size:(i+1)*block_size],sig.modulation_info[dim])
                # write tmp array into return array
                snr_estimate[dim][i*block_size:(i+1)*block_size] = snr_tmp
                
            # call function for remainder block
            if block_rem != 0: 
                snr_tmp = inner(sig.samples[dim][-block_rem:],sig.modulation_info[dim])
                # write remainder SNR estimate into return array
                snr_estimate[dim][-block_rem:] = snr_tmp
        
    return snr_estimate