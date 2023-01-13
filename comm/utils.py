import math

import numpy as np
from scipy import optimize
import scipy.interpolate as sinter
import scipy.special as sspecial
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial
import scipy.signal as ssignal
import scipy.interpolate as sinterp

from . import signal
from . import rx


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


def estimate_osnr_spectrum(power_vector = [], wavelength_vector = [], interpolation_points = [], integration_area = [], resolution_bandwidth = 0.1, polynom_order = 3, plotting = False):

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

def estimate_snr_spectrum(x, y, sig_range, noise_range, order=1, noise_bw=12.5e9, scaling='lin', fit_lin=True, plotting=False):
    """
    Estimate the signal to noise ratio (SNR) from a given power spectrum.
    
    The SNR is estimated from a given power spectrum using the so called interpolation
    method: The signal (and noise (n2)) power (p_sig_n2) is calculated by integration
    of the spectrum between the given spectral signal range. 
    Then, the noise floor is estimated by fitting a polynomial of given order to 
    the spectral points lying in two specified spectral ranges (given by the four
    values in noise_range). In general, one spectral range left and another one right
    from the data signal is specified, which therefore requires some oversampling of 
    the data signal. Two noise powers (p_n1 and p_n2) are calculated from 
    the estimated polynomial by integration within noise_bw (n1) and sig_range (n2), 
    respectively.
    The SNR is then estimated by SNR = (p_sig_n2 - p_n2) / p_n1.
    
    It can be specified if the polynomial is fitted from the spectrum either
    in linear or in logarithmic domain.
    

    Parameters
    ----------
    x : 1D array, float
        "x-axis" values of the spectrum. In general either frequency or wavelength.
        The values in this vector must be monotonically increasing.
    y : 1D array, float
        "y-axis" values of the spectrum. In general either power or power density.
    sig_range : 1D array, float
        Two values (same unit as x) that specify the left and right corner points
        of the data signal, respectively, and thus define the integration limits 
        for p_sig_n2.
    noise_range : 1D array, float
        Four values (same unit as x) indicating the first (noise_range[:2]) and 
        second (noise_range[2:]) spectral ranges used to fit the polynomial of 
        the noise floor. These ranges are generally chosen to lie to the left 
        and right of sig_range, respectively.
    order : int, optional
        Order of the polynomial to be fitted. The default is 1.
    noise_bw : float, optional
        Noise bandwidth (same unit as x) used to calculate the noise power (p_n1) 
        in the SNR formula. For optical SNR, usually set to 0.1 nm or 12.5 GHz. 
        For electrical SNR, usually set to the symmetrical rate of the data signal. 
        The default value is 12.5e9.
    scaling : string, optional
        Is the spectrum (y) given in linear 'lin' (W, V**2, W/Hz, V**2/nm, ...) 
        or in logarithmic 'log' (dBm, dBm/Hz, dBm/m, ...) scale. The default is 'lin'.
    fit_lin : bool, optional
        Is the polynomial fitted from the spectrum linear (True) or logarithmic 
        spectrum (False)? The default is True.
    plotting : boolean, optional
        Shall the spectrum (plus integration regions and noise fit) be plotted?
        Useful for tuning sig_range and noise_range and verifying the fitting results.
        The default is False.

    Returns
    -------
    snr_db : float
        Estimated SNR in dB.

    """
    
    if not (isinstance(x, np.ndarray) and isinstance(y, np.ndarray) and isinstance(noise_range, np.ndarray) and isinstance(sig_range, np.ndarray)):
        raise TypeError('x, y, sig_range and noise range must be numpy arrays')
    
    if (x.ndim != 1) or (y.ndim != 1) or (noise_range.ndim != 1) or (sig_range.ndim != 1):
        raise ValueError('x, y, sig_range and noise range must be 1 dimensional numpy arrays')
    
    if sig_range.size != 2:
        raise ValueError('sig_range must be of size 2')
        
    if noise_range.size != 4:
        raise ValueError('noise_range must be of size 4')
    
    
    if scaling == 'log':
        y = 10**(y/10)
    elif scaling == 'lin':
        pass
    else:
        raise ValueError('parameter scaling must be either "log" or "lin"')
    
    # prepare signal range    
    sig_range.sort()
    sig_range = np.expand_dims(sig_range, axis=-1)
    # find x axis samples closest to given signal range
    sig_range_idx = np.argmin(np.abs(x-sig_range), axis=-1)
    sig_range=np.squeeze(sig_range)

    # prepare noise ranges
    noise_range.sort()
    noise_range = np.expand_dims(noise_range, axis=-1)
    # find x axis samples closest to given noise ranges range
    noise_range_idx = np.argmin(np.abs(x-noise_range), axis=-1)
    noise_range = np.squeeze(noise_range)

    # calc signal and noise power (sig+n2) in signal range numerically from given y
    y_sig_n2 = y[sig_range_idx[0]:sig_range_idx[1]]    
    p_sig_n2 = np.trapz(y_sig_n2, x=x[sig_range_idx[0]:sig_range_idx[1]])
    
    if fit_lin:
        # fit a polynominal with given order from the samples of y which are in noise_range
        # --> supposed to be the noise part of the spectrum
        x_n = np.append(x[noise_range_idx[0]:noise_range_idx[1]], x[noise_range_idx[2]:noise_range_idx[3]])
        y_n = np.append(y[noise_range_idx[0]:noise_range_idx[1]], y[noise_range_idx[2]:noise_range_idx[3]])
        c = Polynomial.fit(x_n, y_n, order)
        
        # find the (coefficients of the) indefinite integral of the fitted polynominal
        C = c.integ()
        
        # calc the power of the fitted polynominal in noise_bw centered around mean signal_range
        # --> integration of fitted polynominal 
        # --> supposed to be the noise power (in noise_bw) under the data signal (n1)
        x_sig_mean = np.mean(sig_range)
        p_n1 = C(x_sig_mean + noise_bw/2) - C(x_sig_mean - noise_bw/2)
        
        # calc the power of the fitted polynominal in signal range
        # --> integration of fitted polynominal 
        # --> supposed to be the noise power which was included in the numverical integration of y (n2)
        p_n2 = C(x[sig_range_idx[1]]) - C(x[sig_range_idx[0]])        
    else:
        # fit a polynominal with given order from the samples of y in logarithmic scale
        # which are in noise_range
        # --> supposed to be the noise part of the spectrum
        x_n = np.append(x[noise_range_idx[0]:noise_range_idx[1]], 
                        x[noise_range_idx[2]:noise_range_idx[3]])
        
        y_n = 10*np. log10(np.append(y[noise_range_idx[0]:noise_range_idx[1]], 
                                     y[noise_range_idx[2]:noise_range_idx[3]]))
        c = Polynomial.fit(x_n, y_n, order)
        
        # calc the power of the fitted polynominal in noise_bw centered around mean signal_range
        # --> numerical integration of fitted polynominal 
        # --> supposed to be the noise power (in noise_bw) under the data signal (n1)
        x_sig_mean = np.mean(sig_range)
        noise_bw_first_idx = np.argmin(np.abs(x - (x_sig_mean - noise_bw/2)), axis=-1)
        noise_bw_last_idx = np.argmin(np.abs(x - (x_sig_mean + noise_bw/2)), axis=-1)
        x_n1 = x[noise_bw_first_idx:noise_bw_last_idx+1]
        y_n1 = 10**(c(x_n1)/10)
        p_n1 = np.trapz(y_n1, x=x_n1)
        
        # calc the power of the fitted polynominal in signal range
        # --> numerical integration of fitted polynominal 
        # --> supposed to be the noise power which was included in the numverical integration of y (n2)
        x_n2 = x[sig_range_idx[0]:sig_range_idx[1]+1]
        y_n2 = 10**(c(x_n2)/10)
        p_n2 = np.trapz(y_n2, x=x_n2)
    
    # calc SNR
    snr = (p_sig_n2 - p_n2) / p_n1
    snr_db = 10 * np.log10(snr)
    
    if plotting:
        # create samples of fitted polynominal within the range of x (for plotting only)
        xx_n, yy_n = c.linspace(n=x.size, domain=[x[0], x[-1]])
        if not fit_lin:
            yy_n = 10**(yy_n/10)
        
        # # linear plot
        # plt.plot(x,y)
        # plt.plot(x[sig_range_idx], y[sig_range_idx], 'o')
        # plt.plot(x[noise_range_idx[0]:noise_range_idx[1]], y[noise_range_idx[0]:noise_range_idx[1]], 'r--')
        # plt.plot(x[noise_range_idx[2]:noise_range_idx[3]], y[noise_range_idx[2]:noise_range_idx[3]], 'r--')
        # plt.plot(xx_n,yy_n,'g-')
        # plt.xlabel('given x value / a.u.')
        # plt.ylabel('given y value / a.u.')
        # plt.title('est. SNR = {:.1f} dB in noise bandwidth of {:.2e}'.format(snr_db, noise_bw))
        # plt.show()
    
        # logarithmic plot
        plt.figure()
        plt.plot(x,10*np.log10(y))
        plt.plot(x[sig_range_idx], 10*np.log10(y[sig_range_idx]), 'o')
        plt.plot(x[noise_range_idx[0]:noise_range_idx[1]], 10*np.log10(y[noise_range_idx[0]:noise_range_idx[1]]), 'ro-')
        plt.plot(x[noise_range_idx[2]:noise_range_idx[3]], 10*np.log10(y[noise_range_idx[2]:noise_range_idx[3]]), 'ro-')
        plt.plot(xx_n,10*np.log10(yy_n),'g-')
        plt.xlabel('given x value / a.u.')
        plt.ylabel('given y value / dB')
        plt.title('est. SNR = {:.1f} dB in noise bandwidth of {:.2e}'.format(snr_db, noise_bw))
        plt.show()
        
    return snr_db

def est_snr_spec_wrapper(sig,roll_off,plotting=False):
    
    if roll_off.size != sig.n_dims:
        raise ValueError('Length of array of roll_off factors and number of signal dimensions be equal.')
    # init snr_dB array
    snr_dB = np.zeros(shape=(sig.n_dims,),dtype='float')
    # loop over signal dimensions
    for dim in range(sig.n_dims):
        #### generate freq axis from signal attributes
        sr = sig.sample_rate[dim]
        n_samples = sig.samples[dim].size
        # frequency delta
        df = sr/n_samples
        faxis = np.linspace((-sr+df)/2,(sr-df)/2,n_samples)
        #### calc linear power spectrum
        pwr_vec = np.abs(np.fft.fftshift(np.fft.fft(sig.samples[dim])))**2
        #### calc bandwidth from symbol rate and roll off factor
        bw_half = sig.symbol_rate[dim]/2 + (1+roll_off[dim])
        sig_range = np.array([-bw_half,bw_half])
        #### calc noise range
        noise_range = np.array([-bw_half-1e9,-bw_half,bw_half,bw_half+1e9])
        #### call estimation function
        snr_dB[dim] = estimate_snr_spectrum(faxis, pwr_vec, sig_range, noise_range,
                                            order=1,noise_bw=sig.symbol_rate[dim],
                                            scaling='lin',plotting=plotting)
        
    return snr_dB
        
    
    
        
def estimate_snr_nda(sig,block_size=-1,bias_comp=True):
    """
    Estimates the SNR per symbol of a noisy signal depending on its modulation 
    format. The noise is assumed to be Gauss-distributed. Different non-data-
    aided algorithms are used in conjunction with empirically determined 
    correction/modification terms depending on whether the modulation format 
    is BPSK, QPSK, or QAM of order 16 or upwards. The function assumes that 
    the 'samples' attribute of the signal class object contains symbols, i.e. 
    that the signal has been downsampled to 1 sample per symbol, and that no 
    sampling phase error is present.

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
        be known, i.e. a string variable must be present in the 'modulation_info'
        attribute per dimension of the signal.
    block_size : int
        The number of symbols to average SNR over. Greater block size means a 
        more accurate estimate. Default value is -1, which treats the entire 
        symbol vector as one block.
    bias_comp: bool
        Flag determining whether the bias of the *BPSK estimation algorithm* is compensated
        according to the data obtained in a Monte Carlo simulation with 10.000 runs for 
        SNR values between 0 dB and 20 dB. Beyond these limits, the bias is extrapolated 
        from the data. Default value is True.

    Returns
    -------
    snr_estimate : numpy array
        Array with the same shape as sig.samples containing the estimated SNR values.

    """

    if type(block_size) != int:
        raise TypeError("Block size must be a positive integer or -1!")
    
    if (block_size < 1) & (block_size != -1):
        raise ValueError("Block size must be a positive integer or -1!")
    
    # nested function that performs the estimation
    def _estimation_helper(samples,mod_info):
        
        # BPSK case
        if mod_info == '2-PSK':
            # ref.: [1], Eqs. 6, 17, valid for SNR values between 0 dB and 20 dB
            # calc 2nd order moment
            m2_hat = np.mean(np.abs(samples)**2)
            # calc novel signal power estimate
            s_hat = ((np.sum(np.abs(np.real(samples)))**2) + (np.sum(np.abs(np.imag(samples)))**2)) / (len(samples)**2)
            # calc SNR estimate
            snr_estimate = 10*np.log10(s_hat / (m2_hat - s_hat))
            
            if bias_comp:
                # estimation bias, obtained through Monte Carlo simulation with 10.000 runs
                # for SNR values between 0 dB and 20 dB - beyond these values, the estimation
                # bias is extrapolated
                bias = np.asarray([4.00640498, 3.46586717, 3.0437816 , 2.69519762, 2.43301419,
                                   2.26609504, 2.1298437 , 2.04872834, 1.99272736, 1.94593644,
                                   1.91146281, 1.89347751, 1.85535157, 1.84488699, 1.8299747 ,
                                   1.81262673, 1.81382646, 1.79700308, 1.79183747, 1.7961362 ,
                                   1.79334496])
                
                # interpolate between LUT of bias values and correct estimation value
                bias_inter = sinter.interp1d(np.arange(0,21),bias,fill_value='extrapolate')
                snr_estimate = snr_estimate - bias_inter(snr_estimate)
            
        # QPSK case
        elif mod_info in ['4-PSK','4-QAM','QPSK']:
            # ref.: [2], Eqs. 11-15, valid for SNR values between -10 dB and 30 dB
            # calc mean
            symb_mean = np.abs(samples).mean()
            # calc variance
            symb_var = np.abs(samples).var()
            # calc SNR estimate
            snr_estimate = 10*np.log10((np.abs(symb_mean)**2)/(2*symb_var))
            # for SNR values below 10 dB, modify the estimate
            if snr_estimate < 10:
                snr_estimate = np.sqrt((snr_estimate-2.5)*39.2)-7
                # exception handling if NaN is returned:
                # assumption: when np.sqrt returns NaN, the SNR is very low - to at least
                # avoid NaNs, replace them with -inf to allow for comparison when combining
                if snr_estimate == np.NaN:
                    snr_estimate = float('-inf')
            
        # QAM case (order between 16 and 256, either square or symmetrical QAM)
        elif mod_info in ['16-QAM','32-QAM','64-QAM','128-QAM','256-QAM']:
            # ref.: [3] Eqs. 11-12, valid for SNR values between -5 dB and 20 dB
            
            # TODO: fix large estimation deviation for 128-QAM exclusively
            if mod_info == '128-QAM':
                raise ValueError("Estimation algorithm for 128-QAM yields unusable results. This particular modulation format is not supported for now.")
            # dictionary with coefficient values for every available QAM order (ref.: [3] Eq. 11)
            coeff_dict = {'qam_order' : np.asarray([16,32,64,128,256]),
                          'coeff' : np.asarray([[0.36060357620798,-1.28034019700542,1.81839487545758,-1.29109195371317,0.45823466671427,-0.06503489292716],
                                                [0.53715056289170,-1.85210885301961,2.55864235321160,-1.76993734603623,0.61298700208470,-0.08502242157078],
                                                [1.81625572448046,-6.24952901412163,8.60050533607873,-5.91706608901663,2.03511551491328,-0.27993710478023],
                                                [0.64033054858630,-2.17678215614423,2.95932583860006,-2.01114864174439,0.68323069211818,-0.09282225372024],
                                                [0.33595278506244,-1.15419807009244,1.58563212231193,-1.08880229086714,0.37369521988006,-0.05128588224013]],dtype='float64')}
            
            # pull coefficient vector that corresponds with signal's QAM order from dictionary and scale according to [Eq. 11]
            coeff_vec = coeff_dict['coeff'][np.argwhere(np.asarray(['16-QAM','32-QAM','64-QAM','128-QAM','256-QAM'])==mod_info)][0][0]
            # scale coeff values
            if mod_info == '256-QAM':
                coeff_vec = coeff_vec*1e7
            else:
                coeff_vec = coeff_vec*1e6
                
            # calc z_hat (ref.: [3] Eq. 12)
            r_kI = np.real(samples)
            r_kQ = np.imag(samples)
            z_hat = (np.mean(r_kI**2)+np.mean(r_kQ**2)) / ((np.mean(np.abs(r_kI))**2)+(np.mean(np.abs(r_kQ))**2))

            # calc SNR estimate
            snr_estimate = coeff_vec[5]*(z_hat**5)+coeff_vec[4]*(z_hat**4)+coeff_vec[3]*(z_hat**3)+coeff_vec[2]*(z_hat**2)+coeff_vec[1]*(z_hat)+coeff_vec[0]
            
        # all other modulation formats
        else:
            raise ValueError('Modulation format is not supported or not specified!')
            
        return snr_estimate
    
    # init array for snr estimates
    # TODO: full_like method fails when signal dimensions have an unequal number
    # of elements! If such a case should occur, preallocation of estimation vec-
    # tor should be performed individually per signal dimension and a separate 
    # object would have to be constructed to output the estimates
    snr_estimate = np.full_like(sig.samples,0,dtype='float')
    
    # loop over signal dimensions
    for dim in range(sig.n_dims):
        # split samples into blocks
        if block_size != -1:
            # identify number of whole blocks
            n_blocks = sig.samples[dim].shape[0] // block_size
            # identify remainder after splitting operation
            block_rem = sig.samples[dim].shape[0] % block_size
            
            # call inner function per block
            for i in range(n_blocks):
                snr_tmp = _estimation_helper(sig.samples[dim][i*block_size:(i+1)*block_size],sig.modulation_info[dim])
                # write tmp array into return array
                snr_estimate[dim][i*block_size:(i+1)*block_size] = snr_tmp
                
            # call function for remainder block
            if block_rem != 0: 
                snr_tmp = _estimation_helper(sig.samples[dim][-block_rem:],sig.modulation_info[dim])
                # write remainder SNR estimate into return array
                snr_estimate[dim][-block_rem:] = snr_tmp
        
        # if block_size = -1, perform estimation on the entire sample vector
        else:
            snr_estimate[dim].fill(_estimation_helper(sig.samples[dim],sig.modulation_info[dim]),dtype='float')
        
    return snr_estimate

def estimate_SNR_m2m4(samples, constellation):
    """
    ref.: Aifen Wang; Hua Xu; Jing Ke, "NDA moment-based SNR estimation for 
    envelope-based QAM", 2012 IEEE 11th International Conference on Signal Processing

    Parameters
    ----------
    samples : TYPE
        DESCRIPTION.
    constellation : TYPE
        DESCRIPTION.

    Returns
    -------
    snr_estimate_m2m4 : TYPE
        DESCRIPTION.

    """
    # normalize samples to have a mean power of 1
    samples_norm = samples / np.sqrt(np.mean(np.abs(samples)**2))
    # normalize constellation to have a mean power of 1 (?)
    const_norm = constellation / np.sqrt(np.mean(np.abs(constellation)**2))
    # find number of unique amplitudes in ideal constellation and how often they appear
    A, A_cnt = np.unique(np.abs(const_norm),return_counts=True)
    # find probability of constellation point with unique amplitude occuring
    p = A_cnt / (const_norm.size)
    # calc constellation moments
    c2 = 1 # due to normalization to mean power of 1
    c4 = np.sum(p*(A**4))
    # c6 = np.sum(p*(A**6))
    # calc sample moments
    M2 = np.mean(np.abs(samples_norm)**2)
    M4 = np.mean(np.abs(samples_norm)**4)
    # calc enumerator and denominator of final SNR estimation eq
    enum = 1 - 2*((M2**2)/M4) - np.sqrt((2-c4)*((2*(M2**4)/(M4**2))-((M2**2)/M4)))
    denom = (c4 * (M2**2) / M4) - 1
    
    snr_estimate_m2m4 = enum/denom
    
    return snr_estimate_m2m4

def estimate_SNR_evm(sig, **kwargs):
    """
    Estimate SNR (in dB) based on the calculated EVM.
    
    The EVM is calculated using the method comm.utils.calc_evm() and from this 
    result, the SNR is derived according to [1].
    
    For more information about the valid parameters see comm.utils.calc_evm().
    
    [1] R. A. Shafik, et al. "On the Extended Relationships Among EVM, BER and 
    SNR as Performance Metrics" https://doi.org/10.1109/ICECE.2006.355657 

    Returns
    -------
    list of floats
        Derived SNR values in dB per signal dimension.

    """
    # calc EVM
    evm = calc_evm(sig, **kwargs)
    # convert to SNR
    return 10*np.log10(1/(np.asarray(evm)**2))


def _evm_helper(scaling_fac, symbols, symbols_ref):
    """
    Helper function for calc_evm to optimize scaling factor.
    
    The symbols are clusterd to the sent reference symbols and the reference symbols
    are subtracted. Then the sum of the absolute values of the individual means in 
    real and imaginary part of all clusters is returned.

    Parameters
    ----------
    scaling_fac : float
        Factor for scaling the symbols.
    symbols : 1D numpy array, real or complex
        Symbols to calc the EVM from.
    symbols_ref : 1D numpy array, real or complex
        Reference (sent) symbols.

    Returns
    -------
    float
        Sum of all means of real and imaginary part of all clusters.

    """
    
    # scale symbols
    symbols_norm = symbols * scaling_fac
    
    means = []
    
    for const_point in np.unique(symbols_ref):
        # cluster symbols to individual constellation points
        symbol_norm = symbols_norm[symbols_ref == const_point]
        symbol_ref = symbols_ref[symbols_ref == const_point]
        # calc mean (of real an imaginary part) of each cluster:
        # the mean will be smallest for the scaling factor which leads to the
        # same means as the original constellations 
        means.append(np.abs(np.mean(symbol_norm - symbol_ref)))
        
    return np.sum(np.asarray(means))   


def calc_evm(sig, norm='max', method='blind', opt=False, dimension=-1):
    """
    Calculate the error vector magnitude (EVM).
    
    The EVM [1] is calculated for given signal considering the received modulation 
    symbols and the ideal constellation points as reference symbols. 
    
    The signal (sig.samples) has to be sampled at one sample per symbol. 
    
    The EVM Normalization Reference [2] can be specivied as constellation maximum
    'max' or as reference RMS 'rms'.
    
    If the input parameter method is given as 'blind' the reference symbols are derived by
    normalization of the samples (sig.samples) to the same power as the ideal 
    constellation points (sig.constellation) and following decision to these 
    ideal constellation points. 
    In case of 'data_aided', the actual sent symbol sequence
    (sig.symbols) is used as reference after normalization. In this case a temporal and
    phase synchronizaiton (comm.rx.symbol_sequence_sync) is performed before
    calculation of the error vector.
    It can further be specified (opt==False), if the normalization of the samples is done
    by scaling the samples to the same mean power as the ideal constellation 
    points (similar to the 'blind' case). In case of opt==True the scaling factor 
    is optimized. Therefore, the samples are 'clustered' to the individual sent
    constellation points and the scaling is adjusted in order to minimize the mean value
    (in real and imaginary part separate) of all clusters.    
     
    NOTE: In case of method=='blind' the error vector is calculated between the 
    received symbols and the DECIDED constellation points and not between the 
    received symbols and the ACTUALLY ("really") sent constellations as in case of
    method=='data_aided'. The former method will therefore lead to an too optimistic 
    EVM in case of low SNR (and many wrong symbol decisions e.g. high BER).  
        
    
    [1] https://rfmw.em.keysight.com/wireless/helpfiles/89600b/webhelp/subsystems/digdemod/Content/digdemod_symtblerrdata_evm.htm
    
    [2] https://rfmw.em.keysight.com/wireless/helpfiles/89600b/webhelp/subsystems/digdemod/Content/dlg_digdemod_comp_evmnormref.htm

    Parameters
    ----------
    sig : comm.signal.Signal
        signal containing the symbols (samples) and the original constellation.    
    norm : string, optional
        Specifies the EVM Normalization Reference [2] and can either be 
        constellation maximum 'max' or reference RMS 'rms'. The default is 'max'.
    method : string, optional
        Specifies if the reference symbols are derived without knowledge of the 
        sent symbols (case 'blind') or not (case 'data_aided'). The default is 
        'blind'.
    opt : bool, optional
        Specifies if the scaling of the symbols (samples) should be optimized.
        The default is False.
    dimension : int, optional
        Which dimension to operate on? -1 operates on all signal dimensions.
        The default is -1.

    Returns
    -------
    evm : list of floats
        calculated EVM value per sigmal dimension as ratio (to convert to 
        percent, the ratio has to be multiplied by 100).

    """
    
    if type(sig) != signal.Signal:
        raise TypeError("input parameter must be of type 'comm.signal.Signal'")
        
    if dimension == -1:
        dims = range(sig.n_dims)
    elif (dimension >= sig.n_dims) or (dimension < -1):
        raise ValueError("-1 <= dimension < sig.n_dims")
    else:
        dims = [dimension]
        
    evm = np.full(len(dims), np.nan)
    
    # iterate over specified signal dimensions
    for dim in dims:    
    
        if norm == 'max':
            evm_norm_ref = np.max(np.abs(sig.constellation[dim]))
        elif norm == 'rms':
            evm_norm_ref = np.sqrt(np.mean(np.abs(sig.constellation[dim])**2))        
        
        # (starting) scaling factor for received constellation
        scaling_fac = np.sqrt(np.mean(np.abs(sig.constellation[dim])**2) / np.mean(np.abs(sig.samples[dim])**2))
        
        # find reference symbols 'blindly', i.e. by decision using tx consellation
        if method == 'blind': 
            # normalize received constellation symbols to match same mean
            # power as the ideal constellation (will be suboptimal for large noise)
            symbols_norm = sig.samples[dim] * scaling_fac
            # decide symbols        
            symbols_ref = rx.decision(symbols_norm, sig.constellation[dim], norm=False)                                    
                
        # use given reference symbols, i.e. sent symbol sequence
        elif method == 'data_aided':            
            # sync samples to sent symbol sequence
            sig_tmp = rx.symbol_sequence_sync(sig)            
            if sig_tmp.symbols[dim].size < sig_tmp.samples[dim].size:
                # repeat sent reference symbols in order to match length of symbols
                ratio_base = sig_tmp.samples[dim].size // sig_tmp.symbols[dim].size
                ratio_rem = sig_tmp.samples[dim].size % sig_tmp.symbols[dim].size
                symbols_ref = np.concatenate((np.tile(sig_tmp.symbols[dim], ratio_base), sig_tmp.symbols[dim][:ratio_rem]), axis=0)
            else:
                symbols_ref = sig_tmp.symbols[dim]           
            
            if opt:                
                # optimize scaling factor
                # the optimal scaling results in the same mean values (in real
                # and imaginary part of each constellation point 'cluster') as the 
                # sent constellation points 
                result = optimize.minimize(_evm_helper, x0=scaling_fac, 
                                            args=(sig.samples[dim], symbols_ref), method='Nelder-Mead')
                scaling_fac_opt = result.x
                
                symbols_norm = sig.samples[dim] * scaling_fac_opt                
            else:                
                # normalize received constellation symbols to match same mean
                # power as the ideal constellation (will be suboptimal for large noise)
                symbols_norm = sig.samples[dim] * scaling_fac                         
        else:
            raise ValueError("unkown method specified, needs to be 'blind' or 'data_aided'")
        
        # calc error
        error = symbols_norm - symbols_ref
        mean_error = np.sqrt(np.mean(np.abs(error)**2))
        # calc evm from (optimized) error        
        evm[dim] = mean_error / evm_norm_ref
        
    return evm


def combine_OSA_traces(x_data, y_data, operator='-', x0=1550e-9):
    """
    Combine multiple spectra of an optical spectrum analyzer.
    
    This function combines multiple spectra which are given as lists of np.arrays
    in x_data and y_data. For each spectrum the method of combination can be 
    specified as addition ('+'), subtraction ('-'), multiplication ('*') or 
    division ('/') using the operator string.
    
    NOTE: This function can only handle spectra with identical x_axis yet. No
    interpolation is performed.

    Parameters
    ----------
    x_data : list od np.ndarray
        Contains the x_axis (either frequency or wavelength) of the spectra to be
        combined.
    y_data : list od np.ndarray
        Contains the y_axis (either dBm W) of the spectra to be combined.
    operator : string
        Each character of the string specifies the combination method for the individual 
        given spectra. Characters can either be '+', '-', '*' or '/'. Note 
        that the number of characters must be by one less than the number of 
        given spectra. The default is '-'.
    x0 : float, optional
        For the specific x0 value the result of the combination is explicitly plotted.
        The default is 1550e-9.    

    Returns
    -------
    comb : np.ndarray
        y values of the combined spectrum.
        
    Examples
    --------
    This simple example calculates the (amplitude) transfer function of an optical device
    by subtraction of two measured spectra (one at the input and another at the
    output of the device). Therefore it assumes that the spectra are given in dBm.
    One specific attenuation at a wavelength of 1550 nm should explicitly be evaluated.
    
    >>> import comm as comm
    >>>
    >>> # get input spectrum to device and save wavelength data to wl1 and power 
    >>> # data to spectrum1 (e.g. by using comm.instrument_control.get_spectrum_HP_71450B_OSA())
    >>>
    >>> # get output spectrum from device and save wavelength data to wl2 and power 
    >>> # data to spectrum2 (e.g. by using comm.instrument_control.get_spectrum_HP_71450B_OSA())
    >>>
    >>> # generate input lists
    >>> x_data = [wl1, wl2]
    >>> y_data = [spectrum1, spectrum2]
    >>> combined = comm.utils.combine_OSA_traces(x_data, y_data, operator='-', x0=1550.15e-9)        

    """
    
    if not (isinstance(x_data, list) and isinstance(y_data, list)):
        raise ValueError('x and y data need to be lists')
    
    if len(x_data) != len(y_data):
        raise ValueError('lengths of x_data and y_data lists must be the same')
    
    if len(x_data) != (len(operator)+1):
        raise ValueError('need one operator less than length of y_data list')
    
    # check if all spectra have the same x axis
    if not (all([all(x_data[0] == elem) for elem in x_data])):
        # TODO: do interpolation to the finest spectrum
        raise ValueError('frequency or wavelength axis of all spetra musst be equal')
    
    if ((x0 < np.min(x_data[0])) or (x0 > np.max(x_data[0]))):
        raise ValueError('x0 needs to be within x axis range')
        
    comb = y_data[0]
    plt.figure(0)
    plt.plot(x_data[0], y_data[0])
    plt.xlabel('wavelength / m or frequency / Hz')
    plt.ylabel('power / dBm or W')
    plt.grid()
        
    for idx, y in enumerate(y_data[1:]):        
        if operator[idx] == '-':
            comb -= y_data[idx+1]
        elif operator[idx] == '+':
            comb += y_data[idx+1]
        elif operator[idx] == '*':
            comb *= y_data[idx+1]
        elif operator[idx] == '/':
            comb /= y_data[idx+1]
        else:
            raise ValueError('operator needs to be "+","-","*" or "/"')
            
        plt.plot(x_data[idx+1], y_data[idx+1])
    
    plt.show()
    att0 = comb[np.argmin(np.abs(x_data[0]-x0))]
        
    plt.figure()
    plt.plot(x_data[0], comb)
    plt.plot(x0, att0, 'ro', label='{:.1f} dB'.format(att0))
    plt.legend()
    plt.grid()
    plt.xlabel('wavelength / m or frequency / Hz')
    plt.ylabel('gain / dB')
    plt.show()
    
    return comb
    
def edfa_model(samples, sample_rate, opt_mid_wl=1550e-9, mode="APC", opt_target=0, opt_noise_figure=4.5, seed=None):
    """
    Tool for simulation of a very simple EDFA model and related noise behavior. 
    Please note: This model is simple and not physically correct in total. The main assumption for this model is that G(ain) >> 1.
    This means that if G == 1, this model will produce physically inaccurate results. 
    The output power generated by the model is given by the output power set with opt_target and additionally the noise power.

    TODO: 
    1. Discuss if this model is enough for related work, in terms of inadequacies in physics.
    2. Do a better validation in a normal python file

    Parameters
    ----------
    samples : np.ndarray
        input signal samples.
    sample_rate : np.ndarray
        sample rate of the input signal samples.
    opt_mid_wl : float, optional
        center wavelengt of the input signal.
        Default is 1550e-9.
    mode : string, optional
        operation mode of the simulated EDFA. Should be "APC" or "AGC", where APC is automatic power control
        and AGC is automatic gain control. 
        Default is "APC".
    opt_target : float, optional
        target output power of the EDFA. If APC modus, this value is given in [dBm], if AGC modus, this value means relative gain in [dB]. 
        Default is 0 [dBm].
    opt_noise_figure : float, optional
        noise figure of the amplifier in [dB]. Specifies the noise behavior.
        Default is 4.5 [dB].
    seed : int, optional
        seed for noise vector generation. 
        Default is "None", for non-fixed seed operation. 

    Returns
    -------
    samples : np.ndarray
        output samples with edfa noise contribution.

    """
    h = 6.62606896e-34 #plancks constant

    sig_pow = np.mean(abs(samples)**2) #calculate signal power

    if mode == "APC" or mode == "power":
        gain = 1e-3*10**(opt_target*0.1)/sig_pow
    elif mode == "AGC" or mode == "gain":
        gain = 10**(opt_target*0.1)
    else:
        raise Exception("EDFA mode not implemented right now. Should be APC/power or AGC/gain")

    #amplify the signal
    samples = samples*np.sqrt(gain)

    #convert wavelength to mid freq
    opt_mid_freq = (299792458/opt_mid_wl)

    #add the noise
    noisePow = abs((10**(opt_noise_figure*0.1)*(gain-1))*h*opt_mid_freq*sample_rate) # formula from VPI Photonic Reference Manual for AmpSysOpt (11) !!!BOTH POLARIZATIONS!!!   
    rng = np.random.default_rng(seed)
    noise = np.sqrt(noisePow/4)*rng.standard_normal((len(samples),2)) # sqrt(noisePow/4) because of noise power calculation for both pols
    samples = samples + (noise.view(dtype=np.complex128)).flatten()

    return samples
