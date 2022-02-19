import numpy as np              # 1.20.1
import matplotlib as plt        # 3.3.4
from numpy.polynomial import Polynomial
from matplotlib import pyplot as plt

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
            Must be same length as wavelength_vector.

        wavelength_vector: numpy array
            Vector with the wavelength values of the OSA trace.
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
        OSNR_val: 
            The calculated OSNR of the integration area.

        OSNR_01nm:
            The calculated OSNR normalized to a noise bandwidth of 0.1nm.

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
        >>> [OSNR,ONSR_1nm] = comm.osnr.osnr(power_vector = power,
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

    except Exception as e:
        print(e)
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
        #closest_interpolation_wavelength_index.append(difference_vector.argmin())
        closest_interpolation_wavelength_index = np.int16(np.append(closest_interpolation_wavelength_index,difference_vector.argmin()))
        #closest_interpolation_wavelength.append(wavelength_vector(closest_interpolation_wavelength_index[idx]))
        closest_interpolation_wavelength = np.append(closest_interpolation_wavelength,wavelength_vector[int(closest_interpolation_wavelength_index[idx])])

    # Correct the input integration area to the nearest wavelength in the wavelength vector
    # - Find the position of these values
    closest_integration_wavelength = np.array([])
    closest_integration_wavelength_index = np.array([])
    for idx,integration_point in enumerate(integration_area):
        difference_vector = np.abs(wavelength_vector - integration_point)
        #closest_integration_wavelength_index.append(difference_vector.argmin())
        closest_integration_wavelength_index = np.int16(np.append(closest_integration_wavelength_index,difference_vector.argmin()))
        #closest_integration_wavelength.append(wavelength_vector(closest_integration_wavelength_index[idx]))   
        closest_integration_wavelength = np.append(closest_integration_wavelength,wavelength_vector[int(closest_integration_wavelength_index[idx])])

    # Getting the wavelengths between lamda 0 and lambda 1
    wavelengths_lambda_0_1 = wavelength_vector[closest_interpolation_wavelength_index[0]:closest_interpolation_wavelength_index[1]+1]

    # Getting the wavelengths between lamda 2 and lambda 3
    wavelengths_lambda_2_3 = wavelength_vector[closest_interpolation_wavelength_index[2]:closest_interpolation_wavelength_index[3]+1]

    # Combine the both wavelengths vectors into one vector.
    sample_point_wavelengths_vector = np.append(wavelengths_lambda_0_1,wavelengths_lambda_2_3)

    # Getting the power between lamda 0 and lambda 1
    power_lambda_0_1 = power_vector[closest_interpolation_wavelength_index[0]:closest_interpolation_wavelength_index[1]+1]

    # Getting the power between lamda 2 and lambda 3
    power_lambda_2_3 = power_vector[closest_interpolation_wavelength_index[2]:closest_interpolation_wavelength_index[3]+1]

    # Combine the both power vectors into one vector.
    sample_point_power_vector = np.append(power_lambda_0_1,power_lambda_2_3)

    # Creation of the polynom
    # - The Polynomial.fit method will give back an scaled version of the the coefficients back. 
    # - To get unscaled values, the convert() method is needed.
    polynom_coeffs = Polynomial.fit(sample_point_wavelengths_vector,sample_point_power_vector,polynom_order).convert().coef

    # Calculate the interpolated noise power values:
    poly = Polynomial(polynom_coeffs)
    interpol_noise_powers_complete = poly(wavelength_vector)


    # For calculation needed span
    interpol_noise_powers = interpol_noise_powers_complete[closest_integration_wavelength_index[0]:closest_integration_wavelength_index[1]+1]

    # To calculate the power the power values must be converted from db to linear.
    delta_lambda = np.diff(wavelength_vector)
    delta_lambda = np.append(delta_lambda,delta_lambda[-1])
    power_vector_lin = 10**np.divide(power_vector,10) * delta_lambda / resolution_bandwidth
    interpol_noise_powers_lin = 10**np.divide(interpol_noise_powers,10) * delta_lambda[closest_integration_wavelength_index[0]:closest_integration_wavelength_index[1]+1] / resolution_bandwidth
    
    # delta_lambda = wavelength_vector[1] - wavelength_vector[0]
    # power_vector_lin = 10**np.divide(power_vector,10) * delta_lambda / resolution_bandwidth
    # interpol_noise_powers_lin = 10**np.divide(interpol_noise_powers,10) * delta_lambda / resolution_bandwidth

    # Calculation noise power
    # pseudo_noise_power = np.sum(interpol_noise_powers_lin)
    pseudo_noise_power = np.trapz(interpol_noise_powers_lin,wavelength_vector[closest_integration_wavelength_index[0]:closest_integration_wavelength_index[1]+1])

    # Calculation signal plus noise power
    # pseudo_signal_noise_power = np.sum(power_vector_lin[closest_integration_wavelength_index[0]:closest_integration_wavelength_index[1]+1]) 
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

    return OSNR_val,OSNR_01nm