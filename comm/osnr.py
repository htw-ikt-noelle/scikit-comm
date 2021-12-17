import numpy as np              # 1.20.1
import matplotlib as plt        # 3.3.4
from numpy.polynomial import Polynomial
from matplotlib import pyplot as plt

def osnr(power_vector = [], wavelength_vector = [], interpolation_points = [], integration_area = [], resolution_bandwidth = 0.1, polynom_order = 3 ):
    #TODO: Add docstring
    #TODO: Check Inputs
    #TODO: Wrong OSNR results. 0.5db difference.

    """
    osnr

    Function to calculate the OSNR from trace data

    """

    # =============================================================================
    #  Check inputs of correctnes
    # ============================================================================= 


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
    delta_lambda = wavelength_vector[1] - wavelength_vector[0]
    power_vector_lin = 10**np.divide(power_vector,10) * delta_lambda / resolution_bandwidth
    interpol_noise_powers_lin = 10**np.divide(interpol_noise_powers,10) * delta_lambda / resolution_bandwidth

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
    OSNR = 10*np.log10(pseudo_signal_power/pseudo_noise_power)
    bandwidth = closest_integration_wavelength[1] - closest_integration_wavelength[0]
    OSNR_01nm = 10*np.log10(pseudo_signal_power / pseudo_noise_power * 0.1 / bandwidth)

    
    plt.plot(wavelength_vector,interpol_noise_powers_complete,'-',
            wavelength_vector,power_vector,'-',
            integration_area,[power_vector[closest_integration_wavelength_index[0]],power_vector[closest_integration_wavelength_index[1]]],'ro',
            wavelengths_lambda_0_1,power_lambda_0_1,'g.',
            wavelengths_lambda_2_3,power_lambda_2_3,'g.',
            )

    plt.gca().legend(('Optical power from OSA','Integration borders','Area for polyfit','Interpolated noise power' ))
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Power density [dBm/{0}nm]'.format(resolution_bandwidth))        
    plt.ylim(np.min(power_vector)-10,np.max(power_vector)+10)
    plt.grid()
    plt.show()

    return OSNR,OSNR_01nm,10*np.log10(pseudo_signal_noise_power/1e-3),10*np.log10(pseudo_noise_power/1e-3)