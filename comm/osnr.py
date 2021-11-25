import numpy as np              # 1.20.1
import matplotlib as plt        # 3.3.4
from numpy.polynomial import Polynomial

def osnr(power_vector = [], wavelength_vector = [], interpolation_points = [], integration_area = [], resolution_bandwidth = 0.1, polynom_order = 3 ):


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
        closest_interpolation_wavelength_index.append(difference_vector.argmin())
        closest_interpolation_wavelength.append(wavelength_vector(closest_interpolation_wavelength_index[idx]))

    # Correct the input integration area to the nearest wavelength in the wavelength vector
    # - Find the position of these values
    closest_integration_wavelength = np.array([])
    closest_integration_wavelength_index = np.array([])
    for idx,integration_point in enumerate(integration_area):
        difference_vector = np.abs(wavelength_vector - integration_point)
        closest_integration_wavelength_index.append(difference_vector.argmin())
        closest_integration_wavelength.append(wavelength_vector(closest_integration_wavelength_index[idx]))   

    # Getting the wavelengths between lamda 0 and lambda 1
    wavelengths_lambda_0_1 = wavelength_vector[closest_interpolation_wavelength_index[0]:closest_interpolation_wavelength_index[1]+1]

    # Getting the wavelengths between lamda 2 and lambda 3
    wavelengths_lambda_2_3 = wavelength_vector[closest_interpolation_wavelength_index[2]:closest_interpolation_wavelength_index[3]+1]

    # Combine the both wavelengths vectors into one vector.
    sample_point_wavelengths_vector = wavelengths_lambda_0_1.append(wavelengths_lambda_2_3)

    # Getting the power between lamda 0 and lambda 1
    wavelengths_lambda_0_1 = wavelength_vector[closest_interpolation_wavelength_index[0]:closest_interpolation_wavelength_index[1]+1]

    # Getting the power between lamda 2 and lambda 3
    wavelengths_lambda_2_3 = wavelength_vector[closest_interpolation_wavelength_index[2]:closest_interpolation_wavelength_index[3]+1]

    # Combine the both power vectors into one vector.
    sample_point_power_vector = wavelengths_lambda_0_1.append(wavelengths_lambda_2_3)

    # Creation of the polynom
    # - The Polynomial.fit method will give back an scaled version of the the coefficients back. 
    # - To get unscaled values, the convert() method is needed.
    polynom_coeffs = Polynomial.fit(sample_point_wavelengths_vector,sample_point_power_vector,polynom_order).convert().coef

    # Calculate the interpolated noise power values:
    poly = Polynomial(polynom_coeffs)
    interpol_noise_powers = poly(wavelength_vector[closest_integration_wavelength_index[0]+1:closest_integration_wavelength_index[0]])
 
    # To calculate the power the power values must be converted from db to linear.
    power_vector_lin = 10**np.divide(power_vector,10)
    interpol_noise_powers_lin = 10**np.divide(interpol_noise_powers,10)

    # Calculation noise power
    # - The power at the 3. and 4. interpolation point will be included
    pseudo_noise_power = np.sum(interpol_noise_powers_lin) + power_vector_lin[closest_integration_wavelength_index[0]] + power_vector[closest_integration_wavelength_index[1]]

    # Calculation signal plus noise power
    pseudo_signal_noise_power = np.sum(power_vector_lin[closest_integration_wavelength_index[0]:closest_integration_wavelength_index[1]+1]) 

    # Calculation signalpower
    pseudo_signal_power = pseudo_signal_noise_power - pseudo_noise_power

    # Calculation OSNR
    OSNR = 10*np.log10(pseudo_signal_power/pseudo_noise_power)
    bandwidth = closest_integration_wavelength[1] - closest_integration_wavelength[0]
    OSNR_01nm = 10*np.log10(pseudo_noise_power / pseudo_noise_power * 0.1 / bandwidth)

    return OSNR,OSNR_01nm