import numpy as np
import scipy.signal as signal
import scipy.special as sspecial
import matplotlib.pyplot as plt
import math
import pickle
import json
import time
import csv
import comm

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


with open('C:/Users/noelle/Desktop/tmp/HP1.pickle', 'rb') as f:
    hp_data1 = pickle.load(f)
hp_data1 = hp_data1['A']

with open('C:/Users/noelle/Desktop/tmp/HP2.pickle', 'rb') as f:
    hp_data2 = pickle.load(f)
hp_data2 = hp_data2['A']

with open('C:/Users/noelle/Desktop/tmp/ID1.pickle', 'rb') as f:
    id_data1 = pickle.load(f)
    
with open('C:/Users/noelle/Desktop/tmp/ID2.pickle', 'rb') as f:
    id_data2 = pickle.load(f)


# x_data = [hp_data1['WL_Vector'], hp_data2['WL_Vector'], hp_data2['WL_Vector']]
# y_data = [hp_data1['Trace_data'], hp_data2['Trace_data'], hp_data2['Trace_data']]
x_data = [hp_data1['WL_Vector'], hp_data2['WL_Vector']]
y_data = [hp_data1['Trace_data'], hp_data2['Trace_data']]

comb1 = combine_OSA_traces(x_data, y_data, operator='-', x0=1550.15e-9)


# x_data = [id_data1['WL_vector_m'], id_data2['WL_vector_m']]
# y_data = [id_data1['Trace_data'], id_data2['Trace_data']]
# # y_data = [10**(hp_data1['Trace_data']/10), 10**(hp_data2['Trace_data']/10)]

# comb2 = combine_OSA_traces(x_data, y_data, operator='-', x0=1550.15e-9)

# plt.plot(hp_data1['WL_Vector'], comb1,x_data[1], comb2)






# f_name = 'C:/Users/noelle/Desktop/tmp/test.csv'
# hdr_content = {'Date':time.ctime(),
#                'resolution BW [m]':hp_data['Resolution_BW'],
#                'wavelength [nm]':'power [dBm]'
#                }



# with open(f_name, 'w', newline='') as f:
#     writer = csv.writer(f, delimiter=';')
#     writer.writerow([f_name])
#     for key, value in hdr_content.items():
#         writer.writerow([key, value])
#     # writer.writerows(hp_data['Trace_data'].tolist())
#     writer.writerows(zip(hp_data['WL_Vector'].tolist(), hp_data['Trace_data'].tolist()))

#     # writer.writeheader()
#     # writer.writerow({'first_name': 'Baked', 'last_name': 'Beans'})
#     # writer.writerow({'first_name': 'Lovely', 'last_name': 'Spam'})
#     # writer.writerow({'first_name': 'Wonderful', 'last_name': 'Spam'})