# -*- coding: utf-8 -*-
"""
Created on Thu Sep 3 15:35:49 2020

@author: UET-Labor

This program shows implementation of real filter in order to filter the incoming signals.

"""

import numpy as np
# import sys
# sys.path.append("..\\comm")
import comm as comm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def filter_arbitrary(samples,FILTER,sample_rate=1.0):
    """
    Filter a given signal with an arbitrary filter defined in a table.
    This filter is implemented in frequency domain (cyclic convolution).
    The filter can be defined as single-sided frequency axis (real filter response is assumed),
    or with double-sided frequency axis (complex filter response is assumed).

    Parameters
    ----------
    samples : 1D numpy array, real or complex
        input signal.
    sample_rate : float
        Sample rate of the signal in Hz. The default is 1.0.
    FILTER : 2D numpy array
        1st col.: frequency in Hz
        2nd col.: magnitude in dB
        3rd col.: phase in deg; optional (zero phase assumed if 3rd col. is omitted)

    Returns
    -------
    samples_out : 1D numpy array, real or complex
        filtered output signal.

    """
    
    # check if FILTER has two or three columns
    
    if isinstance(FILTER,np.ndarray):
    
        if FILTER[0].shape == (3,):
                FILTER = FILTER
        if FILTER[0].shape == (2,):
            #FILTER = np.append(FILTER, np.zeros((FILTER.shape[0], 1), dtype=FILTER.dtype), axis=1)
            FILTER = np.column_stack((FILTER, np.zeros((FILTER.shape[0], 1), float)))  # add the third column with phase in it 
                    
    else:
        raise ValueError('FILTER must be a NX2 or NX3 numpy array')
        
    
        
    ###### Sort the data in which the frequency is in ascending order
    u,udx = np.unique(FILTER[:,0],return_index=True) # unique ascending frequencies
    FILTER_sorted = FILTER[udx,:] # adds sorted and unique 1st row into the array table
    
    f_Hz = FILTER_sorted[:,0] # freq axis from the table
    
    # fill the FILTER-sorted with mag_lin and phase_rad in the 2nd adn 3rd column
    FILTER_sorted[:,1] = 10**(FILTER_sorted[:,1]/20)      # Magnitude dB into linear
    FILTER_sorted[:,2] = np.pi/180 * FILTER_sorted[:,2]   # phase angle from degree to radian
    
    
    if all(f_Hz>=0):# check if all frequencies are positive: condition for real filter
        filter_is_real = True
        FILTER_flip= np.flip(FILTER_sorted,0).copy() # flipping  all the arrays columnwise vertically;
        FILTER_flip[:,0]= -FILTER_flip[:,0]          # set all frequencies to negative and adds to FILTER_flip
        FILTER_flip[:,2]= -FILTER_flip[:,2]           # set all phase to negative and adds to FILTER_flip
        
        H = np.concatenate(((FILTER_flip,FILTER_sorted)))  # H = final freq. dom. filter array; concatenate FILTER_flip and FILTER_sorted, generates table of Nx3 order
        
        
        # if f_Hz[0]==0:
        i,idx = np.unique(H[:,0],return_index=True) # dicard double zero  frequency entries if any
        H = H[idx,:] 
            
    else: # doublesided(complex) filter definition
        filter_is_real = False # check correctly for real filter (complex conj. filter definition)
        H = FILTER_sorted
    
    #phase unwrap 
    H[:,2] = np.unwrap(H[:,2]) 
    
    ##### Interpolator for magnitude and phase 
    # f_mag = interp1d(H[:,0],H[:,1],bounds_error=False,fill_value='extrapolate', kind='linear') # magnitude interpolation
    # f_ph = interp1d(H[:,0],H[:,2],bounds_error=False,fill_value='extrapolate', kind='linear') # phase interpolation
    f_mag = interp1d(H[:,0],H[:,1],bounds_error=False,fill_value=(H[0,1],H[-1,1]), kind='linear') # magnitude interpolation
    f_ph  = interp1d(H[:,0],H[:,2],bounds_error=False,fill_value=(H[0,2],H[-1,2]), kind='linear') # phase interpolation
    
    # Frequency from input signal,FFT frequencies
    f_sig = np.fft.fftshift(np.fft.fftfreq(samples.size,1/sample_rate)) 
    
    #interpolation to input signal frequency axis
    f_mag_ip = f_mag(f_sig)   # with magnitude
    f_ph_ip = f_ph(f_sig) # with phase
    
    # FFT of the input samples
    X_f = np.fft.fft(samples) 
    
    #Calculates Transfer function H(f) and Inverse FFTSHIFT(!!!Important)
    # H_f = (H[:,1])*np.exp(1j*H[:,2])  # H*e^j*phi
    H_interp = np.fft.ifftshift(f_mag_ip * np.exp(1j*f_ph_ip)) 
        
    ##### Multiplying interpolated Transfer function and FFT of the input signal
    Y_f = X_f * H_interp
    
    ### IFFT : back to time domain
    samples_out = np.fft.ifft(Y_f)
    
    ### check for real output samples
    if all(np.isreal(samples)) and filter_is_real:
        samples_out = np.real(samples_out)
    
    return samples_out
