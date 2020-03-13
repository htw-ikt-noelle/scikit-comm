# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 19:01:26 2014

@author: noelle
"""


import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt

    
def plot_spectrum(samples, sample_rate=1, fNum=1, scale='logNorm', tit='spectrum'):
    """ plots the amplitude spectrum of given samples
    
    parameters:
        sample_rate: sample rate of incoming samples
        fNum: figure number to plot into
        scale: plot the spectrum in liear or logarithmic scale, use normalization to max or not 'logNorm' | 'log' | 'linNorm' | 'lin'
        tit: title of figure
    """

    # fft
    fSamples = fft.fft(samples)
    fSamples = fSamples / len(fSamples)
    
    # if real input signal -> mulitply all frequencies, but the DC by 2
    if np.all(np.isreal(samples)):
        fSamples[1:] *= 2        
    
    # calc amplitude and frequency axis
    fSamples = np.abs(fSamples)
    freq = fft.fftfreq(len(fSamples), 1/sample_rate)
    
    # scale spectrum
    if scale == 'logNorm':
        fSamples = 20*np.log10(fSamples / np.max(fSamples))  
        ylabel = "normalized amplitude [dB]"
    elif scale == 'log':
        fSamples = 20*np.log10(fSamples)
        ylabel = "amplitude [dB]"
    elif scale == 'linNorm':
        fSamples = fSamples / np.max(fSamples)
        ylabel = "normalized amplitude [a.u.]"
    elif scale == 'lin':
        fSamples = fSamples
        ylabel = "amplitude [a.u.]"
    else:
        print('plotSpectrum scale must be lin(Norm) or log(Norm)...using "logNorm"')        
        fSamples = 20*np.log10(fSamples / np.max(fSamples))
        ylabel = "normalized amplitude [dB]"    
    
    # plot spectrum
    plt.figure(fNum, facecolor='white', edgecolor='white')
    plt.clf()
    # if signal real -> plot only positive frequencies
    if np.all(np.isreal(samples)):
        plt.plot(freq[0:int(len(fSamples)/2)], fSamples[0:int(len(fSamples)/2)])
    # if signal complex -> plot neg. and pos. frequencies
    else:
        plt.plot(fft.fftshift(freq), fft.fftshift(fSamples))
    plt.title(tit)
    plt.xlabel('frequency [Hz]')
    plt.ylabel(ylabel)
    plt.grid()
    plt.show()


def plot_signal(samples, sample_rate = 1.0, fNum = 1, tit = 'time signal'):
    """ plots the given signals as a funtion of time
    
    parameters:
        sample_rate: sample rate of incoming samples
        fNum: figure number to plot into        
        tit: title of figure
    """
    t = np.linspace(0, (len(samples)-1)/sample_rate, len(samples))
    plt.figure(fNum, facecolor='white', edgecolor='white')
    plt.clf()    
    # if complex input signal -> plot real and imag seperatly
    if np.any(np.iscomplex(samples)):
        plt.subplot(121)
        plt.plot(t, np.real(samples))
        plt.xlabel('time [s]')
        plt.ylabel('amplitude real part')
        plt.subplot(122)
        plt.plot(t, np.imag(samples))
        plt.xlabel('time [s]')
        plt.ylabel('amplitude imaginary part')
        plt.grid()
        plt.show()
    else:
        plt.plot(t, samples)
        plt.xlabel('time [s]')
        plt.ylabel('amplitude')
        plt.grid()
        plt.show()
        
	
def plot_eye(samples, sample_rate = 1, bitRate = 0.5, offset = 0, fNum = 1, tit = 'eye diagramm'):
    """ plots the eye diagramm of a given signals
    
    parameters:
        
    """
         
    sps = sample_rate/bitRate
    t = np.linspace(0, (2 * sps -1) * (1/sample_rate), 2 * sps)
        
    if np.mod(sps, 1):
        raise ValueError('sample_rate must be an integer multiple of bitRate...')
    if np.mod(len(samples), 2*sps):
        raise ValueError('signal must contain an even integer multiple of sps...')
        
    samples = np.reshape(samples, (int(2 * sps), -1), order = 'F')
        
    if np.any(np.iscomplex(samples)):
        plt.subplot(121)
        plt.plot(t, np.real(samples), color = '#1f77b4')
        plt.xlabel('time [s]')
        plt.ylabel('amplitude real part')
        plt.subplot(122)
        plt.plot(t, np.imag(samples), color = '#1f77b4')
        plt.xlabel('time [s]')
        plt.ylabel('amplitude imaginary part')
        plt.grid()
        plt.show()
    else:
        plt.plot(t, samples, color = '#1f77b4')
        plt.xlabel('time [s]')
        plt.ylabel('amplitude')
        plt.grid()
        plt.show()
    
    
	
def plot_hist(samples, nBins=100):
    #TODO: implement automated histogramm
	# check for complex input??
	a = 1
    

def plot_constellation(samples):
    #TODO: implement plot of complex plane
	a = 1