
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
import copy

    
def plot_spectrum(samples, sample_rate=1, fNum=1, scale='logNorm', tit='spectrum'):
    """ plots the amplitude spectrum of given samples
    
    parameters:
        sample_rate: sample rate of incoming samples
        fNum: figure number to plot into
        scale: plot the spectrum in liear or logarithmic scale, use normalization to max or not 'logNorm' | 'log' | 'linNorm' | 'lin'
        tit: title of figure
    """

    isReal = np.all(np.isreal(samples))

    # fft
    fSamples = fft.fft(samples)
    fSamples = fSamples / len(fSamples)
    
    # if real input signal -> mulitply all frequencies, but the DC by 2
    if isReal:
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
    if isReal:
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
        plt.title(tit)
        plt.grid()
        plt.show()
    else:
        plt.plot(t, samples)
        plt.xlabel('time [s]')
        plt.ylabel('amplitude')
        plt.title(tit)
        plt.grid()
        plt.show()
        
	
def plot_eye(samples, sample_rate = 2, bit_rate = 1, offset = 0, fNum = 1, tit = 'eye diagramm'):
    """ plots the eye diagramm of a given signals
    
    parameters:
        
    """
         
    sps = sample_rate/bit_rate
            
    if np.mod(sps, 1):
        raise ValueError('sample_rate must be an integer multiple of bit_rate...')
    if np.mod(len(samples), 2*sps):
        raise ValueError('signal must contain an even integer multiple of sps...')
        
    t = np.linspace(0, (2 * sps -1) * (1/sample_rate), int(2 * sps))
    samples = np.reshape(samples, (int(2 * sps), -1), order = 'F')
    
    plt.figure(fNum, facecolor='white', edgecolor='white')
    
        
    if np.any(np.iscomplex(samples)):
        plt.subplot(121)
        plt.plot(t, np.real(samples), color = '#1f77b4')
        plt.xlabel('time [s]')
        plt.ylabel('amplitude real part')
        plt.grid()        
        plt.title(tit)
        plt.subplot(122)
        plt.plot(t, np.imag(samples), color = '#1f77b4')
        plt.xlabel('time [s]')
        plt.ylabel('amplitude imaginary part')
        plt.grid()
        plt.gcf().tight_layout()
        plt.show()
    else:
        plt.plot(t, samples, color = '#1f77b4')
        plt.xlabel('time [s]')
        plt.ylabel('amplitude')
        plt.title(tit)
        plt.grid()
        plt.show()
    
    
	
def plot_hist(samples, nBins=100):
    #TODO: implement automated histogramm
    # check for complex input??
    pass

def plot_constellation(samples, decimation=1, fNum = 1, tit = 'constellation', hist=False, axMax=None, nBins=128):
    """
    Plot the constellation diagramm (complex plane) of samples.

    Parameters
    ----------
    samples : 1D numpy array, complex
        samples of the input signal.
    decimation : int, optional
        take only every decimations-th sample of the input signal. The default is 1.
    fNum : int, optional
        figure number of the plot to be created. The default is 1.
    tit : string, optional
        title of the plot to be created. The default is 'constellation'.
    hist : bool, optional
        should the constellation diagramm be plotted as 2D histogramm? The default is False.
    axMax : float or None
        maximum abolute axis amplitude (equal in x and y axis)
        if None: 1.1 times maximum absolute value of samples (real and imaginalry part) is used
    nBins: int
        number of bins (in each quadrature) for plotting the 2D histogramm
    """
    
    samples = samples[0::decimation]
    
    if axMax is None:
        axMax = max(np.abs(samples.real).max(), np.abs(samples.imag).max())*1.1
    
    plt.figure(fNum, facecolor='white', edgecolor='white')
    
    if hist:
        bins = nBins
        cm = copy.copy(plt.get_cmap("jet"))
        plt.hist2d(samples.real, samples.imag, bins=bins, cmap=cm, cmin=1, density=False)             
    else:     
        plt.plot(samples.real, samples.imag, 'C0.')      
        plt.grid()
    plt.gca().axis('equal')       
    plt.xlim((-axMax, axMax))
    plt.ylim((-axMax,axMax))  
    plt.title(tit)
    plt.xlabel('real part')
    plt.ylabel('imaginary part')    
    plt.show()
    
    
def plot_poincare_sphere(samplesX, samplesY, decimation=1, fNum = 1, tit = 'Poincaré sphere', simple_plot=False):
    
    samplesX = samplesX[0::decimation]
    samplesY = samplesY[0::decimation]
    
    # calc Stokes parameters
    s0 = np.abs(samplesX)**2 + np.abs(samplesY)**2
    s1 = np.abs(samplesX)**2 - np.abs(samplesY)**2
    s2 = 2 * (samplesX * np.conj(samplesY)).real
    s3 = -2 * (samplesX * np.conj(samplesY)).imag
    
    # prepare figure
    plt.figure(fNum) 
    ax = plt.axes(projection ='3d')    
    plt.axis('Off')
    ax.set_box_aspect([1,1,1])
    plt.title(tit)
    
    # prepare sphere coordinates    
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:20j]
    x = np.cos(u) * np.sin(v)
    y = np.sin(u) * np.sin(v)
    z = np.cos(v)
    
    if simple_plot:
        # plot sphere
        ax.view_init(elev=25, azim=45)        
        ax.plot_wireframe(x, y, z, alpha=0.1, color='k')
        # plot axis
        len = 1.8
        ax.quiver(0, 0, 0, 0, 0, len, color='k')
        ax.quiver(0, 0, 0, len, 0, 0, color='k')
        ax.quiver(0, 0, 0, 0, len, 0, color='k')
        ax.text(len*1.2, 0, 0, 'S1', size=15)
        ax.text(0, len*1.2, 0, 'S2', size=15)
        ax.text(0, 0, len, 'S3', size=15)        
    else:
        # plot sphere
        ax.plot_surface(x, y, z, alpha=0.2)
        ax.view_init(elev=15, azim=-65)          
        ax.set_xlabel('S1')
        ax.set_ylabel('S2')
        ax.set_zlabel('S3')
        # plot three rings
        ph = np.linspace(0, 2*np.pi, 20)
        ax.plot3D(np.zeros_like(ph), np.sin(ph), np.cos(ph), 'k')
        ax.plot3D(np.sin(ph), np.zeros_like(ph), np.cos(ph), 'k')
        ax.plot3D(np.sin(ph), np.cos(ph), np.zeros_like(ph), 'k')
        # plot six points (V, H, 45, -45, RCP, LCP)
        ms = 5
        ax.plot3D(0, 0, 1, 'ko', markersize=ms)
        ax.text(0,0,1.2, 'RCP')
        ax.plot3D(0, 0, -1, 'ko', markersize=ms)
        ax.text(0,0,-1.2, 'LCP')
        ax.plot3D(1, 0, 0, 'ko', markersize=ms)
        ax.text(1.2, 0, 0, 'H')
        ax.plot3D(-1, 0, 0, 'ko', markersize=ms)
        ax.text(-1.2, 0, 0, 'V')
        ax.plot3D(0, 1, 0, 'ko', markersize=ms)
        ax.text(0, 1.2, 0, '45')
        ax.plot3D(0, -1, 0, 'ko', markersize=ms)
        ax.text(0, -1.3, 0, '-45')        
    
    # plot data
    ax.plot3D(s1/s0, s2/s0, s3/s0, '.b')   
    plt.show()