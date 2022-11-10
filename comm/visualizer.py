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


def plot_signal(samples, sample_rate=1.0, fNum=1, boundaries=[None, None], 
                tit='time signal'):
    """
    plot singal as a function of time.
    

    Parameters
    ----------
    samples : 1D numpy array, real or complex
        sampled signal.
    sample_rate : float, optional
        The sample rate of the signal. The default is 1.0.
    fNum : int, optional
        Figure number to plot into. The default is 1.
    boundaries : list of int or None, optional
        The boundaries are given as list with two elements (start and end index).
        The signal is only plotted within these given boundaries. A value of None
        specifies the first and last signal sample, respectively. 
        The default is [None, None] and therefore plots the whole signal.
    tit : string, optional
        Title of the plot. The default is 'time signal'.

    Returns
    -------
    None.

    """
    # generate time axis
    t = np.linspace(0, (len(samples)-1)/sample_rate, len(samples))
    
    # cut signal and time axis if necessary
    t = t[boundaries[0]:boundaries[1]]
    samples = samples[boundaries[0]:boundaries[1]]
    
    # plotting
    plt.figure(fNum, facecolor='white', edgecolor='white')
    plt.clf()    
    # if complex input signal -> plot real and imag seperatly
    if np.any(np.iscomplex(samples)):
        plt.subplot(121)
        plt.plot(t, np.real(samples))
        plt.xlabel('time [s]')
        plt.ylabel('amplitude real part')
        plt.grid()
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
        
	
def plot_eye(samples, sample_rate=2, bit_rate=1, fNum=1, 
             boundaries=[None, None], tit='eye diagramm'):
    """
    Plot eye diagram of sampled signal.

    Parameters
    ----------
    samples : 1D numpy array, real or complex
        sampled signal.
    sample_rate : float, optional
        Sample rate of the signal. Please note that the sample_rate
        must be an integer mulitple of the bit_rate.The default is 2.
    bit_rate : float, optional
        Bit rate (or symbol rate) of hte signal. The default is 1.    
    fNum : int, optional
        Figure number to plot into. The default is 1.
    boundaries : list of int or None, optional
        The boundaries are given as list with two elements (start and end index).
        The eye diagram is only plotted within these given boundaries. A value of None
        specifies the first and last signal sample, respectively. 
        The default is [None, None] and therefore the eye diagram contains
        the whole signal.
    tit : string, optional
        Title of the plot. The default is 'eye diagramm'.

    
    Returns
    -------
    None.

    """
         
    sps = sample_rate/bit_rate
    
    # cut signal and time axis if necessary    
    samples = samples[boundaries[0]:boundaries[1]]
            
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
    
    
def plot_poincare_sphere(samplesX, samplesY, decimation=1, fNum = 1, tit = 'Poincaré sphere', labels=True):
    """
    Plot the signal (given as components of the Jones vector) on the Poincaré sphere.
    
    This function converts the given signal, specified by the two components of
    the Jones vector [1] (array samplesX (e_x) and samplesY (e_y)), into the Stokes
    representation [2] and plots it onto the Poincaré sphere [3] in the three 
    dimensional Stokes space (S1, S2, S3). 
    
    Please note that the Stokes parameters S1, S2 and S3 are normalized to the 
    total instantaneous signal power (S0). Therefore, the Poincaré sphere plot
    in this implementation does not reveal any information on degree of 
    polarization of the signal.
    
    [1] https://en.wikipedia.org/wiki/Jones_calculus#Jones_vector
    [2] https://en.wikipedia.org/wiki/Stokes_parameters
    [3] https://en.wikipedia.org/wiki/Polarization_(waves)#Poincar%C3%A9_sphere

    Parameters
    ----------
    samplesX : 1D numpy array, complex
        samples of the input signal, representing the (time dependent) 
        first component of the Jones vector (e_x). Commonly reffered to as the 
        X (or horizontal (H)) polarization component.
    samplesY : 1D numpy array, complex
        samples of the input signal, representing the (time dependent) 
        second component of the Jones vector (e_y). Commonly reffered to as the 
        Y (or vertical (V)) polarization component.
    decimation : int, optional
        take only every decimations-th sample of the input signal. 
        The default is 1.
    fNum : int, optional
        figure number of the plot to be created. The default is 1.        
    tit : string, optional
        title of the plot to be created. The default is 'Poincaré sphere'.
    labels : bool, optional
        Should the Poincaré sphere be plotted with additional labels indicating
        certain, specific polarization states (like H, V, right circular 
        polarized (RCP), etc. (True) or as a plain sphere only showing the three
        coordinate axes (False)? The default is True.

    Returns
    -------
    handles :  dict containing following keys
        fig : matplotlib.figure.Figure
            Figure object the signal is plotted to.
        ax : matplotlib.axes._subplots.Axes3DSubplot
            Axes object which contains the Poincaré sphere artists
        line : mpl_tooklits.mplot3d.art3d.Line3d
            Line object which contains the Stokes parameters (S1, S2, S3).
    """
    
    # helper function to select the artist labeled with 'SOP'
    def _isSOP(line):
        if line.get_label() == 'SOP':
            return True
        else:
            return False
    
    # decimate signal
    samplesX = samplesX[0::decimation]
    samplesY = samplesY[0::decimation]
    
    # calc Stokes parameters
    s0 = np.abs(samplesX)**2 + np.abs(samplesY)**2
    s1 = np.abs(samplesX)**2 - np.abs(samplesY)**2
    s2 = 2 * (samplesX * np.conj(samplesY)).real
    s3 = -2 * (samplesX * np.conj(samplesY)).imag  
    
    # if figure does not exist: plot all artists (Axis, Sphere, labels, etc.)...
    if fNum not in plt.get_fignums():
    
        # prepare figure
        f = plt.figure(fNum) 
        ax = plt.axes(projection ='3d')    
        plt.axis('Off')
        ax.set_box_aspect([1,1,1])
        plt.title(tit)
        
        # prepare sphere coordinates    
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:20j]
        x = np.cos(u) * np.sin(v)
        y = np.sin(u) * np.sin(v)
        z = np.cos(v)
        
        if labels:
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
        else:
            # plot sphere
            ax.view_init(elev=25, azim=45)        
            ax.plot_wireframe(x, y, z, alpha=0.1, color='k')
            # plot axis
            len = 1.8
            ax.quiver(0, 0,0, 0, 0, len, color='k')
            ax.quiver(0, 0, 0, len, 0, 0, color='k')
            ax.quiver(0, 0, 0, 0, len, 0, color='k')
            ax.text(len*1.2, 0, 0, 'S1', size=15)
            ax.text(0, len*1.2, 0, 'S2', size=15)
            ax.text(0, 0, len, 'S3', size=15)           
        
        # plot data and label line artist as 'SOP'
        li = ax.plot3D(s1/s0, s2/s0, s3/s0, '.b', label='SOP')   
        plt.show()
    # ...otherwise only update data in figure (for speed)
    else:
        # get all required handles
        f = plt.figure(fNum)
        ax = f.axes[0]
        # find all artists labeled as 'SOP'
        li = ax.findobj(_isSOP)
        # update data
        li[0].set_xdata(s1/s0)
        li[0].set_ydata(s2/s0)
        li[0].set_3d_properties(s3/s0)
        # update plot
        f.canvas.draw()
        # wait for figure to update
        plt.pause(0.1)        
    
    # prepare return params
    handles = dict()
    handles['fig'] = f
    handles['ax'] = ax
    handles['line'] = li[0]
    return handles