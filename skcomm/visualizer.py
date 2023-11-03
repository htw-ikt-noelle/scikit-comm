import copy
import tkinter as tk

import numpy as np
import scipy.signal as ssignal
import numpy.fft as fft
import matplotlib.pyplot as plt
import screeninfo

from . import utils

    
def plot_spectrum(samples, sample_rate=1.0, fNum=None, scale='logNorm',
                  resolution_bw=None, ax_lims=[None,None,None,None], tit='spectrum',
                  save_fig=False, ffolder='.', ffname=None, fformat='png',
                  add_timestamp=False):
    """
    Plot the power spectrum of given samples.    
    
    If the samples are real, a onesided spectrum (only positive frequencies) is
    plotted, otherwise a two-sided spectrum (positive and negative frequencies).
    
    Parameters
    ----------
    samples : 1D np.array, float
        time domain samples of the signal.
    sample_rate : float, optional
        sample rate of the signal in Hz. The default is 1.0
    fNum : int, optional
        figure number to be used for plot. The default is None which uses the 
        "next unused figure number".
    scale : string, optional
        scaling of the y axis, can either be 'logNorm', 'log', 'linNorm','lin'.
        The y axis will be shown in linear or logarithmic scale and can either be
        normalized to the maximum y value or not (absolute values). 
        The default is 'logNorm'.
    resolution_bw : float, optional
        Resolution bandwidth of the displayed spectrum given in 'Hz'. This parameter is
        used for the estimation of the power spectrum using Welch's method. If None, no
        averaging is performed (i.e. the intrinsic resolution, defined by the signals
        temporal length is displayed). The defauls is None.
    ax_lims : list of floats, optional
        Specifies the axis limits of the plot as [xmin, xmax, ymin, ymax]. A value of
        None sets no axis limits. The default is [None, None, None, None].
    tit : string, optional
        title of the plot. The default is 'spectrum'.
    save_fig : bool, optional
        should the plot be saved to file? The default is False.
    ffolder : sting, optional
        folder to save figure to. The default is '.'.
    ffname : string, optional
        filename to save figure to. The default is None, which uses the title
        of the plot as filename.
    fformat : string, optional
        format of the saved file, can either be 'png', 'pdf' or 'svg'. 
        The default is 'png'.
    add_timestamp : bool, optional
        should a timestamp be added to the filename? The default is False.   

    Returns
    ---------- 
    results : dict containing following keys
        'power_spectrum' : np.ndarray
            power spectrum of samples in scaling requested by parameter 'scale'
        'freq' : np.ndarray
            frequency bins of the returned power spectrum in Hz
    """
    
    isReal = np.all(np.isreal(samples))
    
    if resolution_bw is None:
        if isReal:
            freq, power_spectrum = ssignal.welch(samples, sample_rate, window='boxcar', nperseg=samples.size, return_onesided=True, scaling='spectrum', detrend=False)
        else:
            freq, power_spectrum = ssignal.welch(samples, sample_rate, window='boxcar', nperseg=samples.size, return_onesided=False, scaling='spectrum', detrend=False)
    else:
        nperseg = int(sample_rate/resolution_bw)        
        if isReal:
            freq, power_spectrum = ssignal.welch(samples, sample_rate, window='boxcar', nperseg=nperseg, return_onesided=True, scaling='spectrum', detrend=False)
        else:
            freq, power_spectrum = ssignal.welch(samples, sample_rate, window='boxcar', nperseg=nperseg, return_onesided=False, scaling='spectrum', detrend=False)
        
    # scale spectrum
    if scale == 'logNorm':
        with np.errstate(divide='ignore'):
            power_spectrum = 10*np.log10(power_spectrum / np.max(power_spectrum))            
        ylabel = "normalized power [dB]"
    elif scale == 'log':
        with np.errstate(divide='ignore'):
            power_spectrum = 10*np.log10(power_spectrum)            
        ylabel = "power [dB]"
    elif scale == 'linNorm':
        power_spectrum = power_spectrum / np.max(power_spectrum)        
        ylabel = "normalized power [a.u.]"
    elif scale == 'lin':
        power_spectrum = power_spectrum       
        ylabel = "power [a.u.]"
    else:
        print('plotSpectrum scale must be lin(Norm) or log(Norm)...using "logNorm"')
        with np.errstate(divide='ignore'):
            power_spectrum = 10*np.log10(power_spectrum / np.max(power_spectrum))            
        ylabel = "normalized amplitude [dB]"    
    
    # plot spectrum
    if fNum:
        fig = plt.figure(fNum, facecolor='white', edgecolor='white')
    else:
        fig = plt.figure(facecolor='white', edgecolor='white')
        
    plt.clf()
    
    if isReal:
        plt.plot(freq, power_spectrum)   
    else:
        plt.plot(fft.fftshift(freq), fft.fftshift(power_spectrum))
    
    plt.title(tit)
    plt.xlabel('frequency [Hz]')
    plt.ylabel(ylabel)
    plt.gca().set(xlim=(ax_lims[:2]), ylim=(ax_lims[2:]))
    plt.grid(visible=True)
    
    if save_fig:
        if not ffname:
            ffname = tit
        utils.save_fig(fig, fformat=fformat, folder=ffolder, f_name=ffname, 
                 add_timestamp=add_timestamp)
    plt.show()

    results = dict()
    results['power_spectrum'] = power_spectrum
    results['freq'] = freq

    return results


def plot_signal(samples, sample_rate=1.0, fNum=None, boundaries=[None, None], 
                tit='time signal', save_fig=False, ffolder='.', ffname=None, 
                fformat='png', add_timestamp=False):
    """
    plot singal as a function of time.
    

    Parameters
    ----------
    samples : 1D numpy array, real or complex
        sampled signal.
    sample_rate : float, optional
        The sample rate of the signal. The default is 1.0.
    fNum : int, optional
        Figure number to plot into. The default is None which uses the 
        "next unused figure number".
    boundaries : list of int or None, optional
        The boundaries are given as list with two elements (start and end index).
        The signal is only plotted within these given boundaries. A value of None
        specifies the first and last signal sample, respectively. 
        The default is [None, None] and therefore plots the whole signal.
    tit : string, optional
        Title of the plot. The default is 'time signal'.
    save_fig : bool, optional
        should the plot be saved to file? The default is False.
    ffolder : sting, optional
        folder to save figure to. The default is '.'.
    ffname : string, optional
        filename to save figure to. The default is None, which uses the title
        of the plot as filename.
    fformat : string, optional
        format of the saved file, can either be 'png', 'pdf' or 'svg'. 
        The default is 'png'.
    add_timestamp : bool, optional
        should a timestamp be added to the filename? The default is False. 

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
    if fNum:
        fig = plt.figure(fNum, facecolor='white', edgecolor='white')
    else:
        fig = plt.figure(facecolor='white', edgecolor='white')
        
    plt.clf()    
    # if complex input signal -> plot real and imag seperatly
    if np.any(np.iscomplex(samples)):
        plt.subplot(121)
        plt.plot(t, np.real(samples))
        plt.xlabel('time [s]')
        plt.ylabel('amplitude real part')
        plt.grid(visible=True)
        plt.subplot(122)
        plt.plot(t, np.imag(samples))
        plt.xlabel('time [s]')
        plt.ylabel('amplitude imaginary part')
        plt.title(tit)
        plt.grid(visible=True)        
    else:
        plt.plot(t, samples)
        plt.xlabel('time [s]')
        plt.ylabel('amplitude')
        plt.title(tit)
        plt.grid(visible=True)
    
    if save_fig:
        if not ffname:
            ffname = tit
        utils.save_fig(fig, fformat=fformat, folder=ffolder, f_name=ffname, 
                 add_timestamp=add_timestamp)
        
    plt.show()
        
	
def plot_eye(samples, sample_rate=2, bit_rate=1, fNum=None, 
             boundaries=[None, None], tit='eye diagramm', save_fig=False, 
             ffolder='.', ffname=None, fformat='png', add_timestamp=False):
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
        Figure number to plot into. The default is None which uses the 
        "next unused figure number".
    boundaries : list of int or None, optional
        The boundaries are given as list with two elements (start and end index).
        The eye diagram is only plotted within these given boundaries. A value of None
        specifies the first and last signal sample, respectively. 
        The default is [None, None] and therefore the eye diagram contains
        the whole signal.
    tit : string, optional
        Title of the plot. The default is 'eye diagramm'.
    save_fig : bool, optional
        should the plot be saved to file? The default is False.
    ffolder : sting, optional
        folder to save figure to. The default is '.'.
    ffname : string, optional
        filename to save figure to. The default is None, which uses the title
        of the plot as filename.
    fformat : string, optional
        format of the saved file, can either be 'png', 'pdf' or 'svg'. 
        The default is 'png'.
    add_timestamp : bool, optional
        should a timestamp be added to the filename? The default is False.
    
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
    
    if fNum:
        fig = plt.figure(fNum, facecolor='white', edgecolor='white')    
    else:
        fig = plt.figure(facecolor='white', edgecolor='white')    
        
    if np.any(np.iscomplex(samples)):
        plt.subplot(121)
        plt.plot(t, np.real(samples), color = '#1f77b4')
        plt.xlabel('time [s]')
        plt.ylabel('amplitude real part')
        plt.grid(visible=True)        
        plt.title(tit)
        plt.subplot(122)
        plt.plot(t, np.imag(samples), color = '#1f77b4')
        plt.xlabel('time [s]')
        plt.ylabel('amplitude imaginary part')
        plt.grid(visible=True)
        plt.gcf().tight_layout()        
    else:
        plt.plot(t, samples, color = '#1f77b4')
        plt.xlabel('time [s]')
        plt.ylabel('amplitude')
        plt.title(tit)
        plt.grid(visible=True)
    
    if save_fig:
        if not ffname:
            ffname = tit
        utils.save_fig(fig, fformat=fformat, folder=ffolder, f_name=ffname, 
                 add_timestamp=add_timestamp)
        
    plt.show()
    
    
	
def plot_hist(samples, nBins=100):
    #TODO: implement automated histogramm
    # check for complex input??
    pass

def plot_constellation(samples, decimation=1, fNum =None, tit='constellation',
                       hist=False, axMax=None, nBins=128, save_fig=False, 
                       ffolder='.', ffname=None, fformat='png', 
                       add_timestamp=False):
    """
    Plot the constellation diagramm (complex plane) of samples.

    Parameters
    ----------
    samples : 1D numpy array, complex
        samples of the input signal.
    decimation : int, optional
        take only every decimations-th sample of the input signal. The default is 1.
    fNum : int, optional
        figure number of the plot to be created. The default is None which uses the 
        "next unused figure number".
    tit : string, optional
        title of the plot to be created. The default is 'constellation'.
    hist : bool, optional
        should the constellation diagramm be plotted as 2D histogramm? The default is False.
    axMax : float or None
        maximum abolute axis amplitude (equal in x and y axis)
        if None: 1.1 times maximum absolute value of samples (real and imaginalry part) is used
    nBins: int
        number of bins (in each quadrature) for plotting the 2D histogramm
    save_fig : bool, optional
        should the plot be saved to file? The default is False.
    ffolder : sting, optional
        folder to save figure to. The default is '.'.
    ffname : string, optional
        filename to save figure to. The default is None, which uses the title
        of the plot as filename.
    fformat : string, optional
        format of the saved file, can either be 'png', 'pdf' or 'svg'. 
        The default is 'png'.
    add_timestamp : bool, optional
        should a timestamp be added to the filename? The default is False.
    """
    
    samples = samples[0::decimation]
    
    if axMax is None:
        axMax = max(np.abs(samples.real).max(), np.abs(samples.imag).max())*1.1
    
    if fNum:
        fig = plt.figure(fNum, facecolor='white', edgecolor='white')
    else:
        fig = plt.figure(facecolor='white', edgecolor='white')
    
    if hist:
        bins = nBins
        cm = copy.copy(plt.get_cmap("jet"))
        plt.hist2d(samples.real, samples.imag, bins=bins, cmap=cm, cmin=1, density=False)             
    else:     
        plt.plot(samples.real, samples.imag, 'C0.')      
    plt.gca().axis('equal')
    plt.gca().set_axisbelow(True) 
    plt.grid(visible=True)      
    plt.xlim((-axMax, axMax))
    plt.ylim((-axMax,axMax))  
    plt.title(tit)
    plt.xlabel('real part')
    plt.ylabel('imaginary part')  
    
    if save_fig:
        if not ffname:
            ffname = tit
        utils.save_fig(fig, fformat=fformat, folder=ffolder, f_name=ffname, 
                 add_timestamp=add_timestamp)
    
    plt.show()
    
    
def plot_poincare_sphere(samplesX, samplesY, decimation=1, fNum=1, 
                         tit = 'Poincaré sphere', labels=True, save_fig=False, 
                         ffolder='.', ffname=None, fformat='png', 
                         add_timestamp=False):
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
    save_fig : bool, optional
        should the plot be saved to file? The default is False.
    ffolder : sting, optional
        folder to save figure to. The default is '.'.
    ffname : string, optional
        filename to save figure to. The default is None, which uses the title
        of the plot as filename.
    fformat : string, optional
        format of the saved file, can either be 'png', 'pdf' or 'svg'. 
        The default is 'png'.
    add_timestamp : bool, optional
        should a timestamp be added to the filename? The default is False.

    Returns
    -------
    handles :  dict containing following keys
        fig : matplotlib.figure.Figure
            Figure object the signal is plotted to.
        ax : matplotlib.axes._subplots.Axes3DSubplot
            Axes object which contains the Poincaré sphere artists
        line : mpl_tooklits.mplot3d.art3d.Line3d
            Line object which contains the Stokes parameters (S1, S2, S3).
            
    References
    ----------
    [1] https://en.wikipedia.org/wiki/Jones_calculus#Jones_vector
    
    [2] https://en.wikipedia.org/wiki/Stokes_parameters
    
    [3] https://en.wikipedia.org/wiki/Polarization_(waves)#Poincar%C3%A9_sphere
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
        fig = plt.figure(fNum) 
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
        fig = plt.figure(fNum)
        ax = fig.axes[0]
        # find all artists labeled as 'SOP'
        li = ax.findobj(_isSOP)
        # update data
        li[0].set_xdata(s1/s0)
        li[0].set_ydata(s2/s0)
        li[0].set_3d_properties(s3/s0)
        # update plot
        fig.canvas.draw()
        # wait for figure to update
        plt.pause(0.1)

    if save_fig:
        if not ffname:
            ffname = tit
        utils.save_fig(fig, fformat=fformat, folder=ffolder, f_name=ffname, 
                 add_timestamp=add_timestamp)
    
    # prepare return params
    handles = dict()
    handles['fig'] = fig
    handles['ax'] = ax
    handles['line'] = li[0]
    return handles


def place_figures(auto_layout=True, monitor_num=0, nc=4, 
                  nr=3, taskbar_offset=35, figure_toolbar=65):
    """
    Place open figure on screen.
    
    Place open figures on screen unsing specified layout. Basic programmatic
    idea taken from [1].    
    
    Parameters
    ----------
    auto_layout : bool, optional
        The layout is chosen automatically depending on the number of opened 
        figures. Number of figures must not exceed 32. The default is True.
    monitor_num : int, optional
        Monitor onto which the figures are placed. Please note that the order 
        is rather random (i.e. the primary monitor is not necessarily number 0).
        The default is 0.
    nc : int, optional
        Number of coloums used for the layout. Only used if auto_layout=False.
        The default is 4.
    nr : int, optional
        Number of rows used for the layout. Only used if auto_layout=False.
        The default is 3.
    taskbar_offset : int, optional
        Height of the (windows) taskbar which should not be covered by the layout. 
        The taskbar is assumed to be on the bottom of the screen. The height
        depends on many parameters (e.g. monitor scaling, layout of taskbar,... )
        and is therefore to determined by the user. The default is 35.
    figure_toolbar : int, optional
        Height of the toolbar of the individual plot windows. The height
        depends on many parameters (e.g. graphical backend, monitor scaling,... )
        and is therefore to determined by the user. The default is 65.

    Returns
    -------
    None.
    
    References
    ----------
    
    [1] JaeJun Lee (2023). automatically arrange figure windows 
    (https://www.mathworks.com/matlabcentral/fileexchange/48480-automatically-arrange-figure-windows), 
    MATLAB Central File Exchange. Retrieved January 23, 2023. 
    """      
    
    # get figure handles
    figHandle = list(map(plt.figure, plt.get_fignums()))   
    n_fig = len(figHandle)

    if n_fig <= 0:
        raise ValueError('no figures found to place')
    
    # https://stackoverflow.com/questions/3129322/how-do-i-get-monitor-resolution-in-python
    monitor = screeninfo.get_monitors()[monitor_num]
    screen_resolution = [monitor.width, monitor.height-taskbar_offset]
    offset = [monitor.x, monitor.y]

    # auto layout?
    if auto_layout:
        grid = [
            [1,1],[1,2],
            [2,2],[2,2],
            [2,3],[2,3],
            [3,3],[3,3],[3,3],
            [3,4],[3,4],[3,4],
            [4,4],[4,4],[4,4],[4,4],
            [4,5],[4,5],[4,5],[4,5],
            [4,6],[4,6],[4,6],[4,6],
            [4,7],[4,7],[4,7],[4,7],
            [4,8],[4,8],[4,8],[4,8]
            ]
       
        if n_fig > len(grid)*2:
            raise ValueError('more figures opened than layout options available')        
        
        # portrait mode
        if screen_resolution[0] < screen_resolution[1]:
            nc = grid[n_fig-1][0]
            nr = grid[n_fig-1][1]
        # landscape mode
        else:
            nc = grid[n_fig-1][1]
            nr = grid[n_fig-1][0]
    # manual layout
    else:
        if (nc * nr) < n_fig:
            raise ValueError(f'more figures opened ({n_fig}) than rows times coloumns given ({nc*nr}): try to increase numbers or switch to auto layout mode')
    
    fig_width = screen_resolution[0]/nc 
    fig_height = screen_resolution[1]/nr-figure_toolbar 

    fig_cnt = 0
    backend = plt.get_backend()
    for r in range(nr):
        for c in range(nc):
            if fig_cnt >= n_fig:
                break
            # move figure to required monitor            
            if backend == 'TkAgg':
                figHandle[fig_cnt].canvas.manager.window.wm_geometry(f"+{offset[0]+1}+{offset[1]+1}")
            elif backend == 'WXAgg':
                figHandle[fig_cnt].canvas.manager.window.SetPosition(offset[0]+1, offset[1]+1)
            else:
                figHandle[fig_cnt].canvas.manager.window.move(offset[0]+1, offset[1]+1)  
            # get DPI (for scaling)
            fig_dpi = figHandle[fig_cnt].get_dpi()  
            # set figure to desired size
            figHandle[fig_cnt].set_figheight(fig_height / fig_dpi)
            figHandle[fig_cnt].set_figwidth(fig_width / fig_dpi)  
            figHandle[fig_cnt].tight_layout() 
            # solution for different backends taken
            # from https://stackoverflow.com/questions/7449585/how-do-you-set-the-absolute-position-of-figure-windows-with-matplotlib
            if r == 0:                
                x_pos = int((fig_width*c)/(fig_dpi/100)+offset[0])
                y_pos = int((fig_height*r)/(fig_dpi/100)+offset[1])
            else:                
                x_pos = int((fig_width*c)/(fig_dpi/100)+offset[0])
                y_pos = int(((fig_height+figure_toolbar)*r)/(fig_dpi/100)+offset[1])
            if backend == 'TkAgg':
                figHandle[fig_cnt].canvas.manager.window.wm_geometry(f"+{x_pos}+{y_pos}")
            elif backend == 'WXAgg':
                figHandle[fig_cnt].canvas.manager.window.SetPosition(x_pos, y_pos)
            else:
                figHandle[fig_cnt].canvas.manager.window.move(x_pos, y_pos)
            fig_cnt += 1