
import visa
import numpy as np
import matplotlib.pyplot as plt
#import math
import time


def get_samples_DLM2034(channels=(1), address='192.168.1.12'):
    """
    Parameters
    ----------
    channels : iterable, int, optional
        iterable containing the channel numbers to fetch data samples from device. The default is (1).
    address : string, optional
        IP Adress of device. The default is '192.168.1.12'.

    Returns
    -------
    sample_rate : float
        actual sample rate of returned samples.
    wfm : list of numpy arrays, float
        each list element constains the samples of a requested channel as numpy float array.

    """

    # create resource 
    rm = visa.ResourceManager('@py')
    # print(rm.list_resources())
    #rm = visa.ResourceManager()
    
    # open connection to scope
    scope = rm.open_resource('TCPIP::' + address + '::INSTR')
    # set number of bytes to retireve at once...TODO: find reasonable value
    scope.chunk_size = 2000
    
    # check instrument IDN
    idn = scope.query('*IDN?')
    print(idn)
    
    # check if device is running or stopped?
    # using Condition register, see. comm. interface manual
    # p. 6-5, CHECK!!!
    # running, if LSB is set -> cond. register is odd
    running = float(scope.query('STATus:CONDition?')) % 2
    
    if running:
        # start a single acquisition
        scope.write('TRIGger:MODE SINGle')
        busy = 1
        while busy:
            busy = float(scope.query('STATus:CONDition?')) % 2
            
    wfm = []
    
    for idx, channel in enumerate(channels):
        # set channel to retrieve
        scope.write('WAVeform:TRACe {:d}'.format(int(channel)))
        
        # set waveform format to int16
        scope.write('WAVeform:FORMat WORD')
        
        # get range
        range = float(scope.query('WAVeform:RANGe?').split(sep=" ")[1])
        
        # get offset
        offset = float(scope.query('WAVeform:OFFSet?').split(sep=" ")[1])
        
        # get waveform
        scope.write('WAVeform:SEND?')
        tmp = scope.read_binary_values(datatype='h', is_big_endian=False, container=np.array)
        
        #scale waveform according to (Range × data ÷ division*) + offset)...
        #division is fix: 3200 for format WORD and 12.5 for BYTE format, see. comm. interface manual
        #p. 5-290
        tmp = (range * tmp / 3200) + offset
        
        wfm.append(tmp)
    
    # get samplerate
    sample_rate = float(scope.query('WAVeform:SRATe?').split()[1])
    
    # set initial state
    if running:
        # reset device condition
        scope.write('TRIGger:MODE NORM')
    
    #print(float(scope.query('WAVeform:OFFSet?').split()[1]))
    
    #t = np.linspace(0, (len(wfm)-1) / sample_rate, len(wfm))
    #plt.plot(t, wfm)
    #plt.grid(True)
    
    # close connection and delete objects
    rm.close()
    del rm
    del scope
    
    # change list of nparrays to single nparray
    wfm = np.asarray(wfm)
    
    return sample_rate, wfm


def write_samples_AWG33522A(samples, ip_address='192.168.1.44', sample_rate=[250e6], offset=[0.0], amp_pp=[1.0], channels=[1], out_filter=['normal']):
    """
    write_samples_AWG33522A
    
    Function for writing samples to an Agilent/Keysight 33500 Series 30MHz Function/Arbitrary Waveform Generator

    Parameters
    ----------
    samples : numpy array, n_outputs x n_samples , float
        samples to output, to be scaled between -1 and 1 (values outside this range are clipped).
    ip_address : string, optional
        DESCRIPTION. The default is '192.168.1.44'. Currently, only LAN connection is supported.
    sample_rate : list of floats, optional
        sample rate of the individual outputs. The default is [250e6]. Range: 1µSa/s to 250 MSa/s, limited to 62.5 MSa/s if out_filter is OFF.
    offset : list of floats,, optional
        output DC offset of individual channels in V. The default is [0.0].
    amp_pp : list of floats, optional
        peak-to-peak output amplitude of individual channels in units of Volt. The default is [1.0].
    channels : list of int, optional
        channels to be programmed and output. The default is [1].
    out_filter : list of strings, optional
        used output filter of each channel ['normal', 'off', 'step']. The default is ['normal'].

    Returns
    -------
    None.

    """
    
    if not (isinstance(sample_rate, list) and isinstance(offset, list) and 
            isinstance(amp_pp, list) and isinstance(channels, list) and 
            isinstance(out_filter, list)):
        raise TypeError('input parameters are not lists...')
        
    if not (len(sample_rate) == len(offset) == len(amp_pp) 
            == len(channels) == len(out_filter)):
        raise TypeError('length of parameter lists are not equal...')
    
    if not isinstance(samples, np.ndarray):
        raise TypeError('samples has to be a numpy array...')
    
    for idx, out_filt in enumerate(out_filter):
        if (sample_rate[idx] > 62.5e6) and (out_filt.lower() == 'off'):
            raise ValueError('If sample rate is above 62.5 MHz, output filter has to be set to "normal" or "step"...')
            
    # TODO: add more input parameter checks

            
    # =============================================================================
    #  importing visa for communication with the device
    # ============================================================================= 
    # create resource 
    rm = visa.ResourceManager('@py')
    # open connection to AWG
    awg = rm.open_resource('TCPIP::' + ip_address + '::INSTR')   

    # selecting byte order , used to make binary data point transfers in the block mode Swapped(LSB) or Normal(MSB)
    # SWAPped byte order,(LSB) of each data point is assumed first. Most computers use the "swapped" byte order.
    awg.write(':FORMat:BORDer %s' % ('SWAPped'))
    
    # clip samples and format to list of int16 numbers
    samples = np.round(np.clip(samples, -1.0, 1.0) * 32767).astype(int)
    # ensure that samples is a nested list, even if ndim == 1
    if samples.ndim == 1:
        samples = samples[np.newaxis,...]
    samples = samples.tolist()
    
    #loop over up to 2 channels
    for ch_idx, ch in enumerate(channels):

        # disable channel coupling
        awg.write(':SOUR{0:d}:VOLT:LEVel:IMMediate:COUP:STAT OFF'.format(ch))
        awg.write(':SOUR{0:d}:RATE:COUP:STAT OFF'.format(ch))

        # output to off is necessary, otherwise the Amplitude is automatically set to 10V, which is dangerous 
        # output set to off/ output will be automatic activated loading up data
        awg.write(':OUTP{0:d} OFF'.format(ch))
        
        # clearing the waveform memory of the specified channel
        awg.write(':SOUR{0:d}:DATA:VOLatile:CLEar'.format(ch))
        
        # writing values representing DAC codes into waveform volatile memory, as binary block data/ list of integer samples from -32767 to +32767.
        # loading data into the AWG as arb%d, where d = 1 or 2 taken from the list of channel
        awg.write_binary_values(':SOUR{0:d}:DATA:ARBitrary:DAC arb{0:d},'.format(ch), samples[ch_idx], datatype='h', is_big_endian=False)
        
        # setting output waveform of channel to ARB
        awg.write(':SOUR{0:d}:FUNC:SHAP:ARB "arb{0:d}"'.format(ch))
        #awg.write(':SOUR%d:FUNC:SHAP:ARBitrary "arb%d"' % (ch, ch))
       
        # applying output filter mode
        awg.write(':SOUR{0:d}:FUNC:SHAP:ARB:FILT {1:s}'.format(ch, out_filter[ch_idx].upper()))
 
        # applying sample rate, amplitude and Offset        
        awg.write(':SOUR{0:d}:APPL:ARB {1:g},{2:g}, {3:g}'.format(ch, sample_rate[ch_idx], amp_pp[ch_idx], offset[ch_idx]))
        #awg.write(':SOURce%d:APPLy:ARBitrary %s,%s,%s' % (ch, sample_rate[ch_idx], amp_pp[ch_idx], offset[ch_idx]))
        
        # wait a moment to have the output to turned on
        time.sleep(0.1)
        
               
    awg.write(':SOUR{0:d}:FUNC:ARB:SYNC'.format(ch))  # synchronising channels
        
    awg.close() # closing AWG
    rm.close()  # closing resource manager 