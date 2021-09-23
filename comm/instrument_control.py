import sys
import visa
import numpy as np
import matplotlib.pyplot as plt
#import math
import time
import logging


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


def get_samples_Tektronix_MSO6B(channels=[1], ip_address='192.168.1.20',number_of_bytes = 1,log_mode = False):   
    """
    get_samples_Tektronix_MSO6B
    
    Function for reading samples from  Tektronix MSO68 Scope

    Parameters
    ----------
    channels : list of integers, optional
        iterable containing the channel numbers to fetch data samples from device. The default is [1].
        For more chennels use [1,2,...]
        Minimum number of channels is 1
        Maximum number of channels is 4
        Ensure that the acquired channels are activated at the scope
    address : string, optional
        IP Adress of device. The default is '192.168.1.20'.
    number_of_bytes: integer, optional
        Defines the length of the requested data from the scope in bytes.
        Allowed are 
        1 Byte (signed char), 2 Bytes (signed short) or 4 Bytes (long)
    log_mode: boolean, optional
        Specifies whether a log file should be created or not

    Returns
    -------
    sample_rate : float
        actual sample rate of returned samples.
    wfm : list of numpy arrays, float
        each list element constains the samples of a requested channel as numpy float array.

    Errors
    -------
    Type Error: 
        Will be raised when a wrong data type is used for the input parameter
        -> Possible errors
            -> Channels is not of type list
            -> Items of channels are not of type integer
            -> ip_address is not of type string
            -> number_of_bytes is not integer

    Value Error:
        Will be raised when the input parameter is in an wrong range
        -> Possible errors
            -> Too much channels are used. Maximum is 4
            -> Too less channels are used. Minimus is 1 
            -> Channel numbers must be between 1 and 4
            -> Wrong number of bytes (1 Byte (signed char), 2 Bytes (signed short) or 4 Bytes (long))

    Exception:
        Will be raised by diverse errors
        -> Possible errors
            -> No connection to the scope
            -> Required channels are not activated at the scope

    """

    # =============================================================================
    #  Create logger which writes to file
    # ============================================================================= 
    # Create logger
    logger = logging.getLogger(__name__)

    # Set the log level
    logger.setLevel(logging.INFO)

    # Create file handler and standard output handler (terminal output)
    file_handler = logging.FileHandler('{0}.log'.format(__name__))
    if log_mode == False:
        file_handler.setLevel(51)
    else:
        file_handler.setLevel(logging.INFO)

    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.ERROR)

    # Set format of the logs with formatter
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(name)s :: Line No %(lineno)d:: %(message)s')

    # Adding formatter to handler
    file_handler.setFormatter(formatter)
    stdout_handler.setFormatter(formatter)

    # Adding handler to logger
    logger.addHandler(file_handler)
    logger.addHandler(stdout_handler)

    # =============================================================================
    #  Check inputs of correctnes
    # ============================================================================= 

    try:
        if not isinstance(channels, list):
            raise TypeError('Type of channels must be list')

        if not isinstance(ip_address, str):
            raise TypeError('Type of ip_address must be string')

        if not all(isinstance(x, int) for x in channels):
            raise TypeError('Type of channels items must be integers')

        if not isinstance(number_of_bytes, int):
            raise TypeError('Type of number_of_bytes must be integer')

        if len(channels) > 4:
            raise ValueError('Too much channels ({0}). The Scope has only 4 input channels'.format(len(channels)))

        if len(channels) < 1:
            raise ValueError('Too less channels ({0}). Use at least one channel'.format(len(channels)))

        if any(ch_number > 4 for ch_number in channels) > 4 or any(ch_number < 1 for ch_number in channels):
            raise ValueError('Channel numbers must be betwenn 1 and 4')

        if number_of_bytes not in [1,2,4]:
            raise ValueError('Wrong number of bytes. Only 1 (signed char), 2 (signed short) or 4 (long) are allowed')

    except Exception as e:
        logger.error('{0}'.format(e))
        exit()
    
    # =============================================================================
    #  importing visa for communication with the AWG device
    # ============================================================================= 
    # create ressource
    rm = visa.ResourceManager('@py')

    # open connection to AWG
    logger.info("Create IP connection with " + str(ip_address))
    try:
        scope = rm.open_resource('TCPIP::' + ip_address + '::INSTR')
    except Exception as e:
        logger.error('No connection possible. Check TCP/IP connection \n  {0}'.format(e))
        exit()
    
    logger.info("Device properties: " + str(scope.query('*IDN?')))


    # =============================================================================
    #  Settings for the Scope
    # =============================================================================  

    # Generate waveform vector
    wfm = []

    # Set the last datapoint of the waveform which will be transmitted. 
    # If this value is bigger than the actual length of the waveform, the Curve? function catches all samples
    scope.write('DATA:STOP 62500000')

    # See if scope is running or not
    if scope.query('ACQuire:STATE?')[0] == '1':
        # Setting the restart condition of the scope after acqusition - True means restart
        is_running = True

        # start a single acquisition
        scope.write('ACQuire:STOPAfter SEQuence')
        scope.write('ACQuire:STATE RUN')  

        # Loop to ascertain that the scope has finished the acquisition process
        while scope.query('ACQuire:STATE?')[0] == '1':
            pass

    else:
        # Setting the restart condition of the scope after acqusition - False means stay by stop
        is_running = False
    
    # Setting outgoing data format (In this case signed binary data)
    scope.write('DATA:ENCDG SRIbinary')

    # Get sample rate of the scope
    sample_rate = scope.query('HORizontal:MODE:SAMPLERate?')

    # Set datatype for acquisition  ( b (signed char),h (signed short) or l (long) )
    # For information see documentation of struct module

    if number_of_bytes == 1:
        acq_data_type = 'b'
    elif number_of_bytes == 2:
        acq_data_type = 'h'
    else:
        acq_data_type = 'l'
        
    logger.info("Used datatype: " + acq_data_type)

    # Read the channels
    for ch in (channels):
        # Select waveform source
        scope.write('DATA:SOURCE CH{0:d}'.format(ch))

        # Setting number of bytes per waveformpoint (sample) 
        scope.write('WFMOutpre:BYT_Nr {0:d}'.format(number_of_bytes))

        # Reading waveform data and write them as numpy array to list
        #scope.write('Curve?')
        try:
            tmp = scope.query_binary_values('Curve?',datatype=acq_data_type ,is_big_endian=False, container=np.array)
        except Exception as e:
            logger.error('Channel {0:d} seems not activated \n '.format(ch))
            exit()

        print(tmp)

        # Reading vertical scaling of the scope (Voltage per div)
        ver_scale = float(scope.query('CH{0:d}:SCAle?'.format(ch),delay = 0.5))

        # Reading vertical position ( Y-Position on scope screen)
        ver_position = float(scope.query('CH{0:d}:POSition?'.format(ch),delay = 0.5))

        # Reading offset from scope
        offset = float(scope.query('CH{0:d}:OFFSet?'.format(ch),delay = 0.5))

        # Scale amplitude
        #wfm.append(5 * ver_scale / (2 ** (float(number_of_bytes) * 7)) * tmp  + ver_offset)
        # if number_of_bytes == 4:
        #     wfm.append(ver_scale * (tmp - ver_position) + offset)
        # else:
        wfm.append(ver_scale * (5 / (2 ** (float(number_of_bytes) * 8 - 1)) * tmp - ver_position) + offset)
    

    # Restart the scope
    if is_running:
        scope.write('ACQuire:STOPAfter RUNSTOP')
        scope.write('ACQuire:STATE RUN')  


    # closing scope connection
    scope.close()
   
    # closing resource manager 
    rm.close()  

    return sample_rate, wfm



def write_samples_Tektronix_AWG70002B(samples, ip_address='192.168.1.21', sample_rate=[250e6], amp_pp=[0.5], channels=[1], out_filter=['normal'],log_mode = False):
    
    # TODO: Write/change the docstring for the method. Find Filter
    # TODO: Change float to integer (Maybe?)

    """
    write_samples_AWG70002B
    
    Function for writing samples to an Tektronix AWG70002B Series 20GHz Function/Arbitrary Waveform Generator

    Parameters
    ----------
    samples : numpy array, n_outputs x n_samples , float
        samples to output, to be scaled between -1 and 1 (values outside this range are clipped).
    ip_address : string, optional
        DESCRIPTION. The default is '192.168.1.21'. Currently, only LAN connection is supported.
    sample_rate : list of floats, optional
        sample rate of the outputs. The default is [250e6]. Must be between 1.49 kSamples/s and 8 GSsamples/s
    amp_pp : list of floats, optional
        peak-to-peak output amplitude of individual channels in units of Volt. The default is [0.5].
    channels : list of int, optional
        channels to be programmed and output. The default is [1]. For two channels input [1,2]
    out_filter : list of strings, optional
        used output filter of each channel ['normal', 'off', 'step']. The default is ['normal'].

    Returns
    -------
    None.

    """
    
    # =============================================================================
    #  Create logger which writes to file
    # ============================================================================= 
    # Create logger
    logger = logging.getLogger(__name__)

    # Set the log level
    logger.setLevel(logging.INFO)

    # Create file handler and standard output handler (terminal output)
    file_handler = logging.FileHandler('{0}.log'.format(__name__))
    if log_mode == False:
        file_handler.setLevel(51)
    else:
        file_handler.setLevel(logging.INFO)

    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.ERROR)

    # Set format of the logs with formatter
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(name)s :: Line No %(lineno)d:: %(message)s')

    # Adding formatter to handler
    file_handler.setFormatter(formatter)
    stdout_handler.setFormatter(formatter)

    # Adding handler to logger
    logger.addHandler(file_handler)
    logger.addHandler(stdout_handler)



    # =============================================================================
    #  Check inputs of correctnes
    # ============================================================================= 
    
    try:
        if not (isinstance(sample_rate, list) and 
            isinstance(amp_pp, list) and isinstance(channels, list) and 
            isinstance(out_filter, list)):
            raise TypeError('Input parameters are not lists...')

        if len(channels) > 2:
            raise ValueError('To much channels ({0}). The AWG has only 2 output channels'.format(
                                                                                         len(channels)))

        if np.iscomplex(samples[0:2]).any():
            raise TypeError('No complex numbers allowed. If you want to use complex values, assign the real part and imaginary part seperately to channel 1 and channel 2')

        if np.isnan(samples[0:2]).any() or np.isinf(samples[0:2]).any():
            raise ValueError('No NaN or Inf values are allowed in the sample vector!')

        if not(len(channels) == len(samples) == len(amp_pp)):
            raise ValueError('Number of channels ({0}), number of signal vectors ({1}) and number of amplitudes ({2}) must be the same!'.format(
                                                                                                                                        len(channels),
                                                                                                                                        len(samples),
                                                                                                                                        len(amp_pp)))

        if sample_rate[0] > 8e9 or sample_rate[0] < 1.49e3:
            raise ValueError('Sample rate must be between 1.49 kSamples/s and 8 GSsamples/s')

        if any(ch_amp > 0.5 for ch_amp in amp_pp) or any(ch_amp < 0.25 for ch_amp in amp_pp):
            raise ValueError('Amplitudes must be between 0.25 and 0.5 (peak to peak)')

        if not isinstance(samples, np.ndarray):
            raise TypeError('Samples has to be from type numpy array. Actual type: {0}'.format(
                                                                                       type(samples)))

        

    except Exception as e:
        logger.error('{0}'.format(e))
        exit()


    # TODO: Search for min and max possible sample rate of the AWG and write an if statement to catch a correct input of sample rate
    # TODO: Find the maximum of the samples vector and write an if statement to catch a correct input
    # TODO: Check for correct input data types 
    # TODO: Writing some query commands to get feedback from the scope


    # =============================================================================
    #  importing visa for communication with the AWG device
    # ============================================================================= 
    # create resource 
    rm = visa.ResourceManager('@py')
    # open connection to AWG
    logger.info("Create IP connection with " + str(ip_address))
    try:
        awg = rm.open_resource('TCPIP::' + ip_address + '::INSTR')
    except Exception as e:
        logger.error('No connection possible. Check TCP/IP connection \n  {0}'.format(e))
        exit()
    
    # Setting timeout
    awg.timeout = 20000

    logger.info("Device properties: " + str(awg.query('*IDN?')))

    # =============================================================================
    #  Clipping the signal vector (Range from -1 to 1)
    # ============================================================================= 
    if np.amax(np.amax(np.abs(samples))) > 1:
        logger.warning("Samples have been clipped")
    else:
        logger.info("Samples have not been clipped")
    
    samples_clipped =(np.clip(samples,-1,1))

    if samples_clipped.ndim == 1:
        samples_clipped = samples_clipped[np.newaxis,...]
    #samples_clipped = samples_clipped.tolist()
    
    logger.debug(type(samples_clipped[0]))

    # =============================================================================
    #  Settings for the AWG
    # =============================================================================  

    # Setting AWG to STOP
    logger.info("Set AWG to stop")
    awg.write('AWGCONTROL:STOP:IMMEDIATE')
    
    # Delete old waveform
    logger.info("Delete old waveform")
    awg.write('WLIST:WAVEFORM:DELETE "Python_waveform_AWG_1"')
    awg.write('WLIST:WAVEFORM:DELETE "Python_waveform_AWG_2"')

    # decoupling of the two channels
    logger.info("Decouple channels")
    awg.write('INSTrument:COUPLe:SOURce OFF')
    
    # Output deactivate
    logger.info("Deactivate output")
    awg.write('OUTPUT1:STATE OFF')
    awg.write('OUTPUT2:STATE OFF')

    # Setting sample rate
    logger.info("Set sample rate to: " + str(sample_rate[0]))
    awg.write('CLOCK:SRATE {0:f}'.format(sample_rate[0]))
    awg.query('*OPC?')[0]


    for ch_idx, ch in enumerate(channels):
        logger.info("\n---Channel {0:d}---".format(ch))
        
        # Output deactivate
        logger.info("Deactivate output")
        awg.write('OUTPUT{0:d}:STATE OFF'.format(ch))

        # Create new waveform
        logger.info("Create new waveform")
        logger.info("Name of new waveform: Python_waveform_AWG_{0:d}".format(ch))
        logger.info("Length of new waveform: {0:d}".format(len(samples_clipped[ch_idx])))

        # Send data to AWG
        length_of_samples = len(samples_clipped[ch_idx])
        awg.write('WLISt:WAVeform:NEW "Python_waveform_AWG_{0:d}",{1:d}'.format(ch,length_of_samples))
        logger.info("Write data to waveform")
            
        awg.write_binary_values('WLIST:WAVEFORM:DATA "Python_waveform_AWG_{0:d}",'.format(ch), samples_clipped[ch_idx], datatype='f')
        print(awg.query('WLIST:WAVEFORM:DATA? "Python_waveform_AWG_{0:d}",0'))
        # WLISt:WAVeform:DATA[:I]? <wfm_name>[,<StartIndex>[,<Size>]])

        # Adding the waveform to an output
        logger.info("Add Python_waveform_AWG_{0:d} to output {0:d}".format(ch))
        awg.write('SOURCE{0:d}:CASSET:WAVEFORM "Python_waveform_AWG_{0:d}"'.format(ch)) 

        # Setting parameters of the waveform
        #  Amplitude (peak to peak)
        logger.info("Set amplitude (peak to peak) of output {0:d} to {1:f}".format(ch,amp_pp[ch_idx]))
        awg.write('SOURCE{0:d}:VOLTAGE:LEVel:IMMediate:AMPLITUDE {1:f}'.format(ch,amp_pp[ch_idx]))

        # Activating outputs
        logger.info("Activate output {0:d}".format(ch))
        awg.write('OUTPUT{0:d}:STATE ON'.format(ch))

    # Starting playback 
    logger.info("\nSet AWG to run")
    awg.write('AWGCONTROL:RUN:IMMEDIATE')

    # closing AWG connection
    awg.close()
   
    # closing resource manager 
    rm.close()  
