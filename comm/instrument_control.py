import sys
import pyvisa as visa
import numpy as np
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
            -> channels is not of type list
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



def write_samples_Tektronix_AWG70002B(samples, ip_address='192.168.1.21', sample_rate=[250e6], amp_pp=[0.5], channels=[1],log_mode = False):


    """
    write_samples_AWG70002B
    
    Function for writing samples to an Tektronix AWG70002B Series 20GHz Function/Arbitrary Waveform Generator

    Parameters
    ----------
    samples : numpy array, n_outputs x n_samples , float
        samples to output, to be scaled between -1 and 1 (values outside this range are clipped).
        Without clipping the AWG would clip the waveform.
        Only real numbers are allowed. To use complex numbers assign the real and imaginray part to different channels.
        Maximum vector length is 234e6.
    ip_address : string, optional
        The default is '192.168.1.21'. Currently, only LAN connection is supported.
    sample_rate : list of floats, optional
        sample rate of the outputs. The default is [250e6]. Must be between 1.49 kSamples/s and 8 GSsamples/s
    amp_pp : list of floats, optional
        peak-to-peak output amplitude of individual channels in units of Volt. The default is [0.5].
        For two channels enter format [x.x,y.y]
    channels : list of int, optional
        channels to be programmed and output. The default is [1]. For two channels input [1,2]
    log_mode : Bool, optional
        When True a log file will be created (Default = False)
        The log file includes error messages and infos about the program flow

    Returns
    -------
    None.

    Errors
    ------
    Type Error: 
        Will be raised when a wrong data type is used for the input parameter
        -> Possible errors
            -> Parameters are not of type list
            -> Items of channels, amp_pp or sample_rate are not of type integer
            -> Items of samples are not of type np.array
            -> Items of samples are of type complex.
            -> ip_address is not a string

    Value Error:
        Will be raised when the input parameter is in an wrong range
        -> Possible errors
            -> The samples np.arrays contains NaN or Inf
            -> The lengths of amp_pp, channels and samples have not the same length
            -> The peak to peak voltage is not between 0.25V and 0.5V
            -> The sampling_rate ist not between 1.49e3 and 8e9
            -> Channel designation is wrong

    Exception:
        Will be raised by diverse errors
        -> Possible errors
            -> No connection to the AWG

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
            isinstance(amp_pp, list) and isinstance(channels, list)):
            raise TypeError('Input parameters are not lists...')

        if not all(isinstance(x, int) for x in amp_pp):
            TypeError('amp_pp items must be of type integer')

        if not all(isinstance(x, int) for x in sample_rate):
            TypeError('sample_rate items must be of type integer')   

        if not all(isinstance(x,int) for x in channels):
            TypeError('channels items must be of type integer')     

        if not isinstance(samples, np.ndarray):
            raise TypeError('Samples has to be from type numpy array. Actual type: {0}'.format(type(samples)))

        if not isinstance(ip_address,str):
            raise TypeError('ip_address must be of type string')  

        if np.iscomplex(samples[0:2]).any():
            raise TypeError('No complex numbers allowed. If you want to use complex values, assign the real part and imaginary part seperately to channel 1 and channel 2')

        if np.isnan(samples[0:2]).any() or np.isinf(samples[0:2]).any():
            raise ValueError('No NaN or Inf values are allowed in the sample vector!')

        # if len(samples) > 234_000_000:
        #     raise ValueError("Maximum length of sample vector is 234e6")


        if len(channels) > 2:
            raise ValueError('To much channels ({0}). The AWG has only 2 output channels'.format(
                                                                                         len(channels)))

        if not(len(channels) == len(samples) == len(amp_pp)):
            raise ValueError('Number of channels ({0}), number of signal vectors ({1}) and number of amplitudes ({2}) must be the same!'.format(
                                                                                                                                        len(channels),
                                                                                                                                        len(samples),
                                                                                                                                        len(amp_pp)))

        if any(ch_num > 2 for ch_num in channels) or any(ch_num < 1 for ch_num in channels):
            raise ValueError('Channel designation must be between 1 and 2')

        if sample_rate[0] > 8e9 or sample_rate[0] < 1.49e3:
            raise ValueError('Sample rate must be between 1.49 kSamples/s and 8 GSsamples/s')

        if any(ch_amp > 0.5 for ch_amp in amp_pp) or any(ch_amp < 0.25 for ch_amp in amp_pp):
            raise ValueError('Amplitudes must be between 0.25 and 0.5 (peak to peak)')


        

    except Exception as e:
        logger.error('{0}'.format(e))
        exit()

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
    awg.timeout = 20_000

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


def get_samples_HP_71450B_OSA (traces = ['A'], GPIB_address='13',log_mode = False, single_sweep = False):


    """
    get_samples_HP_71450B_OSA
    
    Function for reading samples from a HP_71450B optical spectrum analyzer
    
    Parameters
    ----------
        traces: list of stings, optional (default = ['A'])
            Insert here the wanted traces from the OSA as a list of strings. 
            The three traces of the OSA are A, B, and C. It is also possible to use lower case.

        GPIB_address : string, optional (default = '13')
            The address GPIB address of the OSA.
        
        log_mode: boolean, optional (default = False)
            Enables a log file for this method.
            
        single_sweep = boolean, optional (default = False)
            Starts a new sweep and stops after acquisition. Keeps the OSA in Single mode.
            Be careful, because saved traces can be overwritten by this!
            By default the program will acquire the traces, while the OSA is sweeping. The sweeping process is slow enough
            for this.

    Returns
    -------
        trace_information: dict
            Consist of dicts which contains the acquired trace data, amplitude unit and wavelength information.
            To access the dict use:
            >Name of object<[>Name of trace<][>Name of data<]
                -> Name of Trace: 
                    -> A : Trace A
                    -> B : Trace B
                    -> C : Trace C
                -> Name of data:
                    -> Trace_data   : (np.array) Contains numpy array with trace data
                    -> Unit         : (string) Contains the unit of the trace data
                    -> Sensitivity  : (float) Contains the amplitude sensitivity of the spectrum. Is always in dBm
                    -> Start_WL     : (float) Contains the start wavelength of the spectrum (in nm)
                    -> Stop_WL      : (float) Contains the stop wavelength of the spectrum (in nm)
                    -> Resolution_BW: (float) Contains the resolution bandwidth of the spectrum
                    -> WL_vector    : (np.array) Contains a numpy array with an even spaced wavelength vector between Start_WL and Stop_WL (in nm)
            
    Errors
    -------
        Type Error: 
            Will be raised when a wrong data type is used for the input parameter
            -> Possible errors
                -> traces is not of type list
                -> Items of traces are not of type integer
                -> ip_address is not of type string
                -> number_of_bytes is not integer

        Value Error:
            Will be raised when the input parameter is in an wrong range
            -> Possible errors
                -> Too many traces are used. Maximum is 3
                -> Too few traces are used. Minimus is 1 
                -> Trace numbers must be between 1 and 3

        Exception:
            Will be raised by diverse errors
            -> Possible errors
                -> No connection to the scope
                -> Required traces are not activated at the scope

    """
    
    # TODO: Create a way to save the data from the scope to file
    # TODO: File should be similar to LabView file
    
    
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
        if not isinstance(traces, list):
            raise TypeError('Type of traces must be list')

        if not isinstance(GPIB_address, str):
            raise TypeError('Type of GPIB_address must be string')

        if not all(isinstance(x, str) for x in traces):
            raise TypeError('Type of traces items must be strings')

        if len(traces) > 3:
            raise ValueError('Too many traces ({0}). The OSA has maximal 3 traces'.format(len(traces)))

        if len(traces) < 1:
            raise ValueError('Too less traces ({0}). Use at least one trace'.format(len(traces)))

        # Change traces items to upper case
        traces = [each_string.upper() for each_string in traces]

        if any((trace_name not in ['A','B','C']) for trace_name in traces):
            raise ValueError('Wrong trace naming. Traces are named with A, B or C. Lower case is also accepted')

    except Exception as e:
        logger.error('{0}'.format(e))
        return sys.exit(0)

    # =============================================================================
    #  importing visa for communication with the OSA
    # ============================================================================= 

    rm = visa.ResourceManager()

    # open connection to AWG
    logger.info("Create GPIB connection with " + str(GPIB_address))
    try:
        osa = rm.open_resource('GPIB0::' + GPIB_address + '::INSTR')
    except Exception as e:
        logger.error('No connection possible. Check GPIB connection \n  {0}'.format(e))
        return sys.exit()

    # Setting timeout
    # osa.timeout = 20_000


    # =============================================================================
    #  Settings for the analyzer
    # =============================================================================  
    
    ######
    # The page numbers refer to Programmer's Guide of HP 71450B
    ######

    # Query sweep mode 
    # Page 7-477
    # current_sweepmode = osa.query('SWPMODE?').rstrip('\n')
    
    # Check if OSA is sweeping
    # Page 7-476
    # is_running = osa.query('SWEEP?').rstrip('\n')

    if single_sweep:
        # Start a single sweep
        # Page 7-443
        osa.write('SNGLS')
        # Wait till sweep is done
        # Page 7-121
        while not osa.query('DONE?').rstrip('\n') == '1':
            pass

    # if is_running == '1' and current_sweepmode == 'CONTS':
    #     # Start a single sweep
    #     # Page 7-443
    #     osa.write('SNGLS')
    #     # Wait till sweep is done
    #     # Page 7-121
    #     while not osa.query('DONE?').rstrip('\n') == '1':
    #         pass

    # Set datatype of acquisition (Word -> 2 Bytes per sample)
    # Page 7-232 -> 7-234
    osa.write('MDS W')

    # Set type of transmission (I-Block Data field)
    # The I block data field transmit the trace data in binary format
    # Page 7-478 -> 7-480
    osa.write('TDF I')

    # Check amplitude unit
    # Page 7-58 -> 7-59
    amplitude_unit = osa.query('AUNITS?').rstrip('\n') 

    # Check if amplitude unit is logarithmic or linear
    if amplitude_unit in ['V','W']:
        is_log = False
    else:
        is_log = True

    # Create dict with the traces
    trace_information = dict.fromkeys(traces)

    # Read start wave length
    # Page 7-457 -> 7-458
    # Convert from m to nm Page 1-14
    # With restrip(), the terminator \n will be removed
    start_wl = 1e9 * float(osa.query('STARTWL?').rstrip('\n') )

    # Read stop wave length
    # Page 7-464 -> 7-465
    # Convert from m to nm Page 1-14
    stop_wl = 1e9 * float(osa.query('STOPWL?').rstrip('\n') )

    # Loop through traces
    for trace_id,trace in enumerate(traces):

        # Create dictionary for trace data and wave length information
        data_dict = {'Trace_data':[],'Unit':[],'Sensitivity':[],'Start_WL':[],'Stop_WL':[],'Resolution_BW':[], 'WL_Vector':[]}

        # Setting length of Trace
        # Page 7-506 -> 7-507
        osa.write('TRDEF TR{0:s},2048'.format(trace))

        # Read trace
        # Page 7-499 -> 7-502
        # h is 2 bytes (signed short)
        tmp = osa.query_binary_values('TR{0:s}?'.format(trace),datatype='h' ,is_big_endian=True, container=np.array,data_points = 2048)

        # Convert measument units to parameter units
        # Page 2-8
        if is_log:
            # One measurement unit is equal to one hundreth of a dBm
            # To get the dBm the trace data from the scope has to be divided by 100
            data_dict['Trace_data']= tmp / 100
        else:
            # Read reference level
            # For linear the measurment units are between 0 and 10000
            # To convert theme to the real values, the measurment units has to be mapped to the reference level
            reference_level = float(osa.query('RL?').rstrip('\n'))
            data_dict['Trace_data'] = tmp / 10000 * reference_level

        # Write unit infromation to data dict
        data_dict['Unit'] = amplitude_unit

        # Get sensitivity
        # Page 7-438
        data_dict['Sensitivity'] = float(osa.query('SENS?').rstrip('\n'))
        
        # Get resolution bandwidth
        # Page 7-405
        data_dict['Resolution_BW'] = float(osa.query('RB?').rstrip('\n'))
        
        # Write wavelength informations to data_dict
        data_dict['Start_WL'] = start_wl
        data_dict['Stop_WL'] = stop_wl

        # Create wavelength vector
        data_dict['WL_Vector'] = np.linspace(start_wl, stop_wl, data_dict['Trace_data'].shape[0])

        trace_information[trace]=data_dict
        

    # closing OSA connection
    osa.close()
   
    # closing resource manager 
    rm.close()  

    return trace_information


####### HP_8153A lightwave mulitimeter ##############


def get_samples_HP_8153A_lightwave_multimiter (channels = ['1'], GPIB_address='22',power_units = 'DBM', wavelengths=[1500] ,   log_mode = False):

    """

    get_samples_HP_8153A_lightwave_multimiter
    
    Function for reading samples from a HP_8153A_lightwave_multimiter
    
    Parameters
    ----------
        channels: list of strings, optional (default = ['1'])
            Insert here the wanted channel from the lightwave multimiter as a list of strings. 
            The 2 channels of the lightwave multimiter is channel '1' and '2'.channel one is Hp 81533A 
            and channel 2 is Hp 81531A .

        GPIB_address : string, optional (default = '22')
            The address GPIB address of the lightwave_multimiter .
        
        log_mode: boolean, optional (default = False)
            Enables a log file for this method.
            
        power_units : string, optional (default = 'DBM') . 
            Power units  are named with 'DBM' or 'Watt' 

        wavelength : list of floats optional ( default = [1500]).
            The wavelenghts arguments has maximal 2 .
            length of wavelength and channels length must be same . There are 2 channel ( HP 81533A : The range of wavelength
            for Hp 81533A that must be between 450 nm to 1020 nm) and (HP 81531A : the range of wavelength for
            Hp 81531A that must be between 800  nm to 1700 nm) .

    Returns
    -------
        channel_information: dict
            Consist of dicts which contains the acquired channels data,wavelength,power_unit.
            To access the dict use:
            >Name of object<[>Name of channels<][>Name of data<]
                -> Name of channels: 
                    -> 1 : channel1
                    -> 2 : channel2
                -> Name of data:
                    -> power        : (float) contains powerlevel of the channel 
                    -> Unit         : (string) contains power unit types that must be ( 'DBM' , 'watt')
                    -> wavelength   : (float)contains The wavelenghts arguments in nanometers (NM) 

    Errors
    -------
        Type Error: 
            Will be raised when a wrong data type is used for the input parameter
            -> Possible errors
               -> Type of channels must be list .
               -> Type of channel items must be string. 
               -> Type of GPIB_address must be string.
               -> Type of power units must be string.
               -> Type of wawelength must be float.
               -> Type of wawelength must be list.

        Value Error:
            Will be raised when the input parameter is in an wrong range
            -> Possible errors
              ->  Too many channels. The lightwave mulitimeter has maximal 2 channels.
              ->  Too less channels. Use at least one channel.
              ->  Wrong channels naming. Channels are named with 1 , 2.
              ->  Wrong power units type naming. power units  are named with DBM or Watt .
              ->  Too less wavelenght argument. The wavelenght arguments must be at least  1.
              ->  Too much wavelenght arguments. The wavelenghts arguments has maximal 2 
              ->  length of wavelength and channels length must be same.
              ->  The renge of wavelength is not correct for Hp 81533A that must be between 450 nm to 1020 nm.
              ->  The renge of wavelength is not correct for Hp 81531A that must be between 800  nm to 1700 nm.

        Exception:
            Will be raised by diverse errors
            -> Possible errors
                -> No connection to the multimiter
    """

    # TODO: test cases 
    
    

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

        if not isinstance(GPIB_address, str):
            raise TypeError('Type of GPIB_address must be string')

        if not all(isinstance(x, str) for x in channels):
            raise TypeError('Type of channel items must be string')
        
        if not isinstance(power_units, str) :
            raise TypeError('Type of power units must be string')

        if not all(isinstance(x,float) for x in wavelengths):
            raise TypeError('Type of wawelength must be float ')

        if not isinstance(wavelengths, list):
            raise TypeError( 'Tyoe of wawelength must be list ')

        if len(channels) > 2:
            raise ValueError('Too many channels ({0}). The lightwave mulitimeter has maximal 2 channels'.format(len(channels)))

        if len(channels) < 1:
            raise ValueError('Too less channels ({0}). Use at least one channel'.format(len(channels)))

        if any((channel_name not in ['1','2']) for channel_name in channels):
            raise ValueError('Wrong channels naming. Channels are named with 1 , 2. ')

        if any((power_unit not in ['DBM','Watt']) for power_unit in power_units):
            raise ValueError('Wrong power units type naming. power units  are named with DBM or Watt ')
            
        if len (wavelengths) <1:
            raise ValueError('Too less wavelenght arguments ({0}). The wavelenght arguments must be at least  1'.format(len(wavelengths)))
        
        if len (wavelengths) >2 :
            raise ValueError('Too much wavelenght arguments ({0}). The wavelenghts arguments has maximal 2 '.format(len(wavelengths)))

        if len (wavelengths) != len (channels) :
            raise ValueError ('length of wavelength and channels length must be same ')
        
        for ch_idx , ch in enumerate (channels):
            if ch==1 :
                if wavelengths[ch_idx] <450 or wavelengths[ch_idx] >1020 :
                    raise ValueError ('The renge of wavelength is not correct for Hp 81533A that must be between 450 nm to 1020 nm')
            else :
                 if wavelengths[ch_idx] <800 or wavelengths[ch_idx] >1700 :
                    raise ValueError ('The renge of wavelength is not correct for Hp 81531A that must be between 800  nm to 1700 nm') 
        



        
    except Exception as e:
        logger.error('{0}'.format(e))
        return sys.exit(0)

    # =============================================================================
    #  importing visa for communication with the lightwave_multimeter 
    # ============================================================================= 

    rm = visa.ResourceManager()

    # open connection to AWG
    logger.info("Create GPIB connection with " + str(GPIB_address))
    try:
        lwm= rm.open_resource('GPIB0::' + GPIB_address + '::INSTR')
    except Exception as e:
        logger.error('No connection possible. Check GPIB connection \n  {0}'.format(e))
        return sys.exit()
    

    # =============================================================================
    #  Settings for the analyzer
    # ============================================================================= 
    
    # check the unit 
    #page (8-21)
    #power_unit = lwm.query('sense:power:unit?')
    #print(power_unit)
       # Create dict with the the keys 
    channel_information = dict.fromkeys(['channel_'+channels[0],'channel_'+channels[1]])

    for channel,wavelength,power_unit in zip(channels,wavelengths,power_units):
        #This command sets the units in use when an absolute reading is made. This can be dBm (DBM)(0) or Watts (Watt)(1).
        #Page 8-21
        lwm.write('sense{0:s}:power:unit {1:s}'.format(channel,power_unit))
        #check the wavelength 
        old_wavelength = float(lwm.query('sense{0:s}:power:wavelength?'.format(channel)))

        #set new wavelength 
        #nanometers (NM) , micrometers(um), meters (M)
        if wavelength != old_wavelength:
            lwm.write('sense{0:s}:pow:wave {1:d}NM'.format(channel,wavelength))
        
        # we need to get the values
        #page 8-8 , 8-9
        channel_power_level = lwm.query('read{0:s}:power?'.format(channel,))

        #make dictionary for power unit and wavelength 
        data_dict={'power':channel_power_level , 'unit':power_unit , 'wavelength':wavelength }

        #write the data in the dictionary 
        channel_information["channel_{0:s}".format(channel)]=data_dict
    

    # closing lwm connection
    lwm.close()
   
    # closing resource manager 
    rm.close()  

    return channel_information

    



     


