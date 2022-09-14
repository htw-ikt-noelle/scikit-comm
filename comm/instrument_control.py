import logging, struct, sys, time
import pyvisa as visa
import numpy as np
from socket import socket, AF_INET, SOCK_STREAM, SHUT_RDWR


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

def get_screenshot_DLM2034(address='192.168.1.11', folder='.', fname=None, timestamp=True, 
                           tone='COLOR', capture_info='OFF'):
    """
    Get and save screenshot from Yokogawa DLM 2034/3024 sscilloscope as PNG.
    

    Parameters
    ----------
    address : string, optional
        IP adress of the device. The default is '192.168.1.11'.
    folder : string, optional
        Folder to save the screenshot to. The default is '.'.
    fname : string, optional
        Filename to save the screenshot to. If None the standard name 'osci_screenshot'
        is used. The default is None.
    timestamp : bool, optional
        Should a timestamp be added to the filename. The default is True.
    tone : string, optional
        Tone of the output image. Allowed are 'COLOR', 'GRAY', 'OFF', 'REVERSE'. 
        The default is 'COLOR'.
    capture_info : string, optional
        Should oscilloscope settings (e.g. trigger mode) be included in the image?
        'ON' for yes and 'OFF' for no. The default is 'OFF'.    

    """
    
    if fname is None:
        fname = 'osci_screenshot'
    
    if timestamp:
        fname = time.strftime('%Y-%m-%dT%H%M%S_') + fname
        
    fname = folder + '/' + fname
        
    # create resource 
    rm = visa.ResourceManager('@py')    

    # open connection to scope
    scope = rm.open_resource('TCPIP::' + address + '::INSTR', timeout=5000)
    # set number of bytes to retireve at once...TODO: find reasonable value
    scope.chunk_size = 2000

    # check instrument IDN
    idn = scope.query('*IDN?')
    print(idn)
    
    # set type to PNG
    scope.write(':IMAGE:FORMAT PNG')
    # set capture output mode: HARD, NORMAL, WIDE
    scope.write(':IMAGE:MODE HARD')
    # show setting information (e.g. trigger mode, etc.) ON OFF
    scope.write(':IMAGE:INFORMATION ' + capture_info)
    # options are COLOR, GRAY, OFF, REVERSE
    scope.write(':IMAGE:TONE ' + tone)
    # get image data
    scope.write(':IMAGE:SEND?')     
    img_data = scope.read_binary_values(datatype='b', is_big_endian=False, container=np.array)
    
    with open(fname + '.png','wb') as pic_file:
        pic_file.write(img_data)
        pic_file.close()
    
    # close connection and delete objects
    rm.close()
    del rm
    del scope    


def write_samples_Agilent_33522A(samples, ip_address='192.168.1.44', sample_rate=[250e6], offset=[0.0], amp_pp=[1.0], channels=[1], out_filter=['normal'], wait_for_ext_trigger=[False], trig_delay=[0.0]):
    """
    write_samples_Agilent_33522A
    
    Function for writing samples to an Agilent/Keysight 33500 Series 30MHz Function/Arbitrary Waveform Generator

    Parameters
    ----------
    samples : numpy array, n_outputs x n_samples , float
        samples to output, to be scaled between -1 and 1 (values outside this range are clipped).
    ip_address : string, optional
        The default is '192.168.1.44'. Currently, only LAN connection is supported.
    
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
    wait_for_ext_trigger: list of bools
        Should the device wait for an external trigger (raising edge at the device 
        backpanel "Ext trig") to start the ouput? The default is [False].
    trig_delay: list of floats
        The signal is output trig_delay seconds after the trigger event. This parameter
        takes only effect if wait_for_eyt_trigger is True. The default is 0.0.  
        

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
              
        # applying sample rate, amplitude and Offset        
        awg.write(':SOUR{0:d}:APPL:ARB {1:g},{2:g}, {3:g}'.format(ch, sample_rate[ch_idx], amp_pp[ch_idx], offset[ch_idx]))
        #awg.write(':SOURce%d:APPLy:ARBitrary %s,%s,%s' % (ch, sample_rate[ch_idx], amp_pp[ch_idx], offset[ch_idx]))
        
        # applying output filter mode
        awg.write(':SOUR{0:d}:FUNC:SHAP:ARB:FILT {1:s}'.format(ch, out_filter[ch_idx].upper()))
        
        # wait a moment to have the output to turned on
        time.sleep(0.1)
        
        if wait_for_ext_trigger:
            # set external trigger and delay
            awg.write('TRIG{:d}:DELAY {:.1e}'.format(ch, trig_delay[ch_idx]))
            awg.write('TRIG{:d}:SOURCE EXTERNAL'.format(ch))
            # set burst mode (device triggers only in burst or wobble mode)
            # device waits for trigger event and outputs the waveform continously
            # --> only single trigger event needed, further trigger events have no effect
            awg.write('SOUR{:d}:BURST:MODE TRIGGERED'.format(ch))
            awg.write('SOUR{:d}:BURST:NCYCLES INFINITY'.format(ch))
            awg.write('SOUR{:d}:BURST:PHASE 0.0'.format(ch))
            awg.write('SOUR{:d}:BURST:STATE ON'.format(ch))
            
               
    awg.write(':SOUR{0:d}:FUNC:ARB:SYNC'.format(ch))  # synchronising channels
        
    awg.close() # closing AWG
    rm.close()  # closing resource manager 
    
    
def write_samples_TTI_TG5012A(samples=np.asarray([]), ip_address='192.168.1.105', 
                              waveform='SINE', amp_pp=1.0, channel=1, 
                              repetition_freq=1000.0, memory='ARB1', interpolate='OFF', 
                              bit_rate=1000.0, mute_output=False):
    """
    Write samples to TTI arbitrary waveform generator and set waveform parameters.
    
    This funtion allows for setting an output waveform of a TTI TG5012A function and 
    arbitrary waveform generator (AWG). Customized waveforms can be uploaded too.
    
    Parameters
    ----------
    samples : numpy array, float
        samples to be output, have to be scaled between -1 and 1 (values outside this 
        range are clipped). The length of the samples array must be between 2 and 131072.
        The default is np.asarray([]).
    ip_address : string, optional
        IP address of AWG. The default is '192.168.1.105'.
    waveform : string, optional
        What waveform to output. Allowed are 'SINE', 'SQUARE', 'RAMP', 'TRIANG', 
        'PULSE', 'NOISE', 'PRBSPNX', 'ARB', 'DC'. The default is 'SINE'. PRBS lengths
        can be 7, 9, 11, 15, 20 or 23 (e.g. PRBSPN7).
    amp_pp : float, optional
        Peak-to-peak voltage of the output signal in Volts. Be careful: the peak-to-peak
        amplitude has different meanings for different waveform types. For example: 
        in case of 'DC' it is the amplitude of the DC, while for custom waveforms it
        is the voltage difference of the smallest and largest uploaded sample.
        Note that the actual output voltage is also influenced by the setting of the 
        source impedance and the connected load impedance!
    channel : int, optional
        Output channel that is programmed. The default is 1.
    mute_output : boolean, optional
        Mute (switch OFF) the selected output during programming to prevent glitches.
        Note: The output of the selected channel will always be switched to ON 
        after programming.  The default is False.
    repetition_freq : float, optional
        Repetition frequency of the custom waveform in memory in Hz. This 
        is the read-out repetition frequency of the whole memory and therefore, 
        multiplied with the length of the custom waveform, determines the sample rate.
        For other periodic waveforms (e.g. SINE), the parameter specifies the 
        signal frequency. The default is 1000.0.
    memory : string, optional
        Internal memory identifier used to store the custom waveform. Allowed
        names are 'ARB1', 'ARB2', 'ARB3' and 'ARB4' . The default is 'ARB1'.
    interpolate : string, optional
        The AWG is only able to output waveforms of length 2**14 or 2**17. Uploaded
        waveforms shorter than 2**14 are extended to 2**14 and signals of length 
        between 2**14 and 2**17 are extended to 2**17. If interpolate is 'ON', 
        the missing samples are linearly interpolated. If interpolate is 'OFF',
        the missing samples are generated by repeating samples (see p.79 of data 
        sheet of the AWG).The default is 'OFF'.
    bit_rate : float, optional
        Bit rate in case of a PRBS signal in bits/s. The default is 1000.0.    

    Returns
    -------
    None.

    """
    sleeptime = 0.1  # [sec] pause-time between VISA IO operations (increase if driver is unreliable)
    
    if not isinstance(samples, np.ndarray):
        raise TypeError('samples has to be a numpy array...')
        
    if waveform.upper() == 'ARB':
        if (samples.size < 2) or (samples.size > 2**17):
            raise ValueError('length of waveform must be between 2 and 131072 points...')
            
        # clip samples to +-1.0
        samples = np.clip(samples, -1.0, 1.0)
        # shift to positive values and scale to range 0.0...1.0
        samples = (samples + 1.0) / 2.0
        # convert to int between 0...16383
        samples = np.round(samples * 16383.0)
        
    # create resource 
    rm = visa.ResourceManager('@py')
    
    # open connection to AWG
    awg = rm.open_resource('TCPIP::' + ip_address + '::9221::SOCKET', 
                           read_termination='\n', write_termination='\n', timeout=3000)
    
    # select instrument channel for subsequent commands 
    awg.write('CHN {:d}'.format(channel))
    time.sleep(sleeptime)
    
    # switch output X OFF
    if mute_output==True: 
        awg.write('OUTPUT OFF');  
        time.sleep(sleeptime)
    
    # set waveform type
    if waveform.upper() == 'DC':
        awg.write('WAVE ARB');   time.sleep(sleeptime)
        awg.write('ARBLOAD DC'); time.sleep(sleeptime)
        awg.write('ARBDCOFFS {:.6f}'.format(amp_pp))
    else:
        if ('SINE' or 'SQUARE') in waveform.upper(): # workaround, since 'SINE' and 'SQUARE 'sometimes don't load correctly
            awg.write('WAVE TRIANG');
            time.sleep(sleeptime)
        awg.write('WAVE ' + str(waveform.upper())); 
        time.sleep(sleeptime)
    
    if waveform.upper() == 'ARB':
        # load waveform in memory 'ARBX'
        # this is done 'by hand' instead of using "awg.write_binary_values", because
        # the AWG is not compatible with the binary write method of pyvisa...WHY???
        # does not work:
        # awg.write_binary_values('{:s}'.format(memory), samples, datatype='h', header_fmt='ieee', is_big_endian=True)

        # generate IEEE header (see p. 100) of AWG datasheet
        # number of bytes to send (int16 values)
        bytes = samples.size*2
        # number of digits in header
        digits = np.floor(np.log10(bytes))+1
        header = '{:s} #{:d}{:d}'.format(str(memory), int(digits), int(bytes))
        # send header without any termination character
        awg.write('{:s}'.format(header), termination='')
        time.sleep(sleeptime)

        # convert 16bit int ('h') to bytes in big endian order ('>')
        awg.write_raw(samples.astype('>h').tobytes())
        time.sleep(sleeptime*4)
        awg.write('ARBDEF {:s},{:s},{:s}'.format(memory,memory+'u',interpolate))
        time.sleep(sleeptime*3)
        awg.write('ARBLOAD {:s}'.format(memory))
        time.sleep(sleeptime*3)

    
    if not waveform == 'DC':
        # set signal speed
        if 'PRBSPN' in waveform.upper():
            awg.write('PRBSBITRATE {:.6f}'.format(bit_rate))
        else:
            awg.write('FREQ {:.6f}'.format(repetition_freq))
        time.sleep(sleeptime)
        # set ptp amplitude
        awg.write('AMPUNIT VPP')
        time.sleep(sleeptime)
        awg.write('AMPL {:.6f}'.format(np.abs(amp_pp)))
        # align output(s)
        awg.write('ALIGN')        
    time.sleep(sleeptime)    
    
    # set offset to 0 (unknown command)
    # awg.write('OFFSET 0.0')

    # switch output on
    awg.write('OUTPUT ON')
    time.sleep(sleeptime)       
    
    # Returns the instrument to local operation and unlocks the keyboard
    awg.write('LOCAL') 
    
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
    
    
    # Create standard output handler
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.ERROR)
    
    # Set format of the logs with formatter
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(name)s :: Line No %(lineno)d:: %(message)s')
    
    # Adding formatter to handler
    stdout_handler.setFormatter(formatter)
    
    # Adding handler to logger
    logger.addHandler(stdout_handler)
    
    
    if log_mode == True:
        # Create file handler 
        file_handler = logging.FileHandler('{0}.log'.format(__name__))
        file_handler.setLevel(logging.INFO)
        
        # Adding formatter to handler
        file_handler.setFormatter(formatter)
        
        # Adding handler to logger
        logger.addHandler(file_handler)

    # =============================================================================
    #  Check inputs for correctness
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

    # Get sample rate of the scope and convert from string to float
    sample_rate = float(scope.query('HORizontal:MODE:SAMPLERate?'))

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

        # Reading vertical scaling of the scope (Voltage per div)
        ver_scale = float(scope.query('CH{0:d}:SCAle?'.format(ch)))

        # Reading vertical position ( Y-Position on scope screen)
        ver_position = float(scope.query('CH{0:d}:POSition?'.format(ch)))

        # Reading offset from scope
        offset = float(scope.query('CH{0:d}:OFFSet?'.format(ch)))

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
    
    # limitations of specific AWG device    
    MAX_SR = 16e9
    MIN_SR = 1.49e3
    
    
    # =============================================================================
    #  Create logger which writes to file
    # ============================================================================= 
    # Create logger
    logger = logging.getLogger(__name__)
    
    # Set the log level
    logger.setLevel(logging.INFO)
    
    
    # Create standard output handler
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.ERROR)
    
    # Set format of the logs with formatter
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(name)s :: Line No %(lineno)d:: %(message)s')
    
    # Adding formatter to handler
    stdout_handler.setFormatter(formatter)
    
    # Adding handler to logger
    logger.addHandler(stdout_handler)
    
    
    if log_mode == True:
        # Create file handler 
        file_handler = logging.FileHandler('{0}.log'.format(__name__))
        file_handler.setLevel(logging.INFO)
        
        # Adding formatter to handler
        file_handler.setFormatter(formatter)
        
        # Adding handler to logger
        logger.addHandler(file_handler)



    # =============================================================================
    #  Check inputs for correctness
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

        if sample_rate[0] > MAX_SR or sample_rate[0] < MIN_SR:
            raise ValueError('Sample rate not supported by specific AWG device')

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
        # print(awg.query('WLIST:WAVEFORM:DATA? "Python_waveform_AWG_{0:d}",0'))

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
    
def get_spectrum_IDOSA(ip_address='192.168.1.22', new_sweep = False, wl_equidist = False):
    """
    get_spectrum_IDOSA()
    
    Function for reading the optical spectrum from an ID-Photonics ID-OSA.
    OSA settings have to be adjusted externally (using e.g., the ID-OSA GUI) 
    before reading th spectrum (for details, see Manual_IDOSA.pdf).

    Parameters
    ----------
    ip_address : string, optional (default: '192.168.1.22')
        Th IP-address of the device.
        
    new_sweep : boolean or int, optional  (default: False)
        If FALSE or 0, the current spectrum will be fetched from the instrument. 
        If TRUE or 1, a new OSA sweep is initiated and the spectrum is fetched thereafter.
        The device is then reset to the original sweep mode.

    wl_equidist: boolean or int, optional (default: False)
        If FALSE or 0, the returned wavelengths are calculated from an equidistant frequency
        vector which is always used internally by the device.
        If TRUE or 1, the spectrum is interpolated at equidistant wavelengths, where the
        wavelength is interpolated using the same range and number as the original wavelength axis.

    Returns
    -------
    trace: dict
        The dictionary holds the spectrum measurement results in the following key-value pairs:
            key                     value
            'Resolution_BW_Hz'      (float) The OSA resolution bandwidth in units of [Hz].
            'Resolution_BW_m'       (float) The resolution bandwidth in units of [m], referenced to center of the spectrum.
            'frequency'        (np.array) Contains the frequency axis (descending) in units of [Hz].
            'wavelength'           (np.array) Contains the wavelength axis (ascending) in units of [m].
            'spectrum_dBm'            (np.array) Contains the spectrum in log-domain in units of [dBm].
            'Ptotal_dBm_IDOSA'      (float)    The total optical power in [dBm], interally calculated by the instrument.
            'Ptotal_dBm_int'        (float)    The total optical power in [dBm], calculated from the fetched spectrum.
    """
    
    # =============================================================================
    #  Check inputs for correctness
    # ============================================================================= 

    try:
        if not isinstance(ip_address, str):
            raise TypeError('Type of ip_address must be string')
            
        if not isinstance(new_sweep, (bool,int)):
            raise TypeError('Type of new_sweep must be boolean or int')

        if not isinstance(wl_equidist, (bool,int)):
            raise TypeError('Type of wl_equidist must be boolean or int')
    except Exception as e:
        print('{0}'.format(e))
        sys.exit(0)

    
    sleeptime = 0.03  # [sec] pause-time between socket IO operations (increase if driver is unreliable)
    RCV_BUFFSIZE = 32 # size for receive-buffer
    socket_timeout = 2 # [sec] socket timeout (increase if desired)
    tcp_port = 2000  # do not change (fix for ID-OSA)
    c0 = 299792458.0 # [m/s] speed of light
    
    ## Create dictionary for trace data and wavelength information
    trace = {'Resolution_BW_m':np.nan, 'Resolution_BW_Hz':np.nan, 'Ptotal_dBm_IDOSA':np.nan,
             'Ptotal_dBm_int':np.nan, 'wavelength':np.asarray(np.nan),
             'frequency':np.asarray(np.nan), 'spectrum_dBm':np.asarray(np.nan)}
    
    ## connect to socket
    osa = socket(AF_INET,SOCK_STREAM) # https://docs.python.org/3/library/socket.html#socket.socket.connect
    osa.settimeout(socket_timeout) # 2 sec timeout
    try:
        osa.connect( (ip_address,tcp_port) )
    except Exception as e:
        print(e)
        return trace
        sys.exit("Could not connect to socket.")
    
    with osa:
        if bool(new_sweep) == True:
            ## remember current sweep mode
            osa.sendall('smod?;'.encode())
            time.sleep(sleeptime/2);
            sweep_mode = int(osa.recv(RCV_BUFFSIZE).decode().lstrip(';\r\n').rstrip(';\r\n'))
            
            if (sweep_mode==2 or sweep_mode==3): # if repeat or auto -> set to single sweep mode (smod 1)
                osa.sendall('smod 1;'.encode())
                time.sleep(sleeptime); dummy = osa.recv(RCV_BUFFSIZE) # dummy read to empty send buffer (remove \r\n etc.)
            
            ## initiate a single sweep
            osa.sendall('SGL;'.encode());  dummy = osa.recv(RCV_BUFFSIZE)
            #osa.sendall('*WAI;'.encode()); dummy = osa.recv(RCV_BUFFSIZE) # TODO: query OPC? instead of *WAI
            OPC = False
            while not(OPC): # wait for sweep complete
                osa.sendall('*OPC?;'.encode()); time.sleep(sleeptime);
                OPC = bool(int(osa.recv(RCV_BUFFSIZE).decode().lstrip(';\r\n').rstrip(';\r\n')))
            
        ## query the resolution bandwidth (RBW in Hz)
        osa.sendall('step:freq?;'.encode()); time.sleep(sleeptime*2);
        RBW = float(osa.recv(RCV_BUFFSIZE).decode().lstrip(';\r\n').rstrip(';\r\n')); # OSA resolution bandwidth in [Hz]
        #print('ID-OSA RBW: {:3.2f} MHz'.format(RBW/1e6))

        ## query the center frequency (Hz) 
        osa.sendall('cent?;'.encode()); time.sleep(sleeptime*2);
        f_cent = float(osa.recv(RCV_BUFFSIZE).decode().lstrip(';\r\n').rstrip(';\r\n')); # spectrum center in [Hz]
        #print('center frequency (by ID-OSA): {:3.2f} Hz'.format(f_cent))
        
        ## query total optical power from intrument (migth deviate from feteched spectrum if OSA is in RPT scan mode)
        osa.sendall('POW?;'.encode()); time.sleep(sleeptime*2);
        Power_dBm = float(osa.recv(RCV_BUFFSIZE).decode().lstrip(';\r\n').rstrip(';\r\n'));
        #print('Total optical power (by ID-OSA): {:3.3f} dBm'.format(Power_dBm))
        
        ## query OSA trace
        # fetch wavelength axis in ascending order in units of [m] (precedure adapted from LabVIEW ID-OSA driver)
        osa.sendall('FORM REAL,64;'.encode()); time.sleep(sleeptime); dummy = osa.recv(RCV_BUFFSIZE) # double is required for full resolution
        osa.sendall('x?;'.encode()); time.sleep(sleeptime); # wavelength in [m];
        digits = int(osa.recv(2).decode().lstrip(';\r\n').rstrip(';\r\n')[1]); time.sleep(sleeptime) # no. of digits to read
        num_bytes = int(osa.recv(digits).decode().lstrip(';\r\n').rstrip(';\r\n')); time.sleep(sleeptime); # no. of bytes to read
        x = osa.recv(num_bytes+3); # binary data block [{+3 bytes: needed to capture term sequence}
        databytes = bytes(list(zip(*list(struct.iter_unpack("B",x[0:-3]))))[0][::-1]) # [0:-3]: strip off term sequence | unpack and reverse | cast to bytes
        WL_m = np.asarray(list(struct.iter_unpack(">d",databytes))).flatten()[::-1] #">d" for big-endian (8-byte) double (https://docs.python.org/3/library/struct.html)
        f_Hz = c0/WL_m # EQUIDISTANT frequency axis
        
        # fetch y-data (power spectrum) in units of [dBm]
        osa.sendall(':TRACE:DATA:LINLOG LOG;'.encode()); time.sleep(sleeptime); dummy = osa.recv(RCV_BUFFSIZE*3); # set to dBm
        osa.sendall('FORM REAL,32;'.encode()); time.sleep(sleeptime); dummy = osa.recv(RCV_BUFFSIZE) # 32-bit float sufficient for spectrum in dB 
        osa.sendall('y?;'.encode()); time.sleep(sleeptime); # power spectrum in units of (dBm in RBW)
        digits = int(osa.recv(2).decode().lstrip(';\r\n').rstrip(';\r\n')[1]); time.sleep(sleeptime) # no. of digits to read
        num_bytes = int(osa.recv(digits).decode().lstrip(';\r\n').rstrip(';\r\n')); time.sleep(sleeptime); # no. of bytes to read
        y = osa.recv(num_bytes+3); # binary data block [{+3 bytes: needed to capture term sequence}
        sb = bytes(list(zip(*list(struct.iter_unpack("B",y[0:-3]))))[0][::-1]) # y[0:-3]: strip off term sequence
        Spec_dBm = np.asarray(list(struct.iter_unpack(">f",sb))).flatten()[::-1] #">f" for big-endian (4-byte) float

        ## reset OSA to original sweep mode
        if new_sweep == True and (sweep_mode==2 or sweep_mode==3):
            osa.sendall('RPT;'.encode()); time.sleep(sleeptime); dummy = osa.recv(RCV_BUFFSIZE)
            
        ## close socket connection
        osa.shutdown(SHUT_RDWR)
        osa.close()

        ## resample spectrum to equidistant wavelength steps (experimental)
        if bool(wl_equidist) == True:
            WL_m_i = np.linspace(np.min(WL_m),np.max(WL_m),WL_m.size,endpoint=True) # [m]
            Spec_dBm = np.interp(WL_m_i, WL_m, Spec_dBm) #[dBm]
            WL_m = WL_m_i; del WL_m_i
            f_Hz = c0 / WL_m
            
        ## calculate total power from spectrum        
        Pwr_integrate_dBm = 10*np.log10(np.abs(np.trapz(10**(Spec_dBm[::-1]/10-3),x=f_Hz[::-1])/RBW/1e-3))
        #print('Total optical power (spectrum integration): {:3.3f} dBm'.format(Pwr_integrate_dBm))    
    
        
        trace['Resolution_BW_m'] = np.abs(-c0 / f_cent**2 * RBW) # RBW in [m], referenced to center of spectrum
        trace['Resolution_BW_Hz'] = RBW # [Hz]
        trace['wavelength'] = WL_m # [m]
        trace['frequency'] = f_Hz # [Hz]
        trace['spectrum_dBm'] = Spec_dBm # [dBm]
        trace['Ptotal_dBm_IDOSA'] = Power_dBm # [dBm]
        trace['Ptotal_dBm_int'] = Pwr_integrate_dBm # [dBm]
        return trace    


def get_spectrum_HP_71450B_OSA (traces = ['A'], GPIB_bus=0, GPIB_address=13,log_mode = False, single_sweep = False):


    """
    get_spectrum_HP_71450B_OSA
    
    Function for reading spectrum from a HP_71450B optical spectrum analyzer
    
    Parameters
    ----------
        traces: list of stings, optional (default = ['A'])
            Insert here the wanted traces from the OSA as a list of strings. 
            The three traces of the OSA are A, B, and C. It is also possible to use lower case.
        
        GPIB_bus : int, optional (default = 0)
            The GPIB bus number of the OSA.

        GPIB_address : int, optional (default = 13)
            The GPIB address of the OSA.
        
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
                    -> spectrum   : (np.array) Contains numpy array with trace data
                    -> Unit         : (string) Contains the unit of the trace data
                    -> Sensitivity  : (float) Contains the amplitude sensitivity of the spectrum (always in dBm)
                    -> Start_WL     : (float) Contains the start wavelength of the spectrum (in m)
                    -> Stop_WL      : (float) Contains the stop wavelength of the spectrum (in m)
                    -> Resolution_BW: (float) Contains the resolution bandwidth of the spectrum (in m)
                    -> wavelength    : (np.array) Numpy array with evenly spaced wavelengths between Start_WL and Stop_WL (in m)
            
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
    # TODO: set GBIB Bus number
    
    
    # =============================================================================
    #  Create logger which writes to file
    # ============================================================================= 
    # Create logger
    logger = logging.getLogger(__name__)
    
    # Set the log level
    logger.setLevel(logging.INFO)
    
    
    # Create standard output handler
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.ERROR)
    
    # Set format of the logs with formatter
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(name)s :: Line No %(lineno)d:: %(message)s')
    
    # Adding formatter to handler
    stdout_handler.setFormatter(formatter)
    
    # Adding handler to logger
    logger.addHandler(stdout_handler)
    
    
    if log_mode == True:
        # Create file handler 
        file_handler = logging.FileHandler('{0}.log'.format(__name__))
        file_handler.setLevel(logging.INFO)
        
        # Adding formatter to handler
        file_handler.setFormatter(formatter)
        
        # Adding handler to logger
        logger.addHandler(file_handler)

    # =============================================================================
    #  Check inputs for correctness
    # ============================================================================= 

    try:
        if not isinstance(traces, list):
            raise TypeError('Type of traces must be list')
            
        if not isinstance(GPIB_bus, int):
            raise TypeError('Type of GPIB_bus must be int')

        if not isinstance(GPIB_address, int):
            raise TypeError('Type of GPIB_address must be int')

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
        osa = rm.open_resource('GPIB' + str(GPIB_bus) +'::' + str(GPIB_address) + '::INSTR')
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
    start_wl = float(osa.query('STARTWL?').rstrip('\n') )

    # Read stop wave length
    # Page 7-464 -> 7-465
    # Convert from m to nm Page 1-14
    stop_wl = float(osa.query('STOPWL?').rstrip('\n') )

    # Loop through traces
    for trace_id,trace in enumerate(traces):

        # Create dictionary for trace data and wave length information
        data_dict = {'spectrum':[],'Unit':[],'Sensitivity':[],'Start_WL':[],'Stop_WL':[],'Resolution_BW':[], 'wavelength':[]}

        # Setting length of Trace
        # Page 7-506 -> 7-507
        trace_length = int(float(osa.query('TRDEF TR{0:s}?'.format(trace)).rstrip('\n')))

        # Read trace
        # Page 7-499 -> 7-502
        # h is 2 bytes (signed short)
        tmp = osa.query_binary_values('TR{0:s}?'.format(trace),datatype='h' ,is_big_endian=True, container=np.array,data_points = trace_length)

        # Convert measument units to parameter units
        # Page 2-8
        if is_log:
            # One measurement unit is equal to one hundreth of a dBm
            # To get the dBm the trace data from the scope has to be divided by 100
            data_dict['spectrum']= tmp / 100
        else:
            # Read reference level
            # For linear the measurment units are between 0 and 10000
            # To convert theme to the real values, the measurment units has to be mapped to the reference level
            reference_level = float(osa.query('RL?').rstrip('\n'))
            data_dict['spectrum'] = tmp / 10000 * reference_level

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
        data_dict['wavelength'] = np.linspace(start_wl, stop_wl, data_dict['spectrum'].shape[0])

        trace_information[trace]=data_dict
        

    # closing OSA connection
    osa.close()
   
    # closing resource manager 
    rm.close()  

    return trace_information

def get_opt_power_Anritsu_ML910B(GPIB_bus=0, GPIB_address=11):
    """
    Read optical power from Anritsu ML910B optical power meter.
    
    This driver does not set any parameters of the instrument. Make sure that all required
    paremeters e.g. power unit, range, wavelength and selected channel are correctly
    set manually via the device panel.The driver simply queries the power values
    of the selected channels.

    Parameters
    ----------
    GPIB_bus : int, optional
        GPIB bus number to which the instrument is connected. The default is 0.
    GPIB_address : int, optional
        GPIB address set to the instrument. The default is 11.
    

    Returns
    -------
    output : dict of dicts
        A dict containing the keys 'ch1' and 'ch2', respectively. Each entry 
        consists of a dict with keys 'value' and 'unit' containing the measured
        value and corresponding measurement unit.
    """
    
    
    rm = visa.ResourceManager()

    # open connection to AWG
    pm = rm.open_resource('GPIB' + str(GPIB_bus) + '::' + str(GPIB_address) + '::INSTR')
    
    # set termination characters
    pm.read_termination='\r\n'
    
    output = {}
    
    measure_units = {
        'P': 'mW',
        'Q': 'uW',
        'R': 'nW',
        'S': 'pW',
        'T': 'dBm',
        'U': 'dB',
        'Y': 'rel recall',
        'N': 'DIFF'
        }
    
    channel_ind = {
        'L': 'ch1',
        'M': 'ch2',
        '0': 'ch1',
        '1': 'ch2'
        }
    
    # read info from device
    results = pm.read_bytes(count=30,break_on_termchar=True)
    
    # convert to str and strip whitespace characters
    results = str(results, encoding='ASCII').rstrip('\n').rstrip('\r')
    
    # split
    results = results.split(',')  
    
    # only one channel is selected
    if len(results) == 1:
        # generate new dict for each channel
        output[channel_ind[results[0][5]]] = {}
        # write measurement unit into dict
        output[channel_ind[results[0][5]]]['unit'] = measure_units[results[0][1]]
        # is measurement in range?
        if results[0][0] == 'X':
            # write measurement value into dict
            output[channel_ind[results[0][5]]]['value'] = float(results[0][-6:].replace(' ',''))
        else:
            # out of range
            output[channel_ind[results[0][5]]]['value'] = None
    elif len(results) == 2:
        # write results from both channels
        for idx, result in enumerate(results):
            # new dict
            output[channel_ind[str(idx)]] = {}
            # write unit into dict
            output[channel_ind[str(idx)]]['unit'] = measure_units[result[1]]
            # is measurement in range?
            if result[0] == 'X':
                output[channel_ind[str(idx)]]['value'] = float(result[-6:].replace(' ',''))
            else:
                output[channel_ind[str(idx)]]['value'] = None
                
    else:
        raise TypeError('')
        
    pm.close() # closing AWG
    rm.close()  # closing resource manager 
    
    return output
    
    


####### HP_8153A lightwave multimeter ##############
def get_opt_power_HP8153A(channels, GPIB_bus=0, GPIB_address=22 ,power_units = [None], wavelengths = [None] ,verbose_mode = True ,log_mode = False):
    """

    get_opt_pwr_HP8153A
    
    Function for reading power values from a HP8153A lightwave multimeter.
    
    Wavelength ranges of the currently used modules:
        - HP 81533A  850 nm to 1700 nm 
        - HP 81531A  800 nm to 1700 nm
    
    This list can be expanded for additional modules
    
    Parameters
    ----------
        channels: list of strings
            Insert here the required channel of the lightwave multimeter as a list of strings.
            
            - for channel 1:
            - >>> channels = ['1']
            - for channel 2: 
            - >>> channels = ['2']
            - for channels 1 and 2:
            - >>> channels = ['1','2']
                        
            The channel assignment of the parameters power_units and wavelengths corresponds to the elements of this list.

        GPIB_bus : int
            The GPIB bus number of the lightwave multimeter. Use a value between 1 and 30

        GPIB_address : int
            The GPIB address of the lightwave multimeter. Use a value between 1 and 30
        
        log_mode: boolean, optional (default = False)
            Enables a log file for this method.
            
        power_units : list of strings, optional (default = [None]).
            Sets the power unit(s) for the acquired channel(s).
            Available power units are 'DBM' and 'Watt' (case-sensitive).
            Maximum number of list items is 2.
            
            - Examples for one channel :
            >>> power_units = ['DBM']
            >>> power_units = ['Watt']
            - Examples for two channels:
            >>> power_units = ['DBM','Watt']
            >>> power_units = ['DBM','DBM']
            
            The assignment of the units corresponds to the elements in the channels list.
            When the unit list is changed, its length must be the same as the channels list.
            
            If the unit should not be changed, ignore this parameter or set to 'None'.
            
            Special case: The power level of both channels shall be acquired and for one channel the power unit should
            be changed. In this case, a 'None' is used for the channel that should be not changed.
            >>> power_units = ['DBM','None']
         

        wavelengths : list of floats, optional ( default = [None]).
            Sets the calibration wavelength(s) for the acquired channel(s).
            Maximum number of list items is 2.
            
            - Example one channel:
            >>> wavelengths = ['1550']
            - Example two channels:
            >>> wavelengths = ['1550','1500']
            
            The assignment of the wavelengths corresponds to the elements in the channels list.
            When the wavelength list is changed, its length must be the same as the channels list.
            
            If the wavelength(s) should not be changed, ignore this parameter or set it to 'None'.
            
            Special case: The power level of both channels should be acquired and for one channel the wavelength should 
            be changed. In this case, a 'None' is used for the channel that should not changed.

            - Only channel one should change:
            >>> wavelengths = ['1540','None']
            
            Warning: If the wavelength setting of a channel is out of range, the multimeter will ignore the setting.
                 
        verbose_mode : boolean (default = True)
            When this mode is activated, additional information such as the current wavelength, the power unit and the 
            module name are returned from the function. 
            If only the current power is required, this mode should be deactivated (False) to save unnecessary IO operations.

    Returns
    -------
        channel_information: dict
            Returns two versions depending on verbose_mode.
            
            verbose_mode = False:
                Consist of dicts which only contain the acquired channel power level.
            
            verbose_mode = True:
                Consist of dicts which contain the acquired channel power level, wavelength, power unit and modulename.
                
                - To access the dict:
                >>> object[Name of channels][Name of data]
                
                - Name of channels: 
                    - '1': channel 1
                    - '2': channel 2
                - Name of data:
                    - Power: (float) contains the power level of the channel
                    - Unit: (string) contains the power unit
                    - Wavelength: (float) contains the calibration wavelength in nanometers (nm)
                    - Module: (string) contains the module name for the channel
                        
        To access the dict use:
            >>> object[Name of channels][Name of data]
            
                    - Name of channels: 
                        - '1': channel 1
                        - '2': channel 2
                    - Name of data:
                        - Power: (float) contains power level of the channel 
                                       
    Examples
    --------
            >>> import comm as comm
        
        The power of channel 1 should be acquired. The wavelength will be set to 1550 nm and the power unit to dBm.
        GPIB address will be set to 22. In this example, the verbose mode is activated.
            >>> p = comm.instrument_control.get_opt_power_HP8153A(channels = ['1'], GPIB_address = 22, power_units = ['DBM'], wavelengths = [1550.0])
        
        The power of both channels should be acquired. The wavelength of channel 1 should not be changed. For channel 2, 
        the wavelengths will be set to 1550nm. Power unit for channel 1 should be "Watt" and "dBm" for channel 2.
            >>> p = comm.instrument_control.get_opt_power_HP8153A(channels = ['1','2'], GPIB_address = 22, power_units = ['Watt','DBM'], wavelengths = [None,1550.0])
        
        Only the power value should be acquired from both channels. Therfore, the verbose_mode can be deactivated.
            >>> p = comm.instrument_control.get_opt_power_HP8153A(channels = ['1','2'], GPIB_address = 22, verbose_mode = False)
        Note, in this example, no values for wavelengths and power_unit are provided. Hence, the current values of the multimeter will be used.
        
        Access the power level only of channel 1
            >>> power_level_ch1 = p['1']['Power']
            
        Access the wavelength of channel 2
            >>> wavelength_ch2 = p['2']['Wavelength']

    Errors
    -------
        Type Error:
            This will be raised when a wrong data type is used for the input parameter.
            - Possible errors
               - Type of channels must be list.
               - Type of channel items must be string.
               - Type of GPIB_address must be string.
               - Type of power units must be string.
               - Type of power_unit must be list.
               - Type of wavelength must be float.
               - Type of wavelength must be list.
               - Type of verbose_mode must be bool.
               - Type of log_mod must be boolean.

        Value Error:
            This will be raised when the input parameter is in an wrong range.
            - Possible errors
              -  Too many channels ({0}). The lightwave mulitimeter has a maximum of 2 channels
              -  Too few channels. Use at least one channel.
              -  Wrong channel naming. Channels are named with 1 or 2.
              -  Wrong power units naming. Allowed power units are the strings DBM and Watt.
              -  Too few wavelength arguments. The number of wavelength arguments must be at least 1.
              -  Too many wavelength arguments. A maximum of 2 arguments is permitted.
              -  Lengths of wavelength and channels lists must be equal.

        Exception:
            Will be raised by diverse errors
            - Possible errors
                - No connection to the multimeter.
    """
    # =============================================================================
    #  Create logger which writes to file
    # ============================================================================= 
    # Create logger
    logger = logging.getLogger(__name__)
    
    # Set the log level
    logger.setLevel(logging.INFO)
    
    
    # Create standard output handler
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.ERROR)
    
    # Set format of the logs with formatter
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(name)s :: Line No %(lineno)d:: %(message)s')
    
    # Adding formatter to handler
    stdout_handler.setFormatter(formatter)
    
    # Adding handler to logger
    logger.addHandler(stdout_handler)
    
    
    if log_mode == True:
        # Create file handler 
        file_handler = logging.FileHandler('{0}.log'.format(__name__))
        file_handler.setLevel(logging.INFO)
        
        # Adding formatter to handler
        file_handler.setFormatter(formatter)
        
        # Adding handler to logger
        logger.addHandler(file_handler)


    # =============================================================================
    #  Check inputs for correctness
    # ============================================================================= 

    try:
        if not isinstance(channels, list):
            raise TypeError('Type of channels must be list.')

        if not isinstance(GPIB_bus, int):
            raise TypeError('Type of GPIB_bus must be int.')

        if not isinstance(GPIB_address, int):
            raise TypeError('Type of GPIB_address must be int.')
            
        if not isinstance(verbose_mode, bool):
            raise TypeError('Type of verbose_mode must be bool.')

        if not all(isinstance(x, str) for x in channels):
            raise TypeError('Type of channel items must be string.')
        
        if not all((isinstance(x, str) or x == None) for x in power_units):
            raise TypeError('Type of power_units must be string.')
            
        if not isinstance(power_units, list):
            raise TypeError('Type of power_unit must be list.')

        if not all((isinstance(x,float) or x == None) for x in wavelengths):
            raise TypeError('Type of wavelength must be float.')

        if not isinstance(wavelengths, list):
            raise TypeError( 'Tyoe of wavelength must be list.')
        
        if not isinstance(log_mode, bool):
            raise TypeError('Type of log_mod must be boolean.')
            
        # If no parameters are passed for power_units or wavelengths, the list lengths must be adjusted to the length of the channels list
        if all(x == None for x in power_units):
            power_units = [None]*len(channels)
            
        if all(x == None for x in wavelengths):
            wavelengths = [None]*len(channels)

        if len(channels) > 2:
            raise ValueError('Too many channels ({0}). The lightwave mulitimeter has a maximum of 2 channels'.format(len(channels)))

        if len(channels) < 1:
            raise ValueError('Too few channels ({0}). Use at least one channel'.format(len(channels)))

        if any((channel_name not in ['1','2']) for channel_name in channels):
            raise ValueError('Wrong channels naming. Channels are named with 1 or 2.')

        if any(((power_unit not in ['DBM','Watt']) and not power_unit == None) for power_unit in power_units):
            raise ValueError('Wrong power units naming. Allowed power units are the strings DBM and Watt.')
            
        if len (wavelengths) <1:
            raise ValueError('Too few wavelength arguments ({0}). The number of wavelengths arguments must be at least  1'.format(len(wavelengths)))
        
        if len (wavelengths) >2 :
            raise ValueError('Too many wavelength arguments ({0}). A maximum of 2 arguments is permitted.'.format(len(wavelengths)))

        if len (wavelengths) != len (channels) :
            raise ValueError ('Lengths of wavelength and channels lists must be equal.')
        
        # for ch_idx , ch in enumerate (channels):
        #     if ch==1 :
        #         if wavelengths[ch_idx] <850 or wavelengths[ch_idx] >1700 :
        #             raise ValueError ('The renge of wavelength is not correct for Hp 81533A that must be between 850 nm to 1700 nm')
        #     else :
        #          if wavelengths[ch_idx] <800 or wavelengths[ch_idx] >1700 :
        #             raise ValueError ('The renge of wavelength is not correct for Hp 81531A that must be between 800  nm to 1700 nm') 
        

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
        lwm= rm.open_resource('GPIB' + str(GPIB_bus) + '::' + str(GPIB_address) + '::INSTR')
    except Exception as e:
        logger.error('No connection possible. Check GPIB connection \n  {0}'.format(e))
        return sys.exit()
    

    # =============================================================================
    #  Settings for the analyzer
    # ============================================================================= 
    
    # Note: Page numbers refer to the "Operating and Programming Manual HP8153A Lightwave Multimeter".
    # Create dict with the the keys
    channel_information = dict.fromkeys(channels)
    
    # Acquire used modules
    if verbose_mode:
        # page (6-9)
        used_modules = lwm.query('*OPT?').rstrip('\n')

    for channel,wavelength,power_unit in zip(channels,wavelengths,power_units):
        # This command sets the units in use when an absolute reading is made. This can be dBm (DBM|0) or Watts (Watt|1).
        # Page (8-21)
        if not power_unit == None:
            lwm.write('sense{0:s}:power:unit {1:s}'.format(channel,power_unit))
            
        # set new wavelength 
        # nanometers (NM) , micrometers(UM), meters (M)
        # Page (8-21)
        if not wavelength == None:    
            lwm.write('sense{0:s}:pow:wave {1:f}NM'.format(channel,wavelength))
        
        # Acquire power values
        #page (8-8 , 8-9)
        channel_power_level = float(lwm.query('read{0:s}:power?'.format(channel)))
        
        if verbose_mode:
            #check the wavelength 
            # Page (8-22)
            read_wavelength = float(lwm.query('sense{0:s}:power:wavelength?'.format(channel)))
        
            # Get the module name
            module = used_modules.split(',')[int(channel)-1]
            
            # Get power unit
            #page (8-21)
            read_power_unit = lwm.query('SENSe{0:s}:POWer:UNIT?'.format(channel)).rstrip('\n')
            
            if read_power_unit == '+0':
                read_power_unit = 'DBM'
            if read_power_unit == '+1':
                read_power_unit = 'Watt'
        
            #make dictionary for power level ,power, wavelength and the name of the inserted module
            data_dict={'Power':channel_power_level , 'Unit':read_power_unit , 'Wavelength':read_wavelength, 'Module': module }

        else:
            #make dictionary only for power level
            data_dict={'Power':channel_power_level}


        #write the data in the dictionary 
        channel_information[channel]=data_dict
    

    # closing lwm connection
    lwm.close()
   
    # closing resource manager 
    rm.close()

    return channel_information



def get_opt_power_HP8163B(channels=['1'], ip_address='192.168.1.1'):
    
    rm = visa.ResourceManager('@py')
   
    # pm = rsm.open_resource('TCPIP0::' + ip_address + '::INST1::INSTR')
    # pm = rsm.open_resource('TCPIP0::' + ip_address + '::INST1::INSTR')
    pm = rm.open_resource('TCPIP0::' + ip_address + '::5025::SOCKET', 
                          read_termination='\n', write_termination='\n', 
                          timeout=3000)
    
    # used_modules = pm.query('*OPT?').rstrip('\n')

    # =============================================================================
    #  Settings for the analyzer
    # ============================================================================= 
    
    # Note: Page numbers refer to the "Operating and Programming Manual HP8153A Lightwave Multimeter".
    # Create dict with the the keys
    channel_information = dict.fromkeys(channels)
    
     
    for channel in channels:
        
        
        # Acquire power values
        #page (8-8 , 8-9)
        channel_power_level = float(pm.query('fetch:chan{0:s}:power?'.format(channel)))
                       
        read_wavelength = float(pm.query('sense:chan{0:s}:power:wavelength?'.format(channel)))
    
        # Get the module name
        # module = used_modules.split(',')[int(channel)-1]
        
        # Get power unit
        #page (8-21)
        read_power_unit = pm.query('sense:chan{0:s}:POWer:UNIT?'.format(channel)).rstrip('\n')
        
        if read_power_unit == '+0':
            read_power_unit = 'DBM'
        if read_power_unit == '+1':
            read_power_unit = 'Watt'
    
        #make dictionary for power level ,power, wavelength and the name of the inserted module
        data_dict={'Power':channel_power_level , 'Unit':read_power_unit , 
                   'Wavelength':read_wavelength}

        #write the data in the dictionary 
        channel_information[channel]=data_dict
    

    # closing lwm connection
    pm.close()
   
    # closing resource manager 
    rm.close()

    return channel_information


def set_attenuation_MTA_150(cassettes = ['1'], attenuations = [None], offsets = [None], wavelengths = [None], GPIB_address='12', log_mode = False):

    """
    set_attenuation_MTA_150
    
    Function for setting the attenuation of the JDS Uniphase MTA 150 optical attenuator. This method is able to change and read the 
    values of the attenuator. For the write mode, the desired cassettes and the corresponding attenuation, offset and wavelength must
    be specified. For the read mode, only the desired cassettes must be specified. The default setting is to read the values from cassette 1. 
    It is also possible to change only individual parameters. In this case, the parameters that are not to be changed receive a None. 
    
    
    Information:
    To change back the MTA to local mode (Controlling with keys), the LCL key must be pressed.
    
    Information:
    The built-in beam block in each MTA300 cassette is automatically activated when the cassette is 
    powered up.The beam block must be deactivated after power-up so that light can passthrough the attenuator.
    This method will not do this, so it must be done manualy at the device.

    Parameters
    ----------
        cassettes: list of strings, optional (default = ['1'])
            The attenuator has several cassettes, the attenuation of each can be individually adjusted. To select the wanted cassette, the 
            numerical index must be put into the list as string. If several cassettes are used, only the numerical indices must be 
            transferred as a list. For example: ['1','2','3'] (Three cassettes).
            Maximum number of cassettes is 8.
            


            WARNING:
            If a cassette is selected which is not physically available, the last selected cassette is used. There is a risk that 
            values ​​will be overwritten.
            

        attenuations : list of floats or Nones, optional (default = [None])
            Sets the total attenuation to the parameter value by changing the actual attenuation. 
            Value must be between 0dB and 60dB
            If not used, the value of the MTA is unchanged.
            If the value of one cassette is to be changed and the others not, a None can simply be inserted in the vector for the value
            that is not to be changed. For example: [20,None,30] (Value of the second entry will not changed)
        

        offsets: list of floats, optional (default = [None])
            Sets the display offset of the MTA system. The value of the offset has no affecton to the actual attenuation,
            but it does affect the total attenuation.
            The display offset function can be used to include both the insertion loss of theMTA300 cassette and
            connection losses in the attenuation value displayed. 
            Att_total = Att_actual + Offset
            Value must be between -60dB and 60dB.
            If not used, the value of the MTA is unchanged.
            If the value of one cassette is to be changed and the others not, a None can simply be inserted in the vector for the value
            that is not to be changed. For example: [20,None,30] (Value of the second entry will not changed)
            
            
        wavelengths = list of floats, optional (default = [None])
            Sets the calibration wavelength of the MTA system. Because the calibration wavelength is used to account for the wavelength 
            dependence of the attenuation, the calibration wavelength should be set as close as possible to the source wavelength.
            Value must be between 1200nm and 1700nm.
            If not used, the value of the MTA is unchanged.
            If the value of one cassette is to be changed and the others not, a None can simply be inserted in the vector for the value
            that is not to be changed. For example: [1300,None,1200] (Value of the second entry will not changed).
     
        GPIB_address : string, optional (default = '13')
            The address GPIB address of the OSA.
        
        log_mode: boolean, optional (default = False)
            Enables a log file for this method.

    Returns
    -------
        cassette_information : dict
            Consist of dicts which contains the attenuation, offset, wavelength and total attenuation of the sected cassette.
            To access the dict use:
            >Name of object<[>Name of cassette<][>Name of data<]
                -> Name of Trace: (string)
                    -> 1 : Cassette 1
                    -> 2 : Cassette 2
                    ...
                    -> 8 : Cassette 8
                -> Name of data: (string)
                    -> attenuation      : (float) Contains the selected attenuations (Att_actual)
                    -> offset           : (float) Contains the selected offset
                    -> wavelength       : (float) Contains the selected wavelength
                    -> Total attenuation: (float) Contains the total attenuation (Att_total = Att_actual + Offset) Corresponds with the 
                                          showed attenuation of the device.
            
    Examples
    -------

        >>> import comm as comm
        
        Get the data from the cassette 1 and 2. The properties of the device will not be changed.
        >>> a = comm.instrument_control.set_attenuation_MTA_150(cassettes=['1','2'])

        Change values of attenuation, offset and wavelength for cassette 2
        >>> b = comm.instrument_control.set_attenuation_MTA_150(cassettes=['2'],attenuations=[5.0], offsets=[5.0], wavelengths=[1300.5])

        Change attenuation of cassette 1 and wavelength of cassette 2
        >>> c = comm.instrument_control.set_attenuation_MTA_150(cassettes=['1','2'],attenuations=[5.0,None],wavelengths=[None,1550.0])

        Cassettes order can be changed
        >>> d = comm.instrument_control.set_attenuation_MTA_150(cassettes=['2','1'],attenuations=[None,0.0],wavelengths=[1500.0,None])

        Access actual attenuation value of cassette 1
        >>> e = comm.instrument_control.set_attenuation_MTA_150(cassettes=['1','2'])
        >>> attenuation_value = e['1']['attenuation']


    Errors
    -------
        Type Error: 
            Will be raised when a wrong data type is used for the input parameter
            -> Possible errors
                -> cassets, attenuations, offsets or wavelengths are not of type list
                -> Items of cassets are not of type string
                -> Items of attenuations, offsets or wavelengths are not of type float
                -> ip_address is not of type string
                -> number_of_bytes is not integer

        Value Error:
            Will be raised when the input parameter is in an wrong range
            -> Possible errors
                -> Too many cassets are used. Maximum is 8
                -> Too few cassets are used. Minimum is 1 
                -> Cassett numbers must be between 0 and 7
                -> Attenuation value is wrong. Must be between 0 and 60
                -> Offset value is wrong. Must be between -60 and 60
                -> Wavelength value is wrong. Must be bewteen 1200 and 1600

        Exception:
            Will be raised by diverse errors
            -> Possible errors
                -> No connection to the attenuator
                -> Required cassettes are not physically present

    """
    # TODO: Finding a way to check which cassettes are connected
        
	# =============================================================================
    #  Create logger which writes to file
    # ============================================================================= 
    # Create logger
    logger = logging.getLogger(__name__)
    
    # Set the log level
    logger.setLevel(logging.INFO)

	# Create standard output handler
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.ERROR)
    
    # Set format of the logs with formatter
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(name)s :: Line No %(lineno)d:: %(message)s')
    
    # Adding formatter to handler
    stdout_handler.setFormatter(formatter)
    
    # Adding handler to logger
    logger.addHandler(stdout_handler)
	

    if log_mode == True:
		# Create file handler 
        file_handler = logging.FileHandler('{0}.log'.format(__name__))
        file_handler.setLevel(logging.INFO)
		
		# Adding formatter to handler
        file_handler.setFormatter(formatter)
		
		# Adding handler to logger
        logger.addHandler(file_handler)

    # =============================================================================
    #  Check inputs for correctness
    # ============================================================================= 
    
    # Check if a None is in the input parametrs of attennuations, offsets or wavelengths. If the statement is true, the value of the 
    # MTA will not be changed.
    # attenuation_unchanged = all(None in attenuations)
    # offset_unchanged = all(None in offsets)
    # wavelength_unchanged = all(None in wavelengths)

    try:
        if not isinstance(cassettes, list):
            raise TypeError('Type of cassettes must be list')

        if not isinstance(attenuations, list):
            raise TypeError('Type of attenuations must be list')

        if not isinstance(offsets, list):
            raise TypeError('Type of offsetts must be list')

        if not isinstance(wavelengths, list):
            raise TypeError('Type of wavelengths must be list')

        if not isinstance(GPIB_address, str):
            raise TypeError('Type of GPIB_address must be string')

        if not all(isinstance(x, str) for x in cassettes):
            raise TypeError('Type of cassettes items must be strings')

        # If no parameters are passed for attenuations, offsets or wavelengths, the list lengths must be adjusted to the length of the cassette list.
        if all(x == None for x in attenuations):
            attenuations = [None]*len(cassettes)

        if all(x == None for x in offsets):
            offsets = [None]*len(cassettes)

        if all(x == None for x in  wavelengths):
            wavelengths = [None]*len(cassettes)

        if not all((isinstance(x, float) or x == None) for x in attenuations):
            raise TypeError('Type of attenuations items must be floats')

        if not all((isinstance(x, float) or x == None) for x in offsets):
            raise TypeError('Type of offsets items must be floats')    

        if not all((isinstance(x, float) or x == None) for x in wavelengths):
            raise TypeError('Type of wavelengths items must be floats')      

        if len(cassettes) > 8:
            raise ValueError('Too many list items ({0}). The MTA_150 has maximal 8 cassettes'.format(len(cassettes)))

        if len(cassettes) < 1:
            raise ValueError('Too less list items ({0}). Use at least one item'.format(len(cassettes))) 

        if any((cassette_name not in ['1','2','3','4','5','6','7','8']) for cassette_name in cassettes):
            raise ValueError('Wrong cassette naming. Cassettes are named with the numbers 1 to 8.')

        #if not any(attenuation_unchanged,offset_unchanged,wavelength_unchanged):

        if not (len(cassettes) == len(attenuations) == len(offsets) == len(wavelengths)) :
            raise ValueError('List length of cassettes, attenuations, offsetts and wavelengths must be the same')

        for attenuation in attenuations:
            if attenuation != None:
                if (attenuation < 0 or attenuation > 60):
                    raise ValueError('Attenuation must be in range of 0 to 60dB')

        for offset in offsets:
            if offset != None:
                if (offset < -60 or offset > 60):
                    raise ValueError('offset must be in range of -60 to 60dB')
                    
        for wavelength in wavelengths:
            if wavelength != None:
                if (wavelength < 1200 or wavelength > 1700):
                    raise ValueError('Wavelengths must be in range of 1200 nm to 1700 nm')

        # if any(((attenuation < 0 or attenuation > 60) if (attenuation != None)) for attenuation in attenuations):
        #     raise ValueError('Attenuation must be in range of 0 to 60dB')

        # if any(((offset < 0 or offset > 60)if offset != None) for offset in offsets):
        #     raise ValueError('offset must be in range of 0 to 60dB')
            
        # if any(((wavelength < 1200 or wavelength > 1700) if wavelength != None)for wavelength in wavelengths):
        #     raise ValueError('Wavelengths must be in range of 1200nm to 1700nm')

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
        attenuator = rm.open_resource('GPIB0::' + GPIB_address + '::INSTR')
    except Exception as e:
        logger.error('No connection possible. Check GPIB connection \n  {0}'.format(e))
        return sys.exit()

    # =============================================================================
    #  Settings for the analyzer
    # =============================================================================  

    # Check if the selected cassettes are present


    # Create return dictionary
    cassette_information = dict.fromkeys(cassettes)

    for cassette, attenuation, offset, wavelength in zip(cassettes,attenuations,offsets,wavelengths):

        # Choosing cassette
        # Page 51
        attenuator.write(':INSTRUMENT:NSELECT {0:s}'.format(cassette))

        # Set offset
        # Page 49
        if not offset == None:
            attenuator.write(':INPUT:OFFSET {0:f}'.format(offset))

        # Set actual attenuation
        # Page 49
        if not attenuation == None:
            attenuator.write(':INPUT:ATTENUATION {0:f}'.format(attenuation))

        # Set wavelength
        # Page 50
        if not wavelength == None:
            attenuator.write(':INPUT:WAVELENGTH {0:f} nm'.format(wavelength))

        # Read values from device
        actual_total_attenuation = float(attenuator.query(':INPUT:ATTENUATION?'))
        actual_offset = float(attenuator.query(':INPUT:OFFSET?'))
        actual_wavelength = float(attenuator.query(':INPUT:WAVELENGTH?'))
        

        # Calculate total attenuation
        actual_attenuation = actual_total_attenuation - actual_offset

        # Create dictionary with data from the selected cassette
        cassette_data = {'Attenuation': actual_attenuation, 'Offset' : actual_offset, 'Wavelength' : actual_wavelength, 'Total attenuation' : actual_total_attenuation}

        # Putting cassette data into the return dictionary
        cassette_information[cassette] = cassette_data

        # Activate Outputs?

    # Change the MTA back to local mode. Local Mode means that the MTA can be controlled by his key on the front.
    # When the MTA is accessed by PC, it will always change to remote mode and the keys doesn't work anymore.
    # The first parameter is the session id and second parameter is the mode. There are six different modes.
    # With mode 2 the MTA will go to local and deassert REN. For more information to the modes see:
    # https://zone.ni.com/reference/en-XX/help/371361R-01/lvinstio/visa_gpib_control_ren/
    rm.visalib.viGpibControlREN(attenuator.session,3)
    
    # closing OSA connection
    attenuator.close()
   
    # closing resource manager 
    rm.close()  

    return cassette_information   
