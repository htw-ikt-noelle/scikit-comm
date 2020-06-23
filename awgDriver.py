# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 00:01:25 2020

@author: Prabesh
"""

def set_samples_AWG33522A(self,channels=(1), address = '192.168.1.44'):
    import visa 
    import numpy as np
    from scipy import signal as sg
    import comm as comm
        
    rm = visa.ResourceManager('@py')
    awg = rm.open_resource('TCPIP::' + address + '::INSTR')
    idn = awg.query('*IDN?')
    print(idn)
    
    sig = comm.signal.Signal(n_dims=2)
    sig.symbol_rate = 32e9
	
	running = float(awg.query('*STB?')) % 2
    
    if running:
		awg.write('TRIGger:SOUR CH1')
        busy = 1
        while busy:
            busy = float(awg.query('*STB?')) % 2
			
	wfm = np.array([])
    
   # if wfm is SIN, SQUARE, RAMP, TRIANg,PULSE ,PRBSPNX , ARB
   
   if wfm == np.sin: 
       print('This is a sinusoid ')
       return
   if wfm == scipy.signal.square:
       print('This is a square wave')
       return
   if wfm ==  'ARB':
       if len(wfm)<2 or len(wfm)>2**17:
           print('ERROR!! length of waveform must be between 2 and 131072 points...')
           return
       awg.MEM = 'ARB1'
       interpolate = ON # OFF
   
    if ip != address:
        print('Error: IP Adress could not be found, Please check and run again')
        
    #hasattr(x, 'foo')
    
    if waveform == 'ARB':
        if len(sig) < 2 or len(sig) > 2**17:
            print('Error at signal, must be between 2 and 131072 points')
            return
        
        sig = sig - min(sig)
        sig = sig ./ max(sig) .* 2**14
        sig = unit16(sig)
    
    
    for idx, channel in enumerate(channels):
        
        awg.write('')
	
	if len(sig) < 2 or len(sig) > 2**17
		print('Error at signal, must be between 2 and 131072 points')
		
	
	
	
        
    
    
    rm.close()
    del rm
    del awg
    

    
    
    
