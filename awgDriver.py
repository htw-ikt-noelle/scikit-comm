# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 00:01:25 2020

@author: Prabesh
"""

def set_samples_AWG33522A(channels=(1), address = '192.168.1.45'):
    import visa 
    import numpy as np
    import comm as comm
    
    rm = visa.ResourceManager('@py')
    awg = rm.open_resource('TCPIP::' + address + '::INSTR')
    idn = awg.ask('*IDN?')
    print(idn)
    
    sig = comm.signal.Signal(n_dims=2)
    sig.symbol_rate = 32e9
	
	running = float(awg.query('STATus:CONDition?')) % 2
    
    if running:
		awg.write('TRIGger:MODE SINGle')
        busy = 1
        while busy:
            busy = float(awg.query('STATus:CONDition?')) % 2
			
	# sig = np.array([])
	
	if len(sig) < 2 or len(sig) > 2**17
		print('Error at signal, must be between 2 and 131072 points')
		
	
	
	
        
    
    
    rm.close()
    del rm
    del awg
    

    
    
    
