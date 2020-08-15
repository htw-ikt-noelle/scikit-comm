# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 12:55:15 2020

@author: Prabesh 

This file works as python driver for Agilent 33522A a dual channel arbitrary Waveform generator.

The test_driver_AWG33522.py will help you to execute the driver, data can be changed in the same file.

"""
import numpy as np
import math
import matplotlib.pyplot as plt
import visa
import time

# import comm
 
def write_samples_AWG33522A(user_attributes, samples,filter_mode = " "):  

# # =============================================================================
# #    check if parameters are available
# # =============================================================================
    
    # if not hasattr(user_attributes,'ipAdress'):
    #     raise AttributeError ('ipAdress has to be specified...')
    # if not isinstance(user_attributes.ipAdress, str): 
    #     raise TypeError('IP Adress must be a string in 000.000.000.000 format')   
    # if not hasattr(user_attributes, 'channel'):
    #     raise ValueError ('channel has to be specified...')
    # if not isinstance(user_attributes.channel, list)
    
    # if not hasattr(user_attributes, 'Amp_PP'):
    #     raise ValueError ('Amplitude has to be specified...')
    # if not isinstance(user_attributes.Amp_PP, float): 
    #     raise TypeError('Amplitude must be a float')
    
    # if not hasattr(user_attributes, 'Offset'):
    #     raise ValueError ('Offset has to be specified...')
    # if not isinstance(user_attributes.Offset, float): 
    #     raise TypeError('Offset must be a float')
    
    # if not hasattr(user_attributes, 'SampleRate'):
    #     raise ValueError ('Samplerate has to be specified...')
    # if not isinstance(user_attributes.SampleRate, float): 
    #     raise TypeError('Samplerate must be a float')
        
    # if not hasattr(user_attributes, 'nSamples'):
    #     raise ValueError ('Number of samples has to be specified...')
    # if not isinstance(user_attributes.nSamples, int): 
    #     raise TypeError('nSamples must be an integer')
        
    # if user_attributes.SampleRate < 1 or user_attributes.SampleRate > 250e6:
    #     raise ValueError('Device cannot perform sample rates lower than 1 and higher than 250Msa/s ')
        
    # if user_attributes.Amp_PP < 1e-3 or user_attributes.Amp_PP > 10 :
    #     raise ValueError(' Vpp can only be from minimum 1mV to maximum 10V..')
        
    # if user_attributes.Offset < -5 or user_attributes.Offset > 5 :
    #     raise ValueError(' min -5V and max +5V into 50 ohm impedance')
        
    # if user_attributes.SampleRate > 62.5e6 and filter_mode =="OFF":
    #     raise ValueError('For OFF filter mode, sample rate should be lower than 62.5 Million samples per second')
        
    
        
# =============================================================================
#  importing visa for communication with the device
# ============================================================================= 
    # create resource 
    rm = visa.ResourceManager('@py')
    # open connection to AWG
    awg = rm.open_resource('TCPIP::' + getattr(user_attributes,'ipAdress') + '::INSTR')   
# =============================================================================
#    setting up channel 1 
# =============================================================================
    # selecting byte order , used to make binary data point transfers in the block mode Swapped(LSB) or Normal(MSB)
    # SWAPped byte order,(LSB) of each data point is assumed first. Most computers use the "swapped" byte order.
    awg.write(':FORMat:BORDer %s' % ('SWAPped'))
    
    
    for ch_idx, ch in enumerate(user_attributes.channel):  
        
        for idx, el in enumerate(samples[ch_idx]):
            if el >= 1:
                samples[ch_idx][idx] = 1
            if el < -1:
                samples[ch_idx][idx] = -1
            samples[ch_idx][idx] *= 32767
            samples[ch_idx][idx]  = int(samples[ch_idx][idx])
         # TO DO !! may be take np array and do the scaling and change to list 
        
        time.sleep(0.1)
         
        awg.write(':SOURce%d:VOLTage:LEVel:IMMediate:COUPle:STATe OFF' % (ch)) 
        
        # output to off is necessary, otherwise the Amplitude is set automatically to 10V and executed which is not good for the instrument 
        # output set to off/ output will be automatic activated after feeding data
        awg.write(':OUTPut%d OFF' % (ch)) 
        
        # clearing the waveform memory of the specified channel
        awg.write(':SOURce%d:DATA:VOLatile:CLEar' % (ch))
        
        # writing values representing DAC codes into waveform volatile memory, as binary block data/ list of integer samples from -32767 to +32767.
        # loading data into the AWG as arb%d, where d = 1or 2 taken from the list of channel
        awg.write_binary_values(':SOURce%d:DATA:ARBitrary:DAC arb%d,' % (ch, ch),samples[ch_idx], datatype = 'h', is_big_endian = False )
        
        
        # setting function as ARB, data will be saved under arb1 file
        awg.write(':SOUR%d:FUNC:SHAP:ARBitrary "arb%d"' % (ch, ch))
       
        # applying sample rate, amplitude and Offset
        #awg.write(':SOURce1:APPLy:ARBitrary %G,%G,%G' % (SampleRate, user_attributes.Amp_PP, user_attributes.Offset))
        awg.write(':SOURce%d:APPLy:ARBitrary %s,%s,%s' % (ch, user_attributes.SampleRate[ch_idx], user_attributes.Amp_PP[ch_idx], user_attributes.Offset[ch_idx]))
        
        if user_attributes.SampleRate[ch_idx] <= 62.5e6 and filter_mode[ch_idx] == "OFF": 
            filter_gen = awg.write(':SOURce%d:FUNCtion:SHAPe:ARBitrary:FILTer %s' % (ch, filter_mode[ch_idx]))
            
            #filter_gen = awg.write('SOUR1:FUNCtion:ARBitrary:FILTer OFF')
            awg.write(':SOURce%d:FUNCtion:SHAPe:ARBitrary:SRATe %G' % (ch, user_attributes.SampleRate[ch_idx]))             
        
        if filter_mode[ch_idx]  == "NORM":            
            if user_attributes.SampleRate[ch_idx] > 62.5e6 or user_attributes.SampleRate[ch_idx] <= 62.5e6 : 
                awg.write(':SOURce%d:FUNCtion:SHAPe:ARBitrary:SRATe %G' % (ch, user_attributes.SampleRate[ch_idx]))
                filter_gen = awg.write('SOUR%d:FUNCtion:ARBitrary:FILTer NORM ' % (ch))      
                
        elif filter_mode[ch_idx] == "STEP":
            filter_gen = awg.write('SOUR%d:FUNCtion:ARBitrary:FILTer STEP ' % (ch))         
        
        # synchronising channel 2 to 1 :: one channel should be OK!  
        awg.write(':SOUR%d:FUNCtion:ARBitrary:SYNC ' % (ch))   
        # print("Currently applied filter: ",awg.query('SOUR2:FUNCtion:ARBitrary:FILTer?'))
        
        # awg.write(':OUTPut%d ON' % (ch))
                 
    awg.close() # closing AWG
    rm.close()  # closing resource manager
     
 