# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 12:55:15 2020

@author: UET-Labor

"""
import numpy as np
import math
import matplotlib.pyplot as plt
# import comm

# from comm import signal
# from comm import instrument_control

# !! to DO !!

"""
to do is set different parameters such as samplerates, amplitudes for both channels and setup differently for both channels, channels can be selected by giving channesl name"""

#def write_samples_AWG33522A(user_inputs, ipAdress,channels,nSamples,samples,SampleRate,Amp_PP, Offset, filter_settings = None):
def write_samples_AWG33522A(user_attributes,samples, SampleRate = 62e6,channels = 1, Amp_PP = 3.0, Offset= 1.0, filter_settings = "normal"): 
    
    
    # parameters
    
    if hasattr(user_attributes,'ipAdress') == False:
        raise ValueError ('ipAdress has to be specified...')
        
    if hasattr(user_attributes, 'channels') == False:
        raise ValueError ('channels has to be specified...')
    
    if SampleRate < 1 or SampleRate > 250e6:
        raise ValueError('Device cannot perform sample rates lower than 1 and higher than 250Msa/s ')
        
    if Amp_PP < 1e-3 or Amp_PP > 10 :
        raise ValueError(' Vpp can only be from minimum 1mV to maximum 10V..')
        
    if Offset < -5 or Offset > 5 :
        raise ValueError(' min -5V and max +5V into 50 ohm impedance')
        
             
    import visa
    
    # create resource 
    rm = visa.ResourceManager('@py')
    
    # open connection to AWG
    awg = rm.open_resource('TCPIP::' + getattr(user_attributes,'ipAdress') + '::INSTR')
    
        
    # selecting byte order , used to make binary data point transfers in the block mode Swapped(LSB) or Normal(MSB)
    # SWAPped byte order,(LSB) of each data point is assumed first. Most computers use the "swapped" byte order.
   
    awg.write(':FORMat:BORDer %s' % ('SWAPped'))
    
    # clearing the waveform memory of the specified channel, 1 or 2 if any.
    if channels == 1:
        awg.write(':SOURce1:DATA:VOLatile:CLEar')
    
        # writing values representing DAC codes into waveform volatile memory, as binary block data/ list of integer samples from -32767 to +32767.
        # loading data into the AWG as arb1, can be named anything, just need to check source number 1 or 2
        awg.write_binary_values(':SOURce1:DATA:ARBitrary:DAC %s,' % ('arb1' ),samples, datatype = 'h', is_big_endian = False )
        
        # setting function as ARB
        awg.write(':SOURce1:FUNCtion:SHAPe:ARBitrary "%s"' % ('arb1'))
        
        # applying sample rate, amplitude and offset
        awg.write(':SOURce1:APPLy:ARBitrary %G,%G,%G' % (SampleRate, Amp_PP, Offset))
        
        # synchronising both channels to the first point one channel should be OK!
        
        awg.write(':SOUR1:FUNCtion:ARBitrary:SYNC ')
    
    #set offset if needed 
    # awg.write(':SOURce1:FREQuency:COUPle:OFFSet %G' % (30000.0))
    
    if channels == 2:
        awg.write(':SOURce2:DATA:VOLatile:CLEar')
        
        awg.write_binary_values(':SOURce2:DATA:ARBitrary:DAC %s,' % ('arb1' ),samples, datatype = 'h', is_big_endian = False )
        
        # setting function as ARB
        awg.write(':SOURce2:FUNCtion:SHAPe:ARBitrary "%s"' % ('arb1'))
        
        # applying sample rate, amplitude and offset
        awg.write(':SOURce2:APPLy:ARBitrary %G,%G,%G' % (SampleRate, Amp_PP, Offset))
        
        # synchronising both channels to the first point one channel should be OK!
        
        awg.write(':SOUR2:FUNCtion:ARBitrary:SYNC ')
        
    # find which filter has been applied, default is always step filter
    # before_filter = awg.query('SOUR1:FUNCtion:ARBitrary:FILTer?')
    
    # find what is the sample rate
    # sample_rate = awg.query_ascii_values(':SOURce1:FUNCtion:SHAPe:ARBitrary:SRATe?')
    
            
    if SampleRate > 62.5e6 and filter_settings =="off":
        raise ValueError('For OFF filter Settings, sample rate should be lower than 62.5 Million samples per second')
    
    if SampleRate <= 62.5e6 and filter_settings == "off":
        #command = str(input("Do you want to change the STEP filter to OFF? Type Yes or No. ").lower())
        #if filter == 1:
            
        filter_gen = awg.write('SOUR1:FUNCtion:ARBitrary:FILTer OFF ')
        awg.write(':SOURce1:FUNCtion:SHAPe:ARBitrary:SRATe %G' % (SampleRate))
    # if command == "no":
    #     print("The filter is is still step")
    
    if filter_settings  == "normal":            
        if SampleRate > 62.5e6 or SampleRate <= 62.5e6 : 
            # command = str(input("Do you want to change the filter to Normal? Type Yes or No. ").lower())
            # if command == "yes":
            # #if filter_name = awg.query('SOUR1:FUNCtion:ARBitrary:FILTer?') == "NORMAL\n":
            # #if set_filter == "NORMAL":
            awg.write(':SOURce1:FUNCtion:SHAPe:ARBitrary:SRATe %G' % (SampleRate))
            filter_gen = awg.write('SOUR1:FUNCtion:ARBitrary:FILTer NORMAL ') 
           
    ##if filter_name= awg.query('SOUR1:FUNCtion:ARBitrary:FILTer?') == "STEP\n":
    # else:
    
    #     raise ValueError('other filters not allowed yet.., default output will be step')
       
    print("Currently applied filter: ",awg.query('SOUR1:FUNCtion:ARBitrary:FILTer?'))
  
    
    
                   
    awg.close()
    rm.close()
        
        
    
    


        
    
    
    
        
            
    
    
            
    
    
    
    
    
