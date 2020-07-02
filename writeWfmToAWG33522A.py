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
Input pareameters [samples, ipadress, channel, amplitudepp = 1, offset = 1, filter = Normal, samplerate = 62.5e6, MEMname set arb1 or arb2, 
                   delete volatile memory to refeed again, clear volatile meory"""

def write_samples_AWG33522A(ipAdress,channels,nSamples,samples,SampleRate,Amp_PP, Offset, filter_settings = None):
    
    import visa
    
    # create resource 
    rm = visa.ResourceManager('@py')
    
    # open connection to AWG
    awg = rm.open_resource('TCPIP::' + ipAdress + '::INSTR')
    
    # check instrument IDN
    idn = awg.query('*IDN?')
    print(idn)
    
    # representing channels 
    # channels = [1,2,3,4]
    
    # selecting byte order , used to make binary data point transfers in the block mode Swapped(LSB) or Normal(MSB)
    # SWAPped byte order,(LSB) of each data point is assumed first. Most computers use the "swapped" byte order.
   
    awg.write(':FORMat:BORDer %s' % ('SWAPped'))
    
    # clearing the waveform memory of the specified channel, 1 or 2 if any.
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
    
    # find which filter has been applied, default is always step filter
    
    before_filter = awg.query('SOUR1:FUNCtion:ARBitrary:FILTer?')
    
    # find what is the sample rate
    sample_rate = awg.query_ascii_values(':SOURce1:FUNCtion:SHAPe:ARBitrary:SRATe?')
    
    #filters=["OFF\n", "STEP\n", "NORMAL\n "]
    
    # Giving parameters as output, just for knowledge. 
    
    print("Here are your AWG configutations:")
    print("IP Adress :",ipAdress)
    print("Amplitude :",Amp_PP)
    print("Offset :",Offset)
    print("SampleRate:",SampleRate)
       
    #if filter_name = awg.query('SOUR1:FUNCtion:ARBitrary:FILTer?') == "OFF\n":
    filter_settings = str(input("Chosse from one of these filters normal or off : "))

    if SampleRate > 62.5e6 and filter_settings =="off":
        raise ValueError('sample rate should be lower than 62.5 Million samples per second')
    
    if SampleRate <= 62.5e6 and filter_settings == "off":
        command = str(input("Do you want to change the STEP filter to off? Type yes or no. "))
        
        if command == "yes":
        
            filter_gen = awg.write('SOUR1:FUNCtion:ARBitrary:FILTer OFF ')
            awg.write(':SOURce1:FUNCtion:SHAPe:ARBitrary:SRATe %G' % (SampleRate))
        if command == "no":
            print("The filter is is still step")
    
    if filter_settings  == "normal":            
        if SampleRate > 62.5e6 or SampleRate <= 62.5e6 : 
            command = str(input("Do you want to change the filter to Normal? If Yes type yes. "))
            if command == "yes":
            #if filter_name = awg.query('SOUR1:FUNCtion:ARBitrary:FILTer?') == "NORMAL\n":
            #if set_filter == "NORMAL":
                awg.write(':SOURce1:FUNCtion:SHAPe:ARBitrary:SRATe %G' % (SampleRate))
                filter_gen = awg.write('SOUR1:FUNCtion:ARBitrary:FILTer NORMAL ') 
           
    ##if filter_name= awg.query('SOUR1:FUNCtion:ARBitrary:FILTer?') == "STEP\n":
        else:
        
            filter_gen = awg.write('SOUR1:FUNCtion:ARBitrary:FILTer STEP ')
       
    after_filter =awg.query('SOUR1:FUNCtion:ARBitrary:FILTer?')
    
    return   before_filter, after_filter
                   
    awg.close()
    rm.close()
        
        
    
    


        
    
    
    
        
            
    
    
            
    
    
    
    
    
