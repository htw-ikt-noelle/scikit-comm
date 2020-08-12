# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 15:49:10 2020

@author: Prabesh

This is the test file for the AWG33522A driver. data from the user attributes 
can be changed according to the need.

Parameters are set and the generated data are then written in the AWG with the help of awg_driver.py

ipAdress : ip adress of the AWG equipment -- string

Amp_PP : Peak to Peak voltage  ranges from 1e-3 V  to +10 V : source--manual -- float

nSamples : number of samples --integer

SampleRate: sample rate can be used upto 250 millions samples/s --integer

Offset : Offset frequency --float

channel : denotes the channel 1 and 2 in a --list

samples : randomly generated 32 bit  floating point values

filter_mode : 3 Filter modes are available for this AWG 
            
    NORM : provides the flattest frequency response. This effectively smoothes the signal, 
    but sharp transitions will have pre-shoot and overshoot
    
    OFF : steps from point to point at the sample rate. Moves between data points are 
    accomplished as quickly as possible with no smoothing. If the <mode> is set to OFF, 
    the maximum sample rate for the arbitrary waveform is limited to 62.5 MSa/s. 
    
    STEP : the data points in a way that effectively smoothes the signal while minimizing
    the pre-shoot and overshoot. However, this setting has a narrower bandwidth than
    the NORMal setting. 
"""

import numpy as np
import random
from edit_writeWfmToAWG33522A import write_samples_AWG33522A

class user_attributes:
    def __init__(self):  # initialize 
        self.ipAdress = ipAdress
        self.Amp_PP = Amp_PP
        self.Offset = Offset
        self.SampleRate = SampleRate
        self.nSamples = nSamples
        self.Offset = Offset
        self.channel= channel
        
check= user_attributes()
# setting up default values
user_attributes.ipAdress = '192.168.1.44'
user_attributes.nSamples = int(1e6)  # number of samples
user_attributes.SampleRate = 250e6 # not more than 250e6
user_attributes.Offset = 0.0      # range -5V +5V
user_attributes.Amp_PP = 4.0      #range 1e-3 to +10V
user_attributes.channel = [1,2] # channels can be 1 or 2 as specified by the AWG


# random float samples generated, clipped between -1 and 1 and multiplied by 32767 and set type as integer
array_samples= np.round(np.clip(np.random.randn(1,user_attributes.nSamples),-1.0,1.0)*32767).astype(int)
samples = array_samples.tolist().pop() # changes array into list -- shape [[x,y,z...]] and pops out one list

# both channel output, input parameters are used from above, filters can be set to NORM, OFF and STEP according to the manual, only use upper case
write_samples_AWG33522A(user_attributes, samples,filter_mode ="NORM")


