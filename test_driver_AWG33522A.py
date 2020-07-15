# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 15:49:10 2020

@author: UET-Labor
"""

import numpy as np
import random
from writeWfmToAWG33522A import write_samples_AWG33522A

class user_attributes:  
    
    def __init__(self):
        self.ipAdress = ipAdress
        self.Amp_PP = Amp_PP
        self.Offset = Offset
        self.SampleRate = SampleRate
        self.nSamples = nSamples
        self.Offset = Offset
        self.channels = channels

# addresing the attributes
user_attributes.ipAdress = '192.168.1.44'
user_attributes.nSamples = 50000  # number of samples
user_attributes.SampleRate = 72e6 # not more than 250e6
user_attributes.Offset = 3.0      # range -5V +5V
user_attributes.Amp_PP = 2.0      #range 1e-3 to +10V
user_attributes.channels = 1 # channels can be 1 or 2 as specified by the AWG


#generating samples for 1st channel
# random float samples generated, clipped between -1 and 1 and multiplied by 32767 and set type as integer
array_samples1= np.round(np.clip(np.random.randn(1,user_attributes.nSamples),-1.0,1.0)*32767).astype(int)
samples1 = array_samples1.tolist().pop() # changes array into list -- shape [[x,y,z...]]
# samples = list_samples.pop() # pops the list out of the list -- shapes[ x,y,z...]


# generating samples for channel 2
array_samples2= np.round(np.clip(np.random.randn(1,user_attributes.nSamples),-1.0,1.0)*32767).astype(int)
samples2 = array_samples2.tolist().pop() # changes array into list -- shape [[x,y,z...]]

# channel 1 output, input parameters can be changed
write_samples_AWG33522A(user_attributes, samples1, samples2, SampleRate= 240e6 ,channels = 1, Amp_PP = 4.0, Offset= 1.0, filter_settings = "normal")
#channel 2 output
write_samples_AWG33522A(user_attributes, samples1, samples2, SampleRate= 240e6 ,channels = 2, Amp_PP = 3.0, Offset= 0.0, filter_settings = "normal")

