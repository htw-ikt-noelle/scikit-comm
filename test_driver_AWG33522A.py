# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 15:49:10 2020

@author: UET-Labor
"""

import numpy as np
import random
from writeWfmToAWG33522A import write_samples_AWG33522A



channels = [1]
ipAdress = '192.168.1.44'
Amp_PP = 3.0
Offset = 0.0
SampleRate = 52.5e6
nSamples = 100000

#sample generation using numpy array method
# random float samples generated, clipped between -1 and 1 and multiplied by 32767 and set type as integer
array_samples= np.round(np.clip(np.random.randn(1,nSamples),-1.0,1.0)*32767).astype(int)
list_samples = array_samples.tolist() # changes array into list -- shape [[x,y,z...]]
samples = list_samples.pop() # pops the list out of the list -- shapes[ x,y,z...]








