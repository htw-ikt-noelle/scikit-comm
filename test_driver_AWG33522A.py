# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 15:49:10 2020

@author: UET-Labor
"""

import numpy as np
import random
from writeWfmToAWG33522A import write_samples_AWG33522A

channels = [1,2]
ipAdress = '192.168.1.44'
Amp_PP = 3.0
Offset = 1.0
SampleRate = 52.5e6
nSamples = 100

# random float generator , clipped between lower and higher value, set to 32767 blocks, rounded and output as integer

def sample_floats(low, high, k=1):
    """ Returns a k-length list of unique random floats
        in the range of low <= x <= high
    """
    result = []
    seen = set()
    for i in range(k):
        x = random.uniform(low, high)
        while x in seen:
            x = random.uniform(low, high)
        seen.add(x)
        result.append(x)
    return result
result = sample_floats(-1,1,k = nSamples)

# multiplying each float samples from -1 to 1 with 32767 as the AWG accepts list from -32767 to +32767

product = [i *32767 for i in result]

# converting float values to integer using round
samples = [round(i) for i in product]

#samples using numpy array method
# samples = [np.round(np.clip(np.random.randn(1,nSamples),-1.0,1.0)*32767).astype(int)]
# samples = array_sample.tolist()







