import numpy as np
import comm as comm
import matplotlib.pyplot as plt
from edit_writeWfmToAWG33522A import write_samples_AWG33522A
import scipy.signal as signal



# contruct signal
sig = comm.signal.Signal(n_dims=2)
sig.symbol_rate = 5e5

# generate bits
sig.generate_bits(n_bits=2**10)
# set constellation (modualtion format)
sig.generate_constellation(order=4)
sig.modulation_info = 'QPSK'
# create symbols
sig.mapper()
# upsampling and pulseshaping
sig.pulseshaper(upsampling=16, pulseshape='rrc', roll_off=[0.2, 1])

# # add noise
# # sig.set_snr(snr_dB=[35,20])

# sig.plot_constellation(0)
# sig.plot_eye(0)

# # rx filter (matched filtering)
# sig.raised_cosine_filter(roll_off=0.2, root_raised=True)
# sig.plot_eye(0)

# # TODO!!!!
# # sampling
# for i, samples in enumerate(sig.samples):
#     sig.samples[i] = samples[15::16]
#     sig.sample_rate[i] = sig.symbol_rate[i]

class user_attributes:
    def __init__(self):  # initialize 
        self.ipAdress = ipAdress
        self.Amp_PP = Amp_PP
        self.Offset = Offset
        self.SampleRate = SampleRate
        self.nSamples = nSamples
        self.Offset = Offset
        self.channel= channel
        
# check= user_attributes()
# setting up default values
user_attributes.ipAdress = '192.168.1.44'
user_attributes.nSamples = int(1e6) # number of samples
user_attributes.SampleRate = 250e6 # not more than 250e6
user_attributes.Offset = 0.0       # range -5V +5V
user_attributes.Amp_PP = 4.0       # range 1e-3 to +10V
user_attributes.channel = [1,2]    # channels can be 1 or 2 as specified by the AWG

# random float samples generated, clipped between -1 and 1 and multiplied by 32767 and set type as integer
array_samples= np.round(np.clip(np.random.randn(1,user_attributes.nSamples),-1.0,1.0)*32767).astype(int)
samples = array_samples.tolist().pop() # changes array into list -- shape [[x,y,z...]] and pops out one list

# both channel output, input parameters are used from above, filters can be set to NORM, OFF and STEP according to the manual, only use upper case
write_samples_AWG33522A(user_attributes, samples,filter_mode ="NORM")

# downsampling to up times symbolrate
sr_dsp = sig.symbol_rate * 10

import math
IF = 35e6
T = (np.arange(0,np.size(samples),sig.symbol_rate[0])/sig.symbol_rate[0])
# samples_1 = sig.samples*math.e**(1j*IF*T)


# general term 
# r(t) = a(t)* sin(IF + phase + phase noise of laser + static start phase of the laser)

#Sin term r(t) * sin(Wif*t)

samples_Q = sig.samples*math.sin(IF*T)
    

# Cosine Term

samples_I = sig.samples *math.cos(IF *T)

# next LP filter 


# combine sin and cos term_ check result



# sig.plot_eye(0)
# sig.plot_spectrum(dimension=0)



