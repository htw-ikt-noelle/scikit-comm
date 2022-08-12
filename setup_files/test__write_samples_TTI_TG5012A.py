# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 17:12:41 2022

@author: UeT-Labor
"""
import numpy as np
import time
import comm as comm

ipadress = '192.168.1.104'
N = int(2**14)
SAMPLES = np.linspace(-1,1,N,endpoint=True)**3 - 0.5*0


comm.instrument_control.write_samples_TTI_TG5012A(ip_address=ipadress, samples=SAMPLES, channel=1, memory='ARB1',
                    waveform='arb',amp_pp=2, repetition_freq=1e3, mute_output=False, interpolate='Off')

time.sleep(0.5)

SAMPLES = -SAMPLES 
SAMPLES = np.linspace(-1,1,N,endpoint=True)**4 - 0.5

comm.instrument_control.write_samples_TTI_TG5012A(ip_address=ipadress, samples=SAMPLES, channel=2, memory='ARB2',
                    waveform='arb',amp_pp=2, repetition_freq=1e3, mute_output=False, interpolate='Off')