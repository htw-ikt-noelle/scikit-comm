# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 10:57:06 2021

@author: UET-Labor
"""

import sys
import os
if not any(os.path.abspath('..') == p for p in sys.path): 
    print('adding comm module to path...')
    sys.path.insert(0, os.path.abspath('..'))
import numpy as np
import matplotlib.pyplot as plt
import comm as comm

# Get trace from OSA
OSA_trace_dict = comm.instrument_control.get_samples_HP_71450B_OSA()
#print(OSA_trace_dict)
power = OSA_trace_dict['A']['Trace_data']
wl = OSA_trace_dict['A']['WL_Vector']
rbw = OSA_trace_dict['A']['Resolution_BW']*1e9  # Muss in Nanometer umgerechnet werden


plt.plot(wl,power)
plt.show()

# calculate OSNR
OSNR = comm.osnr.osnr(power_vector = power,
                    wavelength_vector = wl,
                    interpolation_points = np.array([1531,1532,1533,1534]),
                    integration_area = np.array([1532.1,1532.69]),
                    resolution_bandwidth = rbw,
                    polynom_order=3)

names = ['OSNR:','OSNR_1nm:','Signal power:','Noise power:']
units = ['dB','dB','dBm','dBm']

for name,element,unit in zip(names,OSNR,units):
    print(str(name) + ' '+ str(element) + ' ' + unit)

