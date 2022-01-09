import sys
import os
if not any(os.path.abspath('.') == p for p in sys.path): 
    print('adding comm module to path...')
    sys.path.insert(0, os.path.abspath('.'))
import numpy as np
import matplotlib.pyplot as plt
import comm as comm

# Create test vector
x = np.linspace(0,2047,2048)

y1 = np.repeat(-12,768)
y2 = np.repeat(0,512)
y3 = np.repeat(-12,768)

y = np.concatenate((y1,y2,y3),axis = 0)

# Integration area
left_inetgra = 768
right_integra = 1279

# Interpolation area
left_interpol_2 = 767 
left_interpol_1 = left_interpol_2 - 50

right_interpol_1 = 1280
right_interpol_2 = right_interpol_1 + 50

# calculate OSNR
OSNR = comm.osnr.osnr(power_vector = y,
                    wavelength_vector = x,
                    interpolation_points = np.array([437,767,1280,1610]),
                    integration_area = np.array([left_inetgra,right_integra]),
                    resolution_bandwidth= 1,
                    polynom_order=3)
print(OSNR)

names = ['OSNR:','OSNR_1nm:','Signal power:','Noise power:']
units = ['dB','dB','dBm','dBm']

for name,element,unit in zip(names,OSNR,units):
    print(str(name) + ' '+ str(element) + ' ' + unit)