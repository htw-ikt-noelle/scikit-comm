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

y1 = np.repeat(-10,768)
y2 = np.repeat(-3,512)
y3 = np.repeat(-10,768)

y = np.concatenate((y1,y2,y3),axis = 0)

# calculate OSNR
OSNR = comm.osnr.osnr(power_vector = y,
                    wavelength_vector = x,
                    interpolation_points = np.array([737,767,1280,1310]),
                    integration_area = np.array([768,1279]),
                    polynom_order=15)
print(OSNR)

