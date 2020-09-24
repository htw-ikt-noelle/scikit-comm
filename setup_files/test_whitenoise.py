# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 12:14:04 2020

@author: UET-Labor
"""

import numpy as np

import matplotlib.pyplot as plt
import comm as comm
import numpy as np
from scipy.interpolate import interp1d


N = 1000
noise = np.random.randn(N,2).view(np.complex128)
noise = noise.reshape(-1,)

sr = 100
# plt.figure(1)
#plt.plot(noise); plt.show()
# plt.xlabel('White noise')
# plt.show()
# plt.figure(2)
comm.visualizer.plot_spectrum(noise)
# comm.visualizer.plot_spectrum(np.real(noise))
# comm.visualizer.plot_spectrum(np.imag(noise))

FILTER = np.array([[-0.4*sr,-0.3*sr,0*sr,0.1*sr,0.2*sr], [0,-100,10,-150,-50]],dtype='float').transpose()

plt.figure(1)
plt.plot(FILTER[:,0],FILTER[:,1],'r-o'),plt.xlim((-sr/2,sr/2)); plt.show()


samples_out = comm.filters.filter_arbitrary(noise,FILTER,sample_rate=sr)

plt.figure(2)
comm.visualizer.plot_spectrum(samples_out)



 