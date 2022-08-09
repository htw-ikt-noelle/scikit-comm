# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 01:09:52 2022

@author: UeT-Labor
"""
import comm
import numpy as np
import matplotlib.pyplot as plt

out = comm.instrument_control.get_spectrum_IDOSA(ip_address='192.168.1.22',new_sweep=False, wl_equidist=True)

print('RBW (by ID-OSA): {:3.2f} MHz / {:3.5f} nm'.format(out['Resolution_BW_Hz']/1e6,out['Resolution_BW_m']/1e-9))
print('total opt. power (ID-OSA / spectrum integration): {:3.3f} / {:3.3f} dBm'.format(out['Ptotal_dBm_IDOSA'],out['Ptotal_dBm_int']))

#plt.plot(out['WL_vector']/1e-9,out['Trace_data']); plt.show()

OSNR_01nm,OSNR_val = comm.utils.estimate_osnr_spectrum(power_vector=out['Trace_data'], wavelength_vector=out['WL_vector_m']/1e-9, 
                                  interpolation_points=np.asarray([1548.3, 1549.6, 1550.7, 1552.0]), 
                                  integration_area=np.asarray([1549.7,1550.5]), 
                                  resolution_bandwidth=out['Resolution_BW_m']/1e-9, polynom_order=1,plotting=True)

print('OSNR (0.1nm): {:3.2f} dB'.format(OSNR_01nm))
print('OSNR: {:3.2f} dB'.format(OSNR_val))
 