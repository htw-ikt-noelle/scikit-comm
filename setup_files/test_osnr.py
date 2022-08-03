import numpy as np
import matplotlib.pyplot as plt
import time
import comm

from ID_OSA_API_python3 import ID_OSA_CTRL_fn

# osa = ID_OSA_CTRL_fn(IP_address=b'192.168.1.22')

# osa.sweep()
tic = time.time()
# output = comm.instrument_control.get_spectrum_IDOSA(ip_address='192.168.1.22', 
#                                                     start_wl=193.3e12, stop_wl=193.5e12, step=1e9)

output = comm.instrument_control.get_spectrum_IDOSA(ip_address='192.168.1.22', 
                                                      start_wl=193.1e12, stop_wl=193.7e12, step=0.3125e9) #start_wl=191.251e12, stop_wl=196.124e12, step=0.3125e9)

# output = comm.instrument_control.get_spectrum_IDOSA(ip_address='192.168.1.22', 
#                                                      start_wl=191.3e12, stop_wl=196e12, step=0.3125e9) #start_wl=191.251e12, stop_wl=196.124e12, step=0.3125e9)

print('Runtime: {:2.5f} s'.format((time.time()-tic)))


wl = output['wavelength']*1e9
power = output['power']
rbw = np.unique(np.diff(wl))[0]

# c0 = 299792458
# wl = c0/freq
# df = np.unique(np.diff(output['frequency']))[0]


# plt.figure()
# plt.plot(freq/1e12, power)
plt.figure()
plt.plot(wl, power)

OSNR_01nm,OSNR_val = comm.utils.estimate_osnr_spectrum(power_vector=power, wavelength_vector=wl, 
                                  interpolation_points=np.asarray([1549.6,1549.8,1550.4,1550.6]), 
                                  integration_area=np.asarray([1550.05,1550.18]), 
                                  resolution_bandwidth=rbw, polynom_order=1,plotting=True)

print(OSNR_01nm) 
print(OSNR_val)
 
    






