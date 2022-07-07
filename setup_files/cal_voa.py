import time

import numpy as np
import matplotlib.pyplot as plt


import comm

voltages = np.arange(0.0, 5.25, 0.05)
attenuation = []

for voltage in voltages:

    comm.instrument_control.write_samples_tti_tg5012a(ip_address='192.168.1.105', waveform='DC',
                                                      amp_pp=voltage, channel=1)
    
    time.sleep(2)
    
    p = comm.instrument_control.get_opt_power_Anritsu_ML910B(GPIB_bus=0, GPIB_address=11)
    
    attenuation.append(-p['ch2']['value'])
    
    print(voltage)
    


plt.plot(voltages, np.asarray(attenuation))