import numpy as np
import matplotlib.pyplot as plt

import comm

# 'SINE', 'SQUARE', 'RAMP', 'TRIANG', 'PULSE', 'NOISE', 'PRBSPNX', 'ARB'

samples = np.sin(np.linspace(0.0, 2*np.pi, 1000, endpoint=False))

plt.plot(samples)

comm.instrument_control.write_samples_tti_tg5012a(samples, ip_address='192.168.1.105', waveform='DC',
                                                  amp_pp=2, channel=1, repetition_freq=50.0, 
                                                  memory='ARB1', interpolate='OFF', bit_rate=500.0)


