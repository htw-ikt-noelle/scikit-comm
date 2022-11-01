import time, copy, sys

import numpy as np
import scipy.signal as ssignal
import matplotlib.pyplot as plt
import scipy.interpolate as sinterp

import comm as comm


sig_range = np.asarray([1549.85,1550.35])*1e-9
noise_range = np.asarray([1549.75, 1549.8, 1550.45, 1550.55])*1e-9

# estimate OSNR
id_out = comm.instrument_control.get_spectrum_IDOSA(new_sweep=False, wl_equidist = False)
hp_out = comm.instrument_control.get_spectrum_HP_71450B_OSA (traces = ['A'], GPIB_bus=0, 
                                                              GPIB_address=13,
                                                              log_mode = False, 
                                                              single_sweep = False)

osnr_id = comm.utils.estimate_snr_spectrum(id_out['WL_vector_m'], id_out['Trace_data'], 
                                  sig_range=sig_range, 
                                  noise_range=noise_range,
                                  order=1, noise_bw=0.1e-9, scaling='log', plotting=True)

osnr_hp = comm.utils.estimate_snr_spectrum(hp_out['A']['WL_Vector'], hp_out['A']['Trace_data'], 
                                  sig_range=sig_range, 
                                  noise_range=noise_range,
                                  order=1, noise_bw=0.1e-9, scaling='log', plotting=True)

print('ID OSNR: {:.1f} dB, HP OSNR: {:.1f} dB'.format(osnr_id, osnr_hp))