import numpy as np
import comm as comm

# print('Hallo Welt!')
# test = np.arange(10)

sig = comm.signal.Signal(n_dims=2)

sig.generate_bits(n_bits=[2**10, 2**10])
sig.generate_constellation(order=[4])
sig.mapper()
sig.pulseshaper(upsampling=[4], pulseshape=['rc'])


# sig.set_snr(snr_dB=[25,10])

sig.plot_constellation(0)
sig.plot_eye(0)
sig.plot_spectrum(dimension=0)


# # comm.visualizer.plot_eye(sig.samples[0], sample_rate=10)

# # sig.plot_spectrum(dimension=0)
# # sig.plot_spectrum(dimension=1)


