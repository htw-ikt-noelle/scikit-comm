import numpy as np
import comm as comm


sig = comm.signal.Signal(n_dims=2)
sig.symbol_rate = 32e9




sig.generate_bits(n_bits=2**10)
sig.generate_constellation(order=4)
sig.mapper()
sig.pulseshaper(upsampling=16, pulseshape='rc', roll_off=[0.2, 1])


# sig.set_snr(snr_dB=[35,20])

# sig.plot_constellation(0)
sig.plot_eye(0)
sig.plot_eye(1)
# sig.plot_spectrum(dimension=0)

