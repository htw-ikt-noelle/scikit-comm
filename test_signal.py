import numpy as np
import comm as comm

# print('Hallo Welt!')
# test = np.arange(10)

sig = comm.signal.Signal(n_dims=2)

sig.generate_bits()
sig.generate_constellation()
sig.mapper()
sig.samples = sig.symbols
sig.plot_spectrum(dimension=0)
sig.plot_spectrum(dimension=1)


