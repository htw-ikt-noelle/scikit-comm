import numpy as np
import comm as comm

# contruct signal
sig = comm.signal.Signal(n_dims=1)
sig.symbol_rate = 32e9

# generate bits
sig.generate_bits(n_bits=2**10)
# set constellation (modualtion format)
sig.generate_constellation(order=4)
sig.modulation_info = 'QPSK'
# create symbols
sig.mapper()
# upsampling and pulseshaping
sig.pulseshaper(upsampling=16, pulseshape='rrc', roll_off=0.2)

# add noise
# sig.set_snr(snr_dB=[35,20])

sig.plot_constellation(0)
sig.plot_eye(0)

# rx filter (matched filtering)
sig.raised_cosine_filter(roll_off=0.2, root_raised=True)
sig.plot_eye(0)

# TODO!!!!
# sampling
for i, samples in enumerate(sig.samples):
    sig.samples[i] = samples[14::16]
    sig.sample_rate[i] = sig.symbol_rate[i]


# sig.plot_eye(0)
# sig.plot_spectrum(dimension=0)
sig.plot_constellation(dimension=0)



