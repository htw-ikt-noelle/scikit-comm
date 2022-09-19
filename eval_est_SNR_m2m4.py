import comm 
import numpy as np
import matplotlib.pyplot as plt
import copy

snr_list = np.arange(0,30,1)
snr_est_list = np.zeros(snr_list.shape)

mod_format = 'QAM'
mod_order = 4

#### TX signal block
sig = comm.signal.Signal(n_dims=int(1))
sig.symbol_rate= 5e9
sig.sample_rate= 5e9
sig.generate_bits(n_bits=int(2**14),seed=None)
sig.generate_constellation(format=mod_format,order=mod_order)
sig.mapper()
# only integer upsampling factors for now
sig.pulseshaper(upsampling=1,pulseshape='rrc',roll_off=0.2)

for i in range(len(snr_list)):
    sig_tx = copy.deepcopy(sig)
    sig_tx.samples[0] = comm.channel.set_snr(samples=sig_tx.samples[0] , snr_dB=snr_list[i], seed=None)
    snr_est_list[i] = comm.utils.estimate_SNR_m2m4(sig_tx.samples[0], sig_tx.constellation[0])

plt.figure()
plt.title("SNR est m2m4 validation")
plt.plot(snr_list, 10*np.log10(snr_est_list), label="est. SNR")
plt.plot(snr_list, snr_list, label="theory", ls="--")
plt.ylabel("est. SNR [dB]")
plt.xlabel("set SNR [dB]")
plt.grid()
plt.legend()
plt.show()