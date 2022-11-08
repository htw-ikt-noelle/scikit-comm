import numpy as np
import matplotlib.pyplot as plt
import comm as comm
import copy

sig_tx = comm.signal.Signal(n_dims=2)
sig_tx.symbol_rate = 12.8e9
amount_of_symbols = 2**12
mod_order=4

# generate bits
n_bits = int((np.log2(mod_order)*amount_of_symbols)-((np.log2(mod_order)*amount_of_symbols)%np.log2(mod_order)))
sig_tx.generate_bits(n_bits=n_bits, seed=1)

# set constellation (modulation format)
sig_tx.generate_constellation(format='QAM', order=mod_order)

# create symbols
sig_tx.mapper()

 # upsampling and pulseshaping
ROLL_OFF = 0.1
sig_tx.pulseshaper(upsampling=2, pulseshape='rrc', roll_off=[ROLL_OFF, ROLL_OFF])

print(sig_tx.samples[0] == sig_tx.samples[1])

sig_range = np.array([-sig_tx.symbol_rate[0]/2-(sig_tx.symbol_rate[0]/2*ROLL_OFF),sig_tx.symbol_rate[0]/2+(sig_tx.symbol_rate[0]/2*ROLL_OFF)])
noise_range = np.array([-sig_tx.symbol_rate[0]/2-(sig_tx.symbol_rate[0]/2*ROLL_OFF)-1e9,-sig_tx.symbol_rate[0]/2-(sig_tx.symbol_rate[0]/2*ROLL_OFF),sig_tx.symbol_rate[0]/2+(sig_tx.symbol_rate[0]/2*ROLL_OFF),sig_tx.symbol_rate[0]/2+1e9+(sig_tx.symbol_rate[0]/2*ROLL_OFF)])

rng = np.random.default_rng()

snr_list = [4,20]
for dim in range(sig_tx.n_dims):
    #print(snr_list[dim])
    #sig_tx.samples[dim] = comm.channel.set_snr(sig_tx.samples[dim], snr_dB=snr_list[dim], sps=sig_tx.sample_rate[dim]/sig_tx.symbol_rate[dim], seed=None)
    sig_tx.samples[dim] = sig_tx.samples[dim] * np.sqrt(10**(snr_list[dim]/10)) + np.sqrt(0.5)*(rng.standard_normal(size=(sig_tx.samples[dim].size,)) + 1j*rng.standard_normal(size=(sig_tx.samples[dim].size,))) 

sig_tx_mrc = copy.deepcopy(sig_tx)
sig_tx_egc = copy.deepcopy(sig_tx)

sig_comb_mrc = comm.rx.combining(sig_tx_mrc, comb_method="MRC", snr=snr_list)
sig_comb_egc = comm.rx.combining(sig_tx_egc, comb_method="EGC")

x_in = np.fft.fftshift(np.fft.fftfreq(len(sig_comb_mrc.samples[0]),1/sig_comb_mrc.sample_rate))
y_in = np.abs(np.fft.fftshift(np.fft.fft(sig_comb_mrc.samples[0])))**2 
snr_mrc = comm.utils.estimate_snr_spectrum(x_in,y_in,sig_range=sig_range, noise_range= noise_range, order=1,noise_bw=sig_comb_mrc.symbol_rate,plotting=False)

x_in = np.fft.fftshift(np.fft.fftfreq(len(sig_comb_egc.samples[0]),1/sig_comb_egc.sample_rate))
y_in = np.abs(np.fft.fftshift(np.fft.fft(sig_comb_egc.samples[0])))**2 
snr_egc = comm.utils.estimate_snr_spectrum(x_in,y_in,sig_range=sig_range, noise_range= noise_range, order=1,noise_bw=sig_comb_egc.symbol_rate,plotting=False)

snr_list_lin = [10**(x/10) for x in snr_list]
snr_mrc_theory = 10*np.log10(sum(snr_list_lin))

snr_egc_theory = 10*np.log10(1/2*sum(abs(np.sqrt(snr_list_lin)))**2)

print("SNR mrc: {}dB, theory: {}dB".format(snr_mrc, snr_mrc_theory))
print("SNR egc: {}dB, theory: {}dB".format(snr_egc, snr_egc_theory))