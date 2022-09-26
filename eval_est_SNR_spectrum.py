import comm 
import numpy as np
import matplotlib.pyplot as plt
import copy

# ##### SNR Performance in General!

# snr_list = np.arange(0,30,1)
# n_blocks = 1
# snr_est_list = np.zeros((n_blocks, len(snr_list)))

# mod_format = 'QAM'
# mod_order = 4
# roll_off = 0.2
# n_symbols = int(10e3)

# #### TX signal block
# sig = comm.signal.Signal(n_dims=int(1))
# sig.symbol_rate= 5e9
# sig.sample_rate= 15e9
# n_bits = comm.utils.symb_to_bits(n_symbols, format=mod_format, order=mod_order)
# sig.generate_bits(n_bits=n_bits,seed=None)
# sig.generate_constellation(format=mod_format,order=mod_order)
# sig.mapper()
# # only integer upsampling factors for now
# sig.pulseshaper(upsampling=sig.sample_rate[0]/sig.symbol_rate[0],pulseshape='rrc',roll_off=roll_off)

# for i in range(len(snr_list)):
#     for j in range(n_blocks):
#         sig_tx = copy.deepcopy(sig)
#         sig_tx.samples[0] = comm.channel.set_snr(samples=sig_tx.samples[0] , snr_dB=snr_list[i], sps=sig.sample_rate[0]/sig.symbol_rate[0], seed=None)

#         sig_range = np.array([-sig.symbol_rate[0]/2-(sig.symbol_rate[0]/2*roll_off),sig.symbol_rate[0]/2+(sig.symbol_rate[0]/2*roll_off)])
#         noise_range = np.array([-sig.symbol_rate[0]/2-(sig.symbol_rate[0]/2*roll_off)-1e9,-sig.symbol_rate[0]/2-(sig.symbol_rate[0]/2*roll_off),sig.symbol_rate[0]/2+(sig.symbol_rate[0]/2*roll_off),sig.symbol_rate[0]/2+1e9+(sig.symbol_rate[0]/2*roll_off)])

#         x_comb_EGC = np.fft.fftshift(np.fft.fftfreq(len(sig_tx.samples[0]),1/sig_tx.sample_rate[0]))
#         y_comb_EGC = np.abs(np.fft.fftshift(np.fft.fft(sig_tx.samples[0])))**2 
#         snr_est_list[j,i] = 10**(comm.utils.estimate_snr_spectrum(x_comb_EGC,y_comb_EGC,sig_range=sig_range, noise_range=noise_range, 
#                                                         order=1,noise_bw=sig.symbol_rate[0],plotting=False)/10)


# plt.figure(1)
# title_string = "Params: {}-{}, {} runs, \n {:.2E} symbols, symbrate: {:.2E}, samplerate: {:.2E}".format(mod_order, mod_format, n_blocks, n_symbols, sig.symbol_rate[0], sig.sample_rate[0])
# plt.title("SNR spectrum est validation"+title_string)
# plt.plot(snr_list, 10*np.log10(np.mean(snr_est_list, axis=0)), label="est. SNR")
# plt.plot(snr_list, snr_list, label="theory", ls="--")
# plt.ylabel("est. SNR [dB]")
# plt.xlabel("set SNR [dB]")
# plt.grid()
# plt.legend()

# plt.figure(2)
# plt.title("Difference SNR spectrum est - set SNR"+title_string)
# plt.plot(snr_list, 10*np.log10(np.mean(snr_est_list, axis=0))-snr_list, label="est. SNR")
# plt.plot(snr_list, np.zeros(len(snr_list)), label="theory", ls="--")
# plt.ylabel("est. SNR (Diff) [dB]")
# plt.xlabel("set SNR [dB]")
# plt.grid()
# plt.legend()

# plt.show()

# ##### SNR Performance vs n_symbols at SNR = 5dB

n_symbols_list = np.linspace(10,4e4,10)
n_symbols_list = [int(x) for x in n_symbols_list]
n_blocks = 100
snr_est_list = np.zeros((n_blocks,len(n_symbols_list)))

mod_format = 'QAM'
mod_order = 4
roll_off = 0.2
snr = 5

for i in range(len(n_symbols_list)):
    #### TX signal block    
    sig = comm.signal.Signal(n_dims=int(1))
    sig.symbol_rate= 1e6
    sig.sample_rate= 2e6
    n_bits = comm.utils.symb_to_bits(n_symbols_list[i], format=mod_format, order=mod_order)
    sig.generate_bits(n_bits=n_bits,seed=None)
    sig.generate_constellation(format=mod_format,order=mod_order)
    sig.mapper()
    # only integer upsampling factors for now
    sig.pulseshaper(upsampling=sig.sample_rate[0]/sig.symbol_rate[0],pulseshape='rrc',roll_off=roll_off)


    for j in range(n_blocks):
        sig_tx = copy.deepcopy(sig)
        sig_tx.samples[0] = comm.channel.set_snr(samples=sig_tx.samples[0] , snr_dB=snr, sps=sig.sample_rate[0]/sig.symbol_rate[0], seed=None)

        sig_range = np.array([-sig.symbol_rate[0]/2-(sig.symbol_rate[0]/2*roll_off),sig.symbol_rate[0]/2+(sig.symbol_rate[0]/2*roll_off)])
        noise_range = np.array([-sig.symbol_rate[0]/2-(sig.symbol_rate[0]/2*roll_off)-1e9,-sig.symbol_rate[0]/2-(sig.symbol_rate[0]/2*roll_off),sig.symbol_rate[0]/2+(sig.symbol_rate[0]/2*roll_off),sig.symbol_rate[0]/2+1e9+(sig.symbol_rate[0]/2*roll_off)])

        x_comb_EGC = np.fft.fftshift(np.fft.fftfreq(len(sig_tx.samples[0]),1/sig_tx.sample_rate[0]))
        y_comb_EGC = np.abs(np.fft.fftshift(np.fft.fft(sig_tx.samples[0])))**2 
        snr_est_list[j,i] = 10**(comm.utils.estimate_snr_spectrum(x_comb_EGC,y_comb_EGC,sig_range=sig_range, noise_range=noise_range, 
                                                        order=1,noise_bw=sig.symbol_rate[0],plotting=False)/10)


plt.figure(1)
title_string = " Params: {}-{}, {} runs, \n snr= {}dB, symbrate: {:.2E}, samplerate: {:.2E}".format(mod_order, mod_format, n_blocks, snr, sig.symbol_rate[0], sig.sample_rate[0])
plt.title("SNR spectrum est validation"+title_string)
plt.plot(n_symbols_list, 10*np.log10(np.mean(snr_est_list, axis=0)), label="est. SNR")
plt.plot(n_symbols_list, np.full(len(n_symbols_list),snr), label="theory", ls="--")
plt.ylabel("est. SNR [dB]")
plt.xlabel("amount of symbols")
plt.grid()
plt.legend()

plt.show()