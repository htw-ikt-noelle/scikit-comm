import numpy as np 
import matplotlib.pyplot as plt 
import comm as comm 
 
# global 
mod_format = 'QAM' 
mod_order = 4 
SNR = 10
        
# construct signal 
sig_tx = comm.signal.Signal(n_dims=1) 
sig_tx.symbol_rate = 1
sig_tx.sample_rate = 1
 
# generate bits 
sig_tx.generate_bits(n_bits=int(np.log2(mod_order))*(2**15)) 
 
# set constellation (modulation format) 
sig_tx.generate_constellation(format=mod_format,order=mod_order) 
 
# create symbols 
sig_tx.mapper() 
sig_tx.samples = sig_tx.symbols[0] 

print(np.mean(np.abs(sig_tx.samples[0])**2))

# add complex noise 
# sig_tx.set_snr(snr_dB=SNR) 
rng = np.random.default_rng(seed=1)
noise = np.sqrt(0.1) * ((rng.standard_normal(size=sig_tx.samples[0].shape[0]) +1j* rng.standard_normal(size=sig_tx.samples[0].shape[0])))
print(np.mean(np.abs(noise)**2))
sig_tx.samples = sig_tx.samples[0] + noise

# calc EVM 
for scaling_fac in np.arange(0.5,1.1,0.1):
    symbols_norm = scaling_fac * sig_tx.samples[0]
    
    # plt.hist(np.abs(symbols_norm)**2, bins=100)
    # plt.xlim([0, 8])
    # plt.show()

    
    error = symbols_norm - sig_tx.symbols[0]
    
    # plt.hist(np.abs(error)**2, bins=100)
    # plt.xlim([0, 2])
    # plt.show()
    
    mean_error = np.sqrt(np.mean(np.abs(error)**2)) 
    
    print(mean_error**2)
    
    evm_norm_ref = np.sqrt(np.mean(np.abs(sig_tx.constellation[0])**2))
    evm = mean_error / evm_norm_ref
    # evm = mean_error
    
    # comm.visualizer.plot_constellation(symbols_norm, hist=True)
        
    print('real SNR: {:.1f} scale: {:.1f}, EVM: {:.1f}%, est. SNR: {:.1f}'.format(SNR, scaling_fac, evm*100, 10*np.log10(1/evm**2)))

 
