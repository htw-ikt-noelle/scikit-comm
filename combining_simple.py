import comm 
import numpy as np
import matplotlib.pyplot as plt
import copy

# signal parameters
n_dims = 4
mod_format = 'QAM'
mod_order = 4

# monte carlo parameters
block_size = 64
n_blocks = 100

# simulation parameters
mean_snr_dB = 15
n_bits = block_size*np.log2(mod_order)
comb_method = 'MRC'

# generate rayleigh-distributed linear SNR values
mean_snr_lin = 10**(mean_snr_dB/10)

rng = np.random.default_rng(seed=None)
rayleigh_snr_lin = rng.rayleigh(scale=mean_snr_lin,size=(n_blocks,n_dims))

# init arrays for SNR values
snr_post_comb_sim = np.zeros(shape=(n_blocks,),dtype='float')
snr_post_comb_theory = np.zeros(shape=(n_blocks,),dtype='float')
mean_snr_pre_comb = np.zeros(shape=(n_blocks,),dtype='float')

# monte carlo loop
for i in range(n_blocks):
    # standard signal generation block
    bit_seed = np.tile(np.random.randint(0,100,size=(1,)),n_dims)
    sig = comm.signal.Signal(n_dims=int(n_dims))
    sig.symbol_rate=5e9
    sig.sample_rate=15e9
    sig.generate_bits(n_bits=int(n_bits),seed=bit_seed)
    sig.generate_constellation(format=mod_format,order=mod_order)
    sig.mapper()
    sig.pulseshaper(upsampling=1,pulseshape=None)
    
    # set SNR per aperture
    sig.set_snr(snr_dB=(10*np.log10(rayleigh_snr_lin[i,:])).tolist(),seed=None)
    
    # combining
    sig_comb = comm.rx.combining(sig,comb_method=comb_method,snr_true=(10*np.log10(rayleigh_snr_lin[i,:])).tolist())
    # normalizing
    sig_comb.samples = sig_comb.samples[0] * (1/np.mean(np.abs(sig_comb.samples[0])))
    
    # snr post-combining
    # (Abi's 'data-aided' method for zero-mean signal and noise)
    # snr_post_comb_sim[i] = 10*np.log10(np.var(sig_comb.symbols[0])/np.var(sig_comb.samples[0]-sig_comb.symbols[0]))
    # M2M4 estimator
    snr_post_comb_sim[i] = 10*np.log10(comm.utils.estimate_SNR_m2m4(sig_comb.samples[0], sig_comb.constellation[0]))
    snr_post_comb_theory[i] = 10*np.log10(np.sum(rayleigh_snr_lin[i,:]))
    mean_snr_pre_comb[i] = 10*np.log10(np.mean(rayleigh_snr_lin[i,:]))

MRC_gain_sim = 10*np.log10(np.mean((10**(snr_post_comb_sim/10))/(10**(mean_snr_pre_comb/10))))
MRC_gain_theory = 10*np.log10(np.mean((10**(snr_post_comb_theory/10))/(10**(mean_snr_pre_comb/10))))
print('MRC gain sim: {}'.format(MRC_gain_sim))
print('MRC gain theory: {}'.format(MRC_gain_theory))

#### M2M4 SNR estimator test
# true_snr = 20
# sig_snr = comm.signal.Signal(n_dims=1)
# sig_snr.generate_bits(n_bits=2**6,seed=None)
# sig_snr.generate_constellation(format='QAM',order=4)
# sig_snr.mapper()
# sig_snr.pulseshaper(upsampling=1,pulseshape=None)
# sig_snr.set_snr(snr_dB=true_snr,seed=None)
# snr_estimate = comm.utils.estimate_SNR_m2m4(sig_snr.samples[0], sig_snr.constellation[0])
# print('true vs. estimated SNR: {} dB vs. {:.2f} dB'.format(true_snr,10*np.log10(snr_estimate)))