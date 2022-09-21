# test file for determining spectral SNR estimation accuracy

import numpy as np
import matplotlib.pyplot as plt
import comm

# signal attributes
n_dims = 1
mod_format = 'QAM'
mod_order = 4
symb_rate = 5e9
samp_rate = 15e9
n_bits = 2**14
n_samples = (n_bits/np.log2(mod_order))*(samp_rate/symb_rate)
roll_off = np.tile(np.array([0.2]),n_dims)

# sim attributes
mc_runs = 100
snr_dB = np.arange(0,21)
plotting = False

rmse_norm_estimation = np.zeros((snr_dB.size,),dtype='float')
# snr loop
for idx, snr in enumerate(snr_dB):
    rmse_norm_estimation_mc = np.zeros((mc_runs,),dtype='float')
    # mc loop
    for j in range(mc_runs):
        # gen signal
        sig = comm.signal.Signal(n_dims=n_dims)
        sig.symbol_rate = symb_rate
        sig.generate_bits(n_bits=n_bits,seed=None)
        sig.generate_constellation(format=mod_format,order=mod_order)
        sig.mapper()
        sig.pulseshaper(upsampling=(samp_rate/symb_rate),pulseshape='rrc',roll_off=0.2)
        # set SNR
        # rng = np.random.default_rng(seed=None)
        # snr_lin = rng.normal(loc=10**(mean_snr_dB/10),scale=1,shape=(n_dims,))
        sig.set_snr(snr_dB=snr.tolist(),seed=None)
        # estimate snr
        snr_est_dB = comm.utils.est_snr_spec_wrapper(sig,roll_off,plotting=plotting)
        # MSE
        # convert snr to lin
        snr_est_lin = 10**(snr_est_dB/10)
        snr_lin = 10**(snr/10)
        # normalized rmse
        rmse_norm_estimation_mc[j] = np.sqrt(np.mean((snr_est_lin-snr_lin)**2))/snr_lin
    rmse_norm_estimation[idx] = np.mean(rmse_norm_estimation_mc)
    
# plot
plt.figure(1)
plt.semilogy(snr_dB,rmse_norm_estimation,color='r')
plt.xlabel('true SNR [dB]')
plt.xticks(snr_dB[::2])
plt.ylabel('normalized RMSE')
plt.title('Normalized RMSE of spectral SNR estimation, {} across {} runs.'.format(str(str(mod_order)+mod_format),n_dims*mc_runs))
plt.grid()
plt.show()