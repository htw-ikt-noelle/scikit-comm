import numpy as np 
import matplotlib.pyplot as plt 
import comm as comm 
 
# global 
N = 100000
rng = np.random.default_rng(seed=1)

samples_orig =  ((rng.standard_normal(size=N) > 0) - 0.5) * 2

# print(np.mean(np.abs(samples)**2))

# add complex noise 
# sig_tx.set_snr(snr_dB=SNR) 

noise = np.sqrt(0.1) * ((rng.standard_normal(size=samples_orig.shape[0])))
# print(np.mean(np.abs(noise)**2))
samples = samples_orig + noise

plt.plot(samples, '.')

# calc EVM 
for scaling_fac in np.arange(0.7,1.1,0.1):
    symbols_norm = scaling_fac * samples
    
    # plt.hist(np.abs(symbols_norm)**2, bins=100)
    # plt.xlim([0, 8])
    # plt.show()

    
    error = symbols_norm - samples_orig
    
    # plt.hist(np.abs(error)**2, bins=100)
    # plt.xlim([0, 2])
    # plt.show()
    
    mean_error = np.sqrt(np.mean(np.abs(error)**2)) 
    
    # print(mean_error**2)
    
    print('scale: {:.1f}, mean: {:.3f}, std: {:.3f}, mean error: {:.3f}'.format(scaling_fac, np.mean(error), np.std(error), mean_error))

 
