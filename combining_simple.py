import comm 
import numpy as np
import matplotlib.pyplot as plt
import copy

# signal parameters
n_apertures = np.arange(2,21)
mod_format = 'QAM'
mod_order = 4
symb_rate = 5e9
sample_rate = 5e9

# monte carlo parameters
block_size = int((2**12)*np.log2(mod_order)*(sample_rate/symb_rate))
n_blocks = 100

# simulation parameters
mean_snr_dB = 15
n_bits = block_size*np.log2(mod_order)/(sample_rate/symb_rate)
SNR_distribution = 'rayleigh'

# MC SNR arrays 
MRC_MC_gain = np.zeros(shape=(n_apertures.size,),dtype='float')
EGC_MC_gain = np.zeros(shape=(n_apertures.size,),dtype='float')
mrc_mc_snr_theory = np.zeros(shape=(n_apertures.size,),dtype='float')
egc_mc_snr_theory = np.zeros(shape=(n_apertures.size,),dtype='float')
# Monte Carlo loop
for idx, n_dims in enumerate(n_apertures):
    # convert mean SNR value to linear 
    mean_snr_lin = 10**(mean_snr_dB/10)
    # init rng
    rng = np.random.default_rng(seed=None)
    #### rayleigh
    if SNR_distribution == 'rayleigh':
        rayleigh_snr_lin = rng.rayleigh(scale=mean_snr_lin/np.sqrt(np.pi/2),size=(n_blocks,n_dims))
    elif SNR_distribution in ['normal','gaussian']:
        rayleigh_snr_lin = rng.normal(loc=(10**(mean_snr_dB/10)),scale=1+(mean_snr_dB/10),size=(n_blocks,n_dims))
    elif SNR_distribution == 'uniform':
        rayleigh_snr_lin = rng.uniform(low=0.25*(10**(mean_snr_dB/10)),high=4*(10**(mean_snr_dB/10)),size=(n_blocks,n_dims))
    else:
        raise ValueError('SNR distribution not implemented at the moment.')
    
    # init arrays for SNR values
    # snr arrays
    snr_post_comb_sim = np.zeros(shape=(n_blocks,),dtype='float')
    snr_post_comb_theory = np.zeros(shape=(n_blocks,),dtype='float')
    mean_snr_pre_comb = np.zeros(shape=(n_blocks,),dtype='float')
    # theory curves
    mrc_snr_theory = np.zeros(shape=(n_blocks,),dtype='float')
    egc_snr_theory = np.zeros(shape=(n_blocks,),dtype='float')
    # arrays for concatenating combined MC sample blocks
    samples_MRC = np.zeros(shape=(int((sample_rate/symb_rate)*n_blocks*block_size),),dtype='complex')
    samples_EGC = np.zeros(shape=(int((sample_rate/symb_rate)*n_blocks*block_size),),dtype='complex')
    # monte carlo loop
    for i in range(n_blocks):
        #### TX signal block
        bit_seed = np.tile(np.random.randint(0,100,size=(1,)),n_dims)
        sig = comm.signal.Signal(n_dims=int(n_dims))
        sig.symbol_rate= symb_rate
        sig.sample_rate= sample_rate
        sig.generate_bits(n_bits=int(n_bits),seed=bit_seed)
        sig.generate_constellation(format=mod_format,order=mod_order)
        sig.mapper()
        # only integer upsampling factors for now
        sig.pulseshaper(upsampling=int(sample_rate/symb_rate),pulseshape='rrc',roll_off=0.2)
        
        #### CH signal block
        #### Abi's method: 
        # normalize each dimension to have a mean power of 1
        for j in range(sig.n_dims):
            sig.samples[j] = sig.samples[j] / np.sqrt(np.mean(np.abs(sig.samples[j])**2))
        # scale samples with np.sqrt(rayleigh_snr_lin), add AWGN with mean=0,std=1 (rng.standard_normal)
        for j in range(sig.n_dims):
            sig.samples[j] = sig.samples[j] * np.sqrt(rayleigh_snr_lin[i,j]) + np.sqrt(0.5)*(rng.standard_normal(size=(sig.samples[j].size,)) + 1j*rng.standard_normal(size=(sig.samples[j].size,)))
        #### set SNR method
        # set SNR per aperture
        # sig.set_snr(snr_dB=(10*np.log10(rayleigh_snr_lin[i,:])).tolist(),seed=None)
        
        #### RX signal block
        # combining
        sig_comb_MRC = comm.rx.combining(sig,comb_method='MRC',snr_true=(10*np.log10(rayleigh_snr_lin[i,:])).tolist())
        sig_comb_EGC = comm.rx.combining(sig,comb_method='EGC',snr_true=(10*np.log10(rayleigh_snr_lin[i,:])).tolist())
    
        #### post-sim evaluation block
        # theory curves for linear SNR values after combining
        mrc_snr_theory[i] = np.sum(rayleigh_snr_lin[i,:]) # ref.: pretty much every source
        egc_snr_theory[i] = (1/n_dims) * (np.sum(np.abs(np.sqrt(rayleigh_snr_lin[i,:])))**2) # ref.: H. Nuszkowski
        # M2M4 estimator
        snr_post_comb_theory[i] = 10*np.log10(np.sum(rayleigh_snr_lin[i,:]))
        mean_snr_pre_comb[i] = np.mean(rayleigh_snr_lin[i,:])
        # write combined samples to array
        samples_MRC[int(i*block_size):int((i+1)*block_size)] = sig_comb_MRC.samples[0]
        samples_EGC[int(i*block_size):int((i+1)*block_size)] = sig_comb_EGC.samples[0]
    
    # SNR post-combining
    MRC_SNR = comm.utils.estimate_SNR_m2m4(samples_MRC, sig_comb_MRC.constellation[0])
    EGC_SNR = comm.utils.estimate_SNR_m2m4(samples_EGC, sig_comb_MRC.constellation[0])
    # mean SNR pre-combining
    mean_SNR_sim = np.mean(mean_snr_pre_comb)
    # SNR combining gain (simulation)
    MRC_MC_gain[idx] = 10*np.log10(MRC_SNR) - 10*np.log10(mean_SNR_sim)
    EGC_MC_gain[idx] = 10*np.log10(EGC_SNR) - 10*np.log10(mean_SNR_sim)
    # mean absolute SNR post-combining (theory)
    mrc_mc_snr_theory[idx] = np.mean(mrc_snr_theory)
    egc_mc_snr_theory[idx] = np.mean(egc_snr_theory)
# end of MC loop

#### plotting
# theory curves
MRC_gain_theory = 10*np.log10(mrc_mc_snr_theory) - 10*np.log10(mean_SNR_sim) # for n_apertures >> 2: MRC_gain_theory = 10*np.log10(n_apertures)
EGC_gain_theory = 10*np.log10(egc_mc_snr_theory) - 10*np.log10(mean_SNR_sim) # for n_apertures >> 2: EGC_gain_theory = 10*np.log10(n_apertures*np.pi/4)
# close all open plots
plt.close('all')
# figures
plt.figure(1)
# curves
plt.plot(n_apertures,MRC_gain_theory,color='salmon',linestyle='dashed')
plt.scatter(n_apertures,MRC_MC_gain,color='r')
plt.plot(n_apertures,EGC_gain_theory,color='steelblue',linestyle='dashed')
plt.scatter(n_apertures,EGC_MC_gain,color='b')
# attributes
plt.grid()
plt.title('Normalized SNR with different combining techniques, '+str(mod_order)+'-'+str(mod_format))
plt.xlabel('Number of apertures')
plt.xticks(n_apertures)
plt.ylabel('Normalized SNR [dB]')
plt.legend(('MRC, theory','MRC, simulation','EGC, theory','EGC, simulation'))
plt.show()

plt.figure(2)
# curves
plt.plot(n_apertures,np.abs(MRC_MC_gain-MRC_gain_theory),color='salmon',linestyle='dashed')
plt.plot(n_apertures,np.abs(EGC_MC_gain-EGC_gain_theory),color='steelblue',linestyle='dashed')
# attributes
plt.grid()
plt.title('Simulation-Theory gap, '+str(mod_order)+'-'+str(mod_format))
plt.xlabel('Number of apertures')
plt.xticks(n_apertures)
plt.ylabel('absolute error of normalized SNR [dB]')
plt.legend(('MRC sim-theory gap','EGC sim-theory gap'))
plt.show()

plt.figure(3)
# curves
plt.plot(n_apertures,(10*np.log10(n_apertures))-MRC_gain_theory,color='salmon',linestyle='dashed')
plt.plot(n_apertures,(10*np.log10(n_apertures*np.pi/4))-EGC_gain_theory,color='steelblue',linestyle='dashed')
# attributes
plt.grid()
plt.title('approximation error of normalized SNR, '+str(mod_order)+'-'+str(mod_format))
plt.xlabel('Number of apertures')
plt.xticks(n_apertures)
plt.ylabel('absolute error ormalized SNR [dB]')
plt.legend(('MRC theory-approximation gap','EGC theory-approximation gap'))
plt.show()