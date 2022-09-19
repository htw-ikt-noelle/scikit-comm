import comm 
import numpy as np
import matplotlib.pyplot as plt
import copy

# signal parameters
n_apertures = np.arange(2,11)
mod_format = 'QAM'
mod_order = 4
symb_rate = 5e9
sample_rate = 5e9

# monte carlo parameters
block_size = int((2**12)*np.log2(mod_order)*(sample_rate/symb_rate))
n_blocks = 200

# simulation parameters
mean_snr_dB = 15
n_bits = block_size*np.log2(mod_order)/(sample_rate/symb_rate)
SNR_distribution = 'rayleigh'

# MC SNR arrays 
MRC_MC_gain = np.zeros(shape=(n_apertures.size,),dtype='float')
EGC_MC_gain = np.zeros(shape=(n_apertures.size,),dtype='float')
mrc_mc_snr_theory = np.zeros(shape=(n_apertures.size,),dtype='float')
egc_mc_snr_theory = np.zeros(shape=(n_apertures.size,),dtype='float')
MRC_MC_real = np.zeros(shape=(n_apertures.size,),dtype='float') 
EGC_MC_real = np.zeros(shape=(n_apertures.size,),dtype='float')
# Monte Carlo loop
for idx, n_dims in enumerate(n_apertures):
    # convert mean SNR value to linear 
    mean_snr_lin = 10**(mean_snr_dB/10)
    # init rng
    rng = np.random.default_rng(seed=None)
    #### rayleigh
    if SNR_distribution == 'rayleigh':
        # rayleigh_amplitude = rng.rayleigh(scale=mean_snr_lin/np.sqrt(np.pi/2),size=(n_blocks,n_dims))
        rayleigh_amplitude = rng.rayleigh(scale=np.sqrt(mean_snr_lin/2),size=(n_blocks,n_dims))
    elif SNR_distribution in ['normal','gaussian']:
        rayleigh_amplitude = rng.normal(loc=(10**(mean_snr_dB/10)),scale=1+(mean_snr_dB/10),size=(n_blocks,n_dims))
    elif SNR_distribution == 'uniform':
        rayleigh_amplitude = rng.uniform(low=0.25*(10**(mean_snr_dB/10)),high=4*(10**(mean_snr_dB/10)),size=(n_blocks,n_dims))
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
    MRC_SNR = np.zeros(shape=(n_blocks,),dtype='float')
    EGC_SNR = np.zeros(shape=(n_blocks,),dtype='float')
    MRC_MC_gain_temp = np.zeros(shape=(n_blocks,),dtype='float')
    EGC_MC_gain_temp = np.zeros(shape=(n_blocks,),dtype='float')
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
        # scale samples with np.sqrt(rayleigh_amplitude), add AWGN with mean=0,std=1 (rng.standard_normal)
        for j in range(sig.n_dims):
            sig.samples[j] = sig.samples[j] * rayleigh_amplitude[i,j] + np.sqrt(0.5)*(rng.standard_normal(size=(sig.samples[j].size,)) + 1j*rng.standard_normal(size=(sig.samples[j].size,)))
        #### set SNR method
        # set SNR per aperture
        # sig.set_snr(snr_dB=(10*np.log10(rayleigh_amplitude[i,:])).tolist(),seed=None)
        
        #### RX signal block
        # combining
        sig_comb_MRC = comm.rx.combining(sig,comb_method='MRC',snr_true=(10*np.log10(rayleigh_amplitude[i,:]**2)).tolist())
        sig_comb_EGC = comm.rx.combining(sig,comb_method='EGC',snr_true=(10*np.log10(rayleigh_amplitude[i,:]**2)).tolist())
    
        #### post-sim evaluation block
        # theory curves for linear SNR values after combining
        mrc_snr_theory[i] = np.sum(rayleigh_amplitude[i,:]**2) # ref.: pretty much every source
        egc_snr_theory[i] = (1/n_dims) * (np.sum(np.abs(np.sqrt(rayleigh_amplitude[i,:]**2)))**2) # ref.: H. Nuszkowski
        # M2M4 estimator
        snr_post_comb_theory[i] = 10*np.log10(np.sum(rayleigh_amplitude[i,:]))
        mean_snr_pre_comb = np.mean(rayleigh_amplitude[i,:]**2)
    
        # SNR post-combining
        MRC_SNR[i] = comm.utils.estimate_SNR_m2m4(sig_comb_MRC.samples[0], sig_comb_MRC.constellation[0])
        EGC_SNR[i] = comm.utils.estimate_SNR_m2m4(sig_comb_EGC.samples[0], sig_comb_MRC.constellation[0])

        # mean SNR pre-combining
        mean_SNR_sim = np.mean(mean_snr_pre_comb)
        # SNR combining gain (simulation)
        MRC_MC_gain_temp[i] = 10*np.log10(MRC_SNR[i]) - 10*np.log10(mean_SNR_sim)
        EGC_MC_gain_temp[i] = 10*np.log10(EGC_SNR[i]) - 10*np.log10(mean_SNR_sim)

    # mean absolute SNR post-combining (theory)
    MRC_MC_gain[idx] = np.mean(MRC_MC_gain_temp)
    EGC_MC_gain[idx] = np.mean(EGC_MC_gain_temp) 
    MRC_MC_real[idx] = np.mean(MRC_SNR)
    EGC_MC_real[idx] = np.mean(EGC_SNR)
    mrc_mc_snr_theory[idx] = np.mean(mrc_snr_theory)
    egc_mc_snr_theory[idx] = np.mean(egc_snr_theory)
# end of MC loop

#### plotting
plt.close("all")
plt.rcParams['text.usetex'] = True

plt.figure(1)
plt.title('SNR with different combining techniques, '+str(mod_order)+'-'+str(mod_format)+"\n mean snr: "+str(mean_snr_dB)+"dB, "+str(SNR_distribution)+" destributed")
plt.plot(n_apertures,10*np.log10(mrc_mc_snr_theory),label="MRC theory", color='salmon',linestyle='dashed')
plt.plot(n_apertures,10*np.log10(MRC_MC_real), label="MRC real", color='r',marker="o", ls="", markersize=5)
plt.plot(n_apertures,10*np.log10(egc_mc_snr_theory),label="EGC theory", color='steelblue',linestyle='dashed')
plt.plot(n_apertures,10*np.log10(EGC_MC_real), label="EGC real", color='b',marker="o", ls="", markersize=5)
plt.xlabel('Number of apertures')
plt.xticks(n_apertures)
plt.ylabel('SNR [dB]')
plt.legend()
plt.grid()

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
textstr = '\n'.join(("Theory MRC:", r"$$\gamma_{R} = \sum_{l=1}^L \gamma_{k}$$", "Theory EGC:", r"$$\gamma_{R} = \frac{1}{L} \cdot \left( \sum_{l=1}^L |\alpha_l| \right )^2$$"))
plt.text(0.6, 0.4, textstr, fontsize=10, verticalalignment='top', bbox=props, transform=plt.gca().transAxes)

plt.figure(2)
plt.title('NORMALIZED SNR with different combining techniques, '+str(mod_order)+'-'+str(mod_format)+"\n mean snr: "+str(mean_snr_dB)+"dB, "+str(SNR_distribution)+" destributed")
plt.plot(n_apertures,10*np.log10(n_apertures),label="MRC theory", color='salmon',linestyle='dashed')
plt.plot(n_apertures,MRC_MC_gain, label="MRC real", color='r',marker="o", ls="", markersize=5)
plt.plot(n_apertures,10*np.log10(1+(n_apertures-1)*np.pi/4),label="EGC theory", color='steelblue',linestyle='dashed')
plt.plot(n_apertures,EGC_MC_gain, label="EGC real", color='b',marker="o", ls="", markersize=5)
plt.xlabel('Number of apertures')
plt.xticks(n_apertures)
plt.ylabel('Normalized SNR [dB]')
plt.legend()
plt.grid()

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
textstr = '\n'.join(("Theory MRC:", r"$$\gamma_{R} = L$$", "Theory EGC:", r"$$\gamma_{R} = 1+(L-1) \cdot \frac{\pi}{4}$$"))
plt.text(0.6, 0.4, textstr, fontsize=10, verticalalignment='top', bbox=props, transform=plt.gca().transAxes)


plt.show()


# plt.plot(n_apertures,np.abs(MRC_MC_gain-MRC_gain_theory),color='salmon',linestyle='dashed')
# plt.plot(n_apertures,np.abs(EGC_MC_gain-EGC_gain_theory),color='steelblue',linestyle='dashed')
# # attributes
# plt.grid()
# plt.title('Simulation-Theory gap, '+str(mod_order)+'-'+str(mod_format))
# plt.xlabel('Number of apertures')
# plt.xticks(n_apertures)
# plt.ylabel('absolute error of normalized SNR [dB]')
# plt.legend(('MRC sim-theory gap','EGC sim-theory gap'))
# plt.show()


# # theory curves
# MRC_gain_theory = 10*np.log10(mrc_mc_snr_theory) - 10*np.log10(mean_SNR_sim) # for n_apertures >> 2: MRC_gain_theory = 10*np.log10(n_apertures)
# EGC_gain_theory = 10*np.log10(egc_mc_snr_theory) - 10*np.log10(mean_SNR_sim) # for n_apertures >> 2: EGC_gain_theory = 10*np.log10(n_apertures*np.pi/4)
# # close all open plots
# plt.close('all')
# # figures
# plt.figure(1)
# # curves
# plt.plot(n_apertures,MRC_gain_theory,color='salmon',linestyle='dashed')
# plt.scatter(n_apertures,MRC_MC_gain,color='r')
# plt.plot(n_apertures,EGC_gain_theory,color='steelblue',linestyle='dashed')
# plt.scatter(n_apertures,EGC_MC_gain,color='b')
# # attributes
# plt.grid()
# plt.title('Normalized SNR with different combining techniques, '+str(mod_order)+'-'+str(mod_format))
# plt.xlabel('Number of apertures')
# plt.xticks(n_apertures)
# plt.ylabel('Normalized SNR [dB]')
# plt.legend(('MRC, theory','MRC, simulation','EGC, theory','EGC, simulation'))
# plt.show()

# plt.figure(2)
# # curves
# plt.plot(n_apertures,np.abs(MRC_MC_gain-MRC_gain_theory),color='salmon',linestyle='dashed')
# plt.plot(n_apertures,np.abs(EGC_MC_gain-EGC_gain_theory),color='steelblue',linestyle='dashed')
# # attributes
# plt.grid()
# plt.title('Simulation-Theory gap, '+str(mod_order)+'-'+str(mod_format))
# plt.xlabel('Number of apertures')
# plt.xticks(n_apertures)
# plt.ylabel('absolute error of normalized SNR [dB]')
# plt.legend(('MRC sim-theory gap','EGC sim-theory gap'))
# plt.show()

# plt.figure(3)
# # curves
# plt.plot(n_apertures,(10*np.log10(n_apertures))-MRC_gain_theory,color='salmon',linestyle='dashed')
# plt.plot(n_apertures,(10*np.log10(n_apertures*np.pi/4))-EGC_gain_theory,color='steelblue',linestyle='dashed')
# # attributes
# plt.grid()
# plt.title('approximation error of normalized SNR, '+str(mod_order)+'-'+str(mod_format))
# plt.xlabel('Number of apertures')
# plt.xticks(n_apertures)
# plt.ylabel('absolute error ormalized SNR [dB]')
# plt.legend(('MRC theory-approximation gap','EGC theory-approximation gap'))
# plt.show()