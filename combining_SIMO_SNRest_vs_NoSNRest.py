import comm 
import numpy as np
import matplotlib.pyplot as plt
import copy

# signal parameters
n_apertures = np.arange(2,11)
mod_format = 'QAM'
mod_order = 4
symb_rate = 5e9
sample_rate = 15e9
roll_off = 0.2
number_of_symbols = 10e3

# monte carlo parameters
n_blocks = 10

# simulation parameters
mean_snr_dB = 15
n_bits = comm.utils.symb_to_bits(number_of_symbols, mod_format, mod_order)
SNR_distribution = 'rayleigh'

# MC SNR arrays 
MRC_MC_gain = np.zeros(shape=(n_apertures.size,),dtype='float')
EGC_MC_gain = np.zeros(shape=(n_apertures.size,),dtype='float')
mrc_mc_snr_theory = np.zeros(shape=(n_apertures.size,),dtype='float')
egc_mc_snr_theory = np.zeros(shape=(n_apertures.size,),dtype='float')
MRC_MC_real = np.zeros(shape=(n_apertures.size,),dtype='float') 
EGC_MC_real = np.zeros(shape=(n_apertures.size,),dtype='float')
MRC_MC_real_est = np.zeros(shape=(n_apertures.size,),dtype='float')
EGC_MC_real_est = np.zeros(shape=(n_apertures.size,),dtype='float')
MRC_MC_gain_est = np.zeros(shape=(n_apertures.size,),dtype='float')
EGC_MC_gain_est = np.zeros(shape=(n_apertures.size,),dtype='float')

# Monte Carlo loop
for idx, n_dims in enumerate(n_apertures):
    # convert mean SNR value to linear 
    mean_snr_lin = 10**(mean_snr_dB/10)
    # init rng
    rng = np.random.default_rng(seed=None)
    #### rayleigh
    if SNR_distribution == 'rayleigh':
        # In a rayleigh channel, amplitudes of channel coefficients are rayleigh-
        # distributed. Linear SNRs are thus chi-square-distributed (with 2 degrees
        # of freedom), since the SNR is proportional to the square of the amplitude 
        # of the coefficients if the noise variance is assumed to be equal to 1.
        # The scale parameter of the rayleigh-distributed amplitudes must be set
        # to sqrt(mean_snr_lin/2), proof of which can be found in the accompanying 
        # PDF doc.
        rayleigh_amplitude = rng.rayleigh(scale=np.sqrt(mean_snr_lin/2),size=(n_blocks,n_dims))
    else:
        raise ValueError('SNR distribution not implemented at the moment.')
    
    # init arrays for SNR values
    # snr arrays
    snr_post_comb_sim = np.zeros(shape=(n_blocks,),dtype='float')
    mean_snr_pre_comb = np.zeros(shape=(n_blocks,),dtype='float')
    # theory curves
    mrc_snr_theory = np.zeros(shape=(n_blocks,),dtype='float')
    egc_snr_theory = np.zeros(shape=(n_blocks,),dtype='float')
    MRC_SNR = np.zeros(shape=(n_blocks,),dtype='float')
    EGC_SNR = np.zeros(shape=(n_blocks,),dtype='float')
    EGC_SNR_est = np.zeros(shape=(n_blocks,),dtype='float')
    MRC_SNR_est = np.zeros(shape=(n_blocks,),dtype='float')
    MRC_MC_gain_temp = np.zeros(shape=(n_blocks,),dtype='float')
    EGC_MC_gain_temp = np.zeros(shape=(n_blocks,),dtype='float')
    MRC_MC_gain_temp_est = np.zeros(shape=(n_blocks,),dtype='float')
    EGC_MC_gain_temp_est = np.zeros(shape=(n_blocks,),dtype='float')
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
        sig.pulseshaper(upsampling=sample_rate/symb_rate,pulseshape='rrc',roll_off=roll_off)
        
        #### CH signal block
        #### Abi's method: 
        sps = np.array(sig.sample_rate)/np.array(sig.symbol_rate)
        # normalize each dimension to have a mean power of 1
        for j in range(sig.n_dims):
            sig.samples[j] = sig.samples[j]/np.sqrt(np.mean(np.abs(sig.samples[j])**2))
        # scale samples with rayleigh_amplitude, add AWGN with mean=0,var=1, and power distributed equally among real and imaginary part
        for j in range(sig.n_dims):
            #sig.samples[j] = sig.samples[j] * rayleigh_amplitude[i,j] + np.sqrt(0.5)*(rng.standard_normal(size=(sig.samples[j].size,)) + 1j*rng.standard_normal(size=(sig.samples[j].size,)))
            n = (rng.standard_normal(size=(sig.samples[j].size,)) + 1j*rng.standard_normal(size=(sig.samples[j].size,)))
            sig.samples[j] = sig.samples[j] * rayleigh_amplitude[i,j]/np.sqrt(sps[j]) + np.sqrt(0.5)*n

        #### RX signal block
        # calulate SNR from sig_samples
        # define sig and noise range of signal
        sig_range = np.array([-sig.symbol_rate[0]/2-(sig.symbol_rate[0]/2*roll_off),sig.symbol_rate[0]/2+(sig.symbol_rate[0]/2*roll_off)])
        noise_range = np.array([-sig.symbol_rate[0]/2-(sig.symbol_rate[0]/2*roll_off)-1e9,-sig.symbol_rate[0]/2-(sig.symbol_rate[0]/2*roll_off),sig.symbol_rate[0]/2+(sig.symbol_rate[0]/2*roll_off),sig.symbol_rate[0]/2+1e9+(sig.symbol_rate[0]/2*roll_off)])

        temp_est_snr = np.zeros(sig.n_dims)
        for j in range(sig.n_dims):
            x_comb_EGC = np.fft.fftshift(np.fft.fftfreq(len(sig.samples[j]),1/sig.sample_rate[j]))
            y_comb_EGC = np.abs(np.fft.fftshift(np.fft.fft(sig.samples[j])))**2 
            temp_est_snr[j] = 10**(comm.utils.estimate_snr_spectrum(x_comb_EGC,y_comb_EGC,sig_range=sig_range, noise_range=noise_range, 
                                                      order=1,noise_bw=sig.symbol_rate[0],plotting=False)/10)
            
        # combining
        sig_comb_MRC = comm.rx.combining(sig,comb_method='MRC',snr_true=(10*np.log10(rayleigh_amplitude[i,:]**2)).tolist())
        sig_comb_EGC = comm.rx.combining(sig,comb_method='EGC',snr_true=(10*np.log10(rayleigh_amplitude[i,:]**2)).tolist())
        sig_comb_MRC_est = comm.rx.combining(sig,comb_method='MRC',snr_true=(10*np.log10(temp_est_snr)).tolist())
        sig_comb_EGC_est = comm.rx.combining(sig,comb_method='EGC',snr_true=(10*np.log10(temp_est_snr)).tolist())
    
        #### post-sim evaluation block
        # ! from here on out, linear SNR is equal to square of the rayleigh-distributed
        # channel coefficient amplitudes, since noise variance is equal to 1!
        
        # theory curves for linear SNR values after combining
        mrc_snr_theory[i] = np.sum(rayleigh_amplitude[i,:]**2) # ref.: pretty much every source
        egc_snr_theory[i] = (1/n_dims) * (np.sum(np.abs(np.sqrt(rayleigh_amplitude[i,:]**2)))**2) # ref.: H. Nuszkowski
        mean_snr_pre_comb = np.mean(rayleigh_amplitude[i,:]**2)

        # snr comb calc
        x_comb_EGC = np.fft.fftshift(np.fft.fftfreq(len(sig_comb_EGC.samples[0]),1/sig_comb_EGC.sample_rate[0]))
        y_comb_EGC = np.abs(np.fft.fftshift(np.fft.fft(sig_comb_EGC.samples[0])))**2 
        EGC_SNR[i] = 10**(comm.utils.estimate_snr_spectrum(x_comb_EGC,y_comb_EGC,sig_range=sig_range, noise_range=noise_range, 
                                                      order=1,noise_bw=sig_comb_EGC.symbol_rate[0],plotting=False)/10)

        x_comb_MRC = np.fft.fftshift(np.fft.fftfreq(len(sig_comb_MRC.samples[0]),1/sig_comb_MRC.sample_rate[0]))
        y_comb_MRC = np.abs(np.fft.fftshift(np.fft.fft(sig_comb_MRC.samples[0])))**2 
        MRC_SNR[i] = 10**(comm.utils.estimate_snr_spectrum(x_comb_MRC,y_comb_MRC,sig_range=sig_range, noise_range= noise_range, 
                                                      order=1,noise_bw=sig_comb_EGC.symbol_rate[0],plotting=False)/10)

        x_comb_EGC = np.fft.fftshift(np.fft.fftfreq(len(sig_comb_EGC_est.samples[0]),1/sig_comb_EGC_est.sample_rate[0]))
        y_comb_EGC = np.abs(np.fft.fftshift(np.fft.fft(sig_comb_EGC_est.samples[0])))**2 
        EGC_SNR_est[i] = 10**(comm.utils.estimate_snr_spectrum(x_comb_EGC,y_comb_EGC,sig_range=sig_range, noise_range=noise_range, 
                                                      order=1,noise_bw=sig_comb_EGC.symbol_rate[0],plotting=False)/10)

        x_comb_MRC = np.fft.fftshift(np.fft.fftfreq(len(sig_comb_MRC_est.samples[0]),1/sig_comb_MRC_est.sample_rate[0]))
        y_comb_MRC = np.abs(np.fft.fftshift(np.fft.fft(sig_comb_MRC_est.samples[0])))**2 
        MRC_SNR_est[i] = 10**(comm.utils.estimate_snr_spectrum(x_comb_MRC,y_comb_MRC,sig_range=sig_range, noise_range= noise_range, 
                                                      order=1,noise_bw=sig_comb_EGC.symbol_rate[0],plotting=False)/10)

        # mean SNR pre-combining
        mean_SNR_sim = np.mean(mean_snr_pre_comb)
        # SNR combining gain (simulation)
        MRC_MC_gain_temp[i] = 10*np.log10(MRC_SNR[i]) - 10*np.log10(mean_SNR_sim)
        EGC_MC_gain_temp[i] = 10*np.log10(EGC_SNR[i]) - 10*np.log10(mean_SNR_sim)
        MRC_MC_gain_temp_est[i] = 10*np.log10(MRC_SNR_est[i]) - 10*np.log10(mean_SNR_sim)
        EGC_MC_gain_temp_est[i] = 10*np.log10(EGC_SNR_est[i]) - 10*np.log10(mean_SNR_sim)

    # mean absolute SNR post-combining (theory)
    MRC_MC_gain[idx] = np.mean(MRC_MC_gain_temp)
    EGC_MC_gain[idx] = np.mean(EGC_MC_gain_temp)
    MRC_MC_gain_est[idx] = np.mean(MRC_MC_gain_temp_est)
    EGC_MC_gain_est[idx] = np.mean(EGC_MC_gain_temp_est) 
    MRC_MC_real[idx] = np.mean(MRC_SNR)
    EGC_MC_real[idx] = np.mean(EGC_SNR)
    MRC_MC_real_est[idx] = np.mean(MRC_SNR_est)
    EGC_MC_real_est[idx] = np.mean(EGC_SNR_est)
    mrc_mc_snr_theory[idx] = np.mean(mrc_snr_theory)
    egc_mc_snr_theory[idx] = np.mean(egc_snr_theory)
# end of MC loop

#### plotting
plt.close("all")
plt.rcParams['text.usetex'] = True

title_string = "\n Params: {}-{} with mean snr {}dB, {} destributed, {:.2E} runs, \n {:.2E} symbols, symbrate: {:.2E}, samplerate: {:.2E}".format(mod_order, mod_format, mean_snr_dB, SNR_distribution, n_blocks, number_of_symbols, symb_rate, sample_rate)

plt.figure(1)
plt.title("SNR with different combining techniques"+title_string)
plt.plot(n_apertures,10*np.log10(mrc_mc_snr_theory),label="MRC theory", color='salmon',linestyle='dashed')
plt.plot(n_apertures,10*np.log10(MRC_MC_real), label="MRC sim", color='r',marker="o", ls="", markersize=5)
plt.plot(n_apertures,10*np.log10(MRC_MC_real_est), label="MRC (est values) sim", color='firebrick',marker="^", ls="", markersize=5)
plt.plot(n_apertures,10*np.log10(egc_mc_snr_theory),label="EGC theory", color='steelblue',linestyle='dashed')
plt.plot(n_apertures,10*np.log10(EGC_MC_real), label="EGC sim", color='b',marker="o", ls="", markersize=5)
plt.plot(n_apertures,10*np.log10(EGC_MC_real_est), label="EGC (est values) sim", color='skyblue',marker="^", ls="", markersize=5)
plt.xlabel('Number of apertures')
plt.xticks(n_apertures)
plt.ylabel('SNR [dB]')
plt.legend()
plt.grid()
plt.savefig("1_baseline.png")

plt.figure(2)
plt.title("NORMALIZED SNR with different combining techniques"+title_string)
plt.plot(n_apertures,10*np.log10(n_apertures),label="MRC theory", color='salmon',linestyle='dashed')
plt.plot(n_apertures,MRC_MC_gain, label="MRC sim", color='r',marker="o", ls="", markersize=5)
plt.plot(n_apertures,MRC_MC_gain_est, label="MRC (est values) sim", color='firebrick',marker="^", ls="", markersize=5)
plt.plot(n_apertures,10*np.log10(1+(n_apertures-1)*np.pi/4),label="EGC theory", color='steelblue',linestyle='dashed')
plt.plot(n_apertures,EGC_MC_gain, label="EGC sim", color='b',marker="o", ls="", markersize=5)
plt.plot(n_apertures,EGC_MC_gain_est, label="EGC (est values) sim", color='skyblue',marker="^", ls="", markersize=5)
plt.xlabel('Number of apertures')
plt.xticks(n_apertures)
plt.ylabel('Normalized SNR [dB]')
plt.legend()
plt.grid()
plt.savefig("2_baseline.png")

plt.figure(3)
plt.title("Difference between combining with est SNR and without"+title_string)
plt.plot(n_apertures, np.zeros(len(n_apertures)), label="ideal", color="dimgray", ls="dashed")
plt.plot(n_apertures, 10*np.log10(MRC_MC_real)-10*np.log10(MRC_MC_real_est), label="MRC - MRC_est", color="maroon", marker="o")
plt.plot(n_apertures, 10*np.log10(EGC_MC_real)-10*np.log10(EGC_MC_real_est), label="EGC - EGC_est", color="dodgerblue", marker="o")
plt.xlabel('Number of apertures')
plt.xticks(n_apertures)
plt.ylabel('combined SNR difference [dB]')
plt.legend()
plt.grid()
plt.savefig("3_baseline.png")

plt.show()