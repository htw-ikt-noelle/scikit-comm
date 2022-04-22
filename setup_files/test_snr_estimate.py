import numpy as np 
import matplotlib.pyplot as plt 
import comm as comm 
 
# global 
mod_format = 'QAM' 
mod_order = 16 
TX_UPSAMPLE_FACTOR = 5 
# SNR = 10 #in dB 
SNR = np.arange(0,30,dtype='float') 
# number of runs in the Monte Carlo trial 
MC_runs = 1
 
# init vectors of SNR estimates 
snr_vec = np.zeros((len(SNR),MC_runs),dtype='float') 
evm_vec = np.zeros((len(SNR),MC_runs),dtype='float') 
evm_vec_opt = np.zeros((len(SNR),MC_runs),dtype='float') 
evm_ref_vec = np.zeros((len(SNR),MC_runs),dtype='float') 
evm_ref_vec_opt = np.zeros((len(SNR),MC_runs),dtype='float') 

plt.close()
 
# SNR loop 
for i in range(len(SNR)): 
 
    # MC loop 
    for j in range(MC_runs): 
         
        # construct signal 
        sig_tx = comm.signal.Signal(n_dims=1) 
        sig_tx.ndims = 1 
        sig_tx.symbol_rate = 50e6  
        # sig_tx.sample_rate = 50e6  
         
        # generate bits 
        sig_tx.generate_bits(n_bits=int(np.log2(mod_order))*(2**15)) 
         
        # set constellation (modulation format) 
        sig_tx.generate_constellation(format=mod_format,order=mod_order) 
         
        # create symbols 
        sig_tx.mapper() 
         
        # # upsampling and pulseshaping 
        ROLL_OFF = 0.2 
        sig_tx.pulseshaper(upsampling=TX_UPSAMPLE_FACTOR, pulseshape='rrc', roll_off=[ROLL_OFF]) 
        # sig_tx.samples = sig_tx.symbols[0] 
        
        # add complex noise 
        sig_tx.set_snr(snr_dB=SNR[i]) 
         
        # # RX matched filter 
        sig_tx.samples = comm.filters.raised_cosine_filter(sig_tx.samples[0],root_raised=True,roll_off=ROLL_OFF, 
                                                            symbol_rate=sig_tx.symbol_rate[0], 
                                                            sample_rate=sig_tx.sample_rate[0]) 
         
         
        # # downsample to 1 sps 
        sig_tx.samples = sig_tx.samples[0][::TX_UPSAMPLE_FACTOR]         
         
        # calc EVM 
        evm = comm.rx.calc_evm(sig_tx, norm='rms', method='blind', opt=False, dimension=-1) 
        evm_vec[i,j] = evm 
        
        # evm = comm.rx.calc_evm(sig_tx, norm='rms', method='blind', opt=True, dimension=-1) 
        # evm_vec_opt[i,j] = evm 
        
        evm_ref = comm.rx.calc_evm(sig_tx, norm='rms', method='data_aided', opt=False, dimension=-1) 
        evm_ref_vec[i,j] = evm_ref
        
        evm_ref = comm.rx.calc_evm(sig_tx, norm='rms', method='data_aided', opt=True, dimension=-1) 
        evm_ref_vec_opt[i,j] = evm_ref
        
        
         
#%% visualization  
 
evm_mean = np.nanmean(evm_vec,axis=1) 
# evm_mean_opt = np.nanmean(evm_vec_opt,axis=1) 
evm_ref_mean = np.nanmean(evm_ref_vec,axis=1) 
evm_ref_mean_opt = np.nanmean(evm_ref_vec_opt,axis=1) 

# from EVM to SNR: "On the Extended Relationships Among EVM, BER and SNR as  
# Performance Metrics" https://doi.org/10.1109/ICECE.2006.355657 
snr_evm_mean = 10*np.log10(1/evm_mean**2) 
snr_evm_std = np.nanstd(10*np.log10(1/evm_vec**2),axis=1) 
# snr_evm_mean_opt = 10*np.log10(1/evm_mean_opt**2) 
snr_evm_ref_mean = 10*np.log10(1/evm_ref_mean**2) 
snr_evm_ref_std = np.nanstd(10*np.log10(1/evm_ref_vec**2),axis=1) 
snr_evm_ref_mean_opt = 10*np.log10(1/evm_ref_mean_opt**2) 
 
plt.figure(1) 
# plot ideal/unbiased SNR curve for reference 
plt.plot(SNR,SNR,color='b') 
# plot mean 
plt.plot(SNR,snr_evm_mean,color='k') 
# plt.plot(SNR,snr_evm_mean_opt,color='g') 

plt.plot(SNR,snr_evm_ref_mean,color='r') 
plt.plot(SNR,snr_evm_ref_mean_opt,color='y') 

# # # plot std deviation as grayscale 
# # plt.fill_between(SNR,snr_evm_mean+snr_evm_std,snr_evm_mean-snr_evm_std,color='grey') 
# # plt.fill_between(SNR,snr_evm_ref_mean+snr_evm_ref_std,snr_evm_ref_mean-snr_evm_ref_std,color='grey') 
# # plt.legend(('linear','est. SNR from EVM','est. SNR from EVM w ref symb','standard deviation across ' + str(MC_runs) + ' runs', 'std. EVM ref')) 
# plt.xticks(SNR[::2]) 
# plt.yticks(SNR[::2]) 
# plt.xlim(np.min(SNR), np.max(SNR)) 
# plt.ylim(np.min(SNR), np.max(SNR)) 
 
# plt.title('actual vs. estimated mean SNR') 
# plt.grid() 
 
# plt.figure(2) 
# plt.plot(SNR,np.abs(SNR-snr_mean),color='r') 
# plt.plot(SNR,np.abs(SNR-snr_evm_mean),color='b') 
# plt.plot(SNR,np.full_like(SNR,1),color='grey') 
# plt.title('mean SNR estimation error') 
# plt.legend(('absolute mean SNR estimation error in dB','absolute mean EVM-based SNR estimation error in dB','1 dB error threshold')) 
# plt.grid() 