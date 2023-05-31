import copy

import numpy as np
import scipy.signal as ssignal
import matplotlib.pyplot as plt
from scipy.stats import rayleigh

import skcomm as skc

def tx_sim(config):
    """
    TX simulation
    """
    # contruct signal
    sig_tx = skc.signal.Signal(n_dims=config.n_dims)
    sig_tx.symbol_rate = config.symbol_rate

    tx_upsample_factor = config.DAC_sr / config.symbol_rate

    # generate bits
    sig_tx.generate_bits(n_bits=int(config.n_symbols*np.log2(config.modulation_order)), seed=config.TX_seed_bits)

    # set constellation (modulation format)
    sig_tx.generate_constellation(format=config.modulation_format, order=config.modulation_order)

    # create symbols
    sig_tx.mapper()

    # upsampling and pulseshaping
    sig_tx.pulseshaper(upsampling=tx_upsample_factor, pulseshape=config.pulseshape, roll_off=config.roll_off)

    return sig_tx

def ch_sim(sig_ch, config, SNR_arr):
    """
    CH simulation
    """
    # add snr
    for i in range(config.n_dims):
        sig_ch.samples[i] = skc.channel.set_snr(sig_ch.samples[i], snr_dB=SNR_arr[i], sps=sig_ch.sample_rate[0]/sig_ch.symbol_rate[0], seed=None)

    return sig_ch

def rx_sim(sig_rx, config, which_dim=0):
    """
    RX simulation (for 1 channel of samples)
    """
    # set which_dim param to int
    which_dim = int(which_dim)

    # if there is only one dimension, there are some indexing errors. 
    try:
        _ = sig_rx.sample_rate[which_dim]
    except TypeError:
        for key in vars(sig_rx):
            if key == '_n_dims' or key == '_samples':
                pass
            else:
                vars(sig_rx)[key] = [vars(sig_rx)[key]]

    sps_new = 2
    sps = sig_rx.sample_rate[which_dim]/sig_rx.symbol_rate[which_dim]
    new_length = int(sig_rx.samples[which_dim].size/sps*sps_new)
    sig_rx.samples[which_dim] = ssignal.resample(sig_rx.samples[which_dim], new_length, window='boxcar')
    sig_rx.sample_rate[which_dim] = sps_new*sig_rx.symbol_rate[which_dim]

    ## expected FOE?
    fo = 0

    # # estimate SNR
    sig_range = np.asarray([-1, 1])*sig_rx.symbol_rate[which_dim]/2*(1+config.roll_off[which_dim]) + fo 
    noise_range = np.asarray([-1.1, -1.05, 1.05, 1.1]) * sig_rx.symbol_rate[which_dim]/2 * (1+config.roll_off[which_dim]) + fo

    spec = np.abs(np.fft.fftshift(np.fft.fft(sig_rx.samples[which_dim])))**2
    freq = np.fft.fftshift(np.fft.fftfreq(sig_rx.samples[which_dim].size, 1/sig_rx.sample_rate[which_dim]))
    snr = skc.utils.estimate_snr_spectrum(freq, spec, sig_range=sig_range, 
                                        noise_range=noise_range, order=1, 
                                        noise_bw=sig_rx.symbol_rate[which_dim], 
                                        scaling='lin', plotting=False)
    
    return snr, sig_rx

if __name__ == "__main__":
    config = skc.simulations.config_sim()

    Nr_list = [1,2,4] # set Nr of apertures we want to simulate

    # set config data
    config.symbol_rate = 10e9
    config.DAC_sr = 20e9
    config.n_dims = int(np.max(Nr_list))  # set num of max. Antenna to simulate
    config.TX_seed_bits = [1] * config.n_dims
    config.roll_off = [0.1] * config.n_dims
    config.crop_factor_matched_filter = 0
    config.n_symbols = 5e5

    # do TX simulation
    sig_tx = tx_sim(config)

    # check if sig_tx hast two dimension with same samples
    print("equality between samples 0 and 1: {}".format((sig_tx.samples[0] == sig_tx.samples[1]).any()))

    # set amount of MC-Runs
    mc_runs = 10
    save_data = True
    # installing rng & generate some random SNR data
    rng = np.random.default_rng()
    #EsN0_dB = rng.integers(0,20,mc_runs)
    #SNR_arr = rng.integers(0,20,size=(config.n_dims,mc_runs))

    # For Rayleigh distributed Amplitudes:
    rayleigh_amp = np.abs(rng.normal(0,np.sqrt(0.5), size=(config.n_dims,mc_runs))+1j*rng.normal(0,np.sqrt(0.5), size=(config.n_dims,mc_runs)))
    mean_SNR_values = rng.integers(0,20,size=mc_runs)
    SNR_arr = np.full(rayleigh_amp.shape, mean_SNR_values)+10*np.log10((rayleigh_amp**2))

    # generate some array for results
    snr_single_channels = np.zeros((config.n_dims,mc_runs))
    snr_comb_channel = np.zeros((len(Nr_list), mc_runs))
    snr_comb_theory = np.zeros((len(Nr_list), mc_runs))

    # doing MC_runs:
    run_count = 1
    for run in range(mc_runs):

        # do CH simulation
        sig_ch = copy.deepcopy(sig_tx)
        sig_ch = ch_sim(sig_ch, config, SNR_arr[:, run])

        # do RX processing single channels
        for dim in range(config.n_dims):
            sig_ch_temp = copy.deepcopy(sig_ch)
            snr_single_channels[dim, run], _ = rx_sim(sig_ch_temp,config,which_dim=dim)

        for Nr in range(len(Nr_list)): 
            # calc expected SNR
            snr_comb_theory[Nr, run] = 10*np.log10(np.sum(10**(snr_single_channels[:Nr_list[Nr], run]/10)))
            
            # do combining
            sig_ch_temp = copy.deepcopy(sig_ch)
            sig_ch_comb_mrc = skc.rx.combining(sig_ch_temp, comb_method="MRC", weights=(10**(snr_single_channels[:, run]/10)), combine_upto=Nr_list[Nr])

            # do RX processing comb channel
            snr_comb_channel[Nr, run], _ = rx_sim(sig_ch_comb_mrc,config,which_dim=0)

        # make some sort of progress bar 
        if run_count%5 == 0: 
            print("{}/{}".format(run_count, mc_runs))        
        run_count += 1

    # save as pickle file
    if save_data == True:
        save_dict = {
            "snr_single_channels": snr_single_channels, 
            "snr_comb_channel" : snr_comb_channel,
            "snr_comb_theory" : snr_comb_theory,
            "SNR_arr" : SNR_arr,
            "Nr_list": Nr_list
        }
        skc.utils.save_pickle(save_dict, folder='.', f_name='sim_data', add_timestamp=False)
        # save config as well
        skc.utils.save_pickle(config, folder='.', f_name='config_obj', add_timestamp=False)

    # plotting stuff
    # plot raw runs
    plt.figure(1)
    for dim in range(config.n_dims):
        plt.plot(snr_single_channels[dim,:], label="dim {}".format(dim), linestyle="", marker="o")
    for Nr in range(len(Nr_list)):
        plt.plot(snr_comb_theory[Nr,:], label="MRC comb theory $N_R$="+str(Nr_list[Nr]), linestyle="--")
        plt.plot(snr_comb_channel[Nr,:], label="MRC comb $N_R$="+str(Nr_list[Nr]), linestyle="", marker="o")
    plt.legend()
    plt.ylabel("Es/N0 [dB]")
    plt.xlabel("raw MC run")

    diff_theo_sim = snr_comb_theory-snr_comb_channel
    plt.figure(2)
    for Nr in range(len(Nr_list)):
        plt.plot(diff_theo_sim[Nr,:], label="Nr of branches: {}".format(Nr_list[Nr]))
    plt.ylim([-0.2,0.2])
    plt.legend()
    plt.ylabel("Diff theory measurement [dB]")
    plt.xlabel("raw MC run")

    # check if rayleigh is rayleigh as expected
    if mc_runs < 9:
        amount_rayleigh_samp = int(1e6)
        rayleigh_amp = np.abs(rng.normal(0,np.sqrt(0.5), size=(config.n_dims,amount_rayleigh_samp))+1j*rng.normal(0,np.sqrt(0.5), size=(config.n_dims,amount_rayleigh_samp))).flatten()
    plt.figure(3)
    x = np.linspace(rayleigh.ppf(0.01), rayleigh.ppf(0.99), 100)
    plt.plot(x,rayleigh.pdf(x, scale=np.sqrt(0.5)), label="rayleigh pdf with $\sigma^2 = 0.5$")
    plt.hist(rayleigh_amp.flatten(), density=True, label="sim rayleigh amplitude")
    mean_ray_power = np.mean(np.abs(rayleigh_amp.flatten()**2))
    plt.plot([mean_ray_power, mean_ray_power],[0,1], linestyle="--", color="k", label="$\overline{|X|^2}$ should be 1 (Nuszkowski)")
    plt.legend()
    plt.grid()
    plt.ylabel("PDF")
    plt.xlabel("Rayleigh-Amp")

    plt.figure(4)
    plt.hist(SNR_arr.flatten(), density=True, label="SNR")
    plt.legend()
    plt.grid()
    plt.ylabel("PDF")
    plt.xlabel("SNR [dB]")

    # # check if there is a huge difference between snr set and snr est
    # diff_snr_est_snr_set = snr_single_channels - EsN0_dB
    # plt.figure(3)
    # for dim in range(config.n_dims):
    #     plt.plot(diff_snr_est_snr_set[dim,:])
    # plt.xlabel("raw MC run")
    # plt.ylabel("diff set/est SNR [dB]")

    plt.show()
    print("done")

           