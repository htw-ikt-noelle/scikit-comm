from tkinter import DOTBOX
import comm as comm
import numpy as np
import matplotlib.pyplot as plt
import sim_routines as routs
from sim_config import sim_param
import copy
import multiprocessing as mp
from multiprocessing import Process, Manager
import time

sim_par = sim_param()
sim_par.amount_of_symbols = 1e3
sim_par.format = "QAM"
sim_par.order = 4
sim_par.DAC_SR = 50e6
sim_par.symbol_rate = sim_par.DAC_SR/1
sim_par.roll_off = 0.1

n_dims_list = [2,3,4,5,6,7,8,9,10]
MC_runs = 10000
snr_mean = 20

##########################################


snr_middle_samples_lin = np.zeros((len(n_dims_list),MC_runs),dtype=float)
snr_middle_samples_dB = np.zeros((len(n_dims_list),MC_runs),dtype=float)
snr_mrc_dB = np.zeros((len(n_dims_list),MC_runs),dtype=float)
snr_egc_dB = np.zeros((len(n_dims_list),MC_runs),dtype=float)

sig_tx = routs.sim_SISO_TX(sim_par, seed_bits=1)

def worker(sig_tx, n_dims, MC_runs, snr_mean, result_dict):
    time_tic = time.time()
    result_dict[str(n_dims)+"_snr_mrc_dB"] = np.zeros(MC_runs, dtype=float)
    result_dict[str(n_dims)+"_snr_egc_dB"] = np.zeros(MC_runs, dtype=float)
    result_dict[str(n_dims)+"_snr_middle_samples_dB"] = np.zeros(MC_runs, dtype=float)
    result_dict[str(n_dims)+"_snr_middle_samples_lin"] = np.zeros(MC_runs, dtype=float)

    temp_list_snr_mrc_dB = np.zeros(MC_runs, dtype=float)
    temp_list_snr_egc_dB = np.zeros(MC_runs, dtype=float)
    temp_list_snr_middle_samples_lin = np.zeros(MC_runs, dtype=float)
    temp_list_snr_middle_samples_dB = np.zeros(MC_runs, dtype=float)

    #snr_list = np.full(n_dims,fill_value=30)
    #snr_list = np.random.default_rng().uniform(10,30,n_dims)
    
    for j in range(MC_runs):
        #setting samples array
        sample_array = np.zeros((n_dims,len(sig_tx.samples[0])),dtype=complex)
        samples_mrc_signal = np.zeros(len(sample_array[0,:]), dtype=complex)
        samples_egc_signal = np.zeros(len(sample_array[0,:]), dtype=complex)

        snr_list = 10*np.log10(np.random.default_rng().rayleigh(scale=(10**(snr_mean/10))/(np.sqrt(np.pi/2)), size=n_dims))

        for i in range(n_dims):
            #create Signal
            sample_array[i,:] = sig_tx.samples[0]

            #setting snr
            sample_array[i,:] = comm.channel.set_snr(sample_array[i,:],snr_dB=snr_list[i],sps=int(sim_par.DAC_SR/sim_par.symbol_rate))

            #combining
            samples_mrc_signal += sample_array[i,:]*np.sqrt(10**(snr_list[i]/10))
            samples_egc_signal += sample_array[i,:]
        
        #normalizing
        mag_const = np.mean(abs(sig_tx.constellation[0])) 
        mag_samples = np.mean(abs(samples_mrc_signal)) 
        samples_mrc_signal = samples_mrc_signal * (mag_const / mag_samples)

        mag_samples = np.mean(abs(samples_egc_signal)) 
        samples_egc_signal = samples_egc_signal * (mag_const / mag_samples)

        #calulate SNR of combined vectors
        S = np.var(sig_tx.samples[0])
        N_mrc = np.var(samples_mrc_signal - sig_tx.samples[0])
        N_egc = np.var(samples_egc_signal - sig_tx.samples[0])
        temp_list_snr_mrc_dB[j] = 10*np.log10(S/N_mrc)
        temp_list_snr_egc_dB[j] = 10*np.log10(S/N_egc)

        #calculate middle SNR of normal samples
        temp_list_snr_middle_samples_lin[j] = np.mean([10**(x/10) for x in snr_list])
        temp_list_snr_middle_samples_dB[j] = 10*np.log10(temp_list_snr_middle_samples_lin[j])
    
    result_dict[str(n_dims)+"_snr_mrc_dB"] = temp_list_snr_mrc_dB
    result_dict[str(n_dims)+"_snr_egc_dB"] = temp_list_snr_egc_dB
    result_dict[str(n_dims)+"_snr_middle_samples_lin"] =  temp_list_snr_middle_samples_lin 
    result_dict[str(n_dims)+"_snr_middle_samples_dB"] =  temp_list_snr_middle_samples_dB

    time_toc = time.time()
    print("Dimension {} done in {:.2f}s.".format(n_dims, (time_toc-time_tic)))

if __name__ == "__main__":
    manager = Manager()
    result_dict = manager.dict()

    jobs = []
    for n_dims in n_dims_list:
        p = mp.Process(target=worker, args=(sig_tx, n_dims, MC_runs, snr_mean, result_dict))
        jobs.append(p)
        p.start()

    for proc in jobs:
        proc.join()

    #print(result_dict)

    plt.figure(1)
    mean_mrc_snr_list = np.zeros(len(n_dims_list))
    mean_egc_snr_list = np.zeros(len(n_dims_list))
    plt.title("Simulation with {:.1E} symbols and {:.1E} MC runs, \n SNR mean overall = {:.2f} dB".format(sim_par.amount_of_symbols, MC_runs, snr_mean))

    for i in range(len(n_dims_list)):
        plt.plot(np.full(MC_runs,fill_value=n_dims_list[i]), result_dict[str(n_dims_list[i])+"_snr_mrc_dB"]-result_dict[str(n_dims_list[i])+"_snr_middle_samples_dB"], linestyle="", marker="o", color="cornflowerblue", alpha=0.01)
        plt.plot(np.full(MC_runs,fill_value=n_dims_list[i]), result_dict[str(n_dims_list[i])+"_snr_egc_dB"]-result_dict[str(n_dims_list[i])+"_snr_middle_samples_dB"], linestyle="", marker="o", color="greenyellow", alpha=0.01)
        mean_mrc_snr_list[i] = np.mean(result_dict[str(n_dims_list[i])+"_snr_mrc_dB"]-result_dict[str(n_dims_list[i])+"_snr_middle_samples_dB"])
        mean_egc_snr_list[i] = np.mean(result_dict[str(n_dims_list[i])+"_snr_egc_dB"]-result_dict[str(n_dims_list[i])+"_snr_middle_samples_dB"])
        
    plt.plot(n_dims_list, mean_mrc_snr_list, linestyle="--", marker="o", color="blue", label="MRC means", lw=2)
    plt.plot(n_dims_list, mean_egc_snr_list, linestyle="--", marker="o", color="lime", label="EGC means", lw=2)
    plt.plot(n_dims_list,10*np.log10(n_dims_list),color='slateblue', linestyle="-.", lw=2, label="MRC Theory")
    plt.plot(n_dims_list,10*np.log10([x*np.pi/4 for x in n_dims_list]),color='seagreen', lw=2, linestyle="-.", label="EGC Theory")
    plt.ylim(0,15)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    diff_mrc = 10*np.log10(n_dims_list[2])-mean_mrc_snr_list[2]
    diff_egc = 10*np.log10([x*np.pi/4 for x in n_dims_list])[2]-mean_egc_snr_list[2]
    textstr = '\n'.join(("Values at order 4:", "Diff MRC/Theory: %.2f dB" % (diff_mrc,), "Diff EGC/Theory: %.2f dB" % (diff_egc,)))
    plt.text(2.05, 14, textstr, fontsize=10, verticalalignment='top', bbox=props)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    diff_mrc = 10*np.log10(n_dims_list[-1])-mean_mrc_snr_list[-1]
    diff_egc = 10*np.log10([x*np.pi/4 for x in n_dims_list])[-1]-mean_egc_snr_list[-1]
    textstr = '\n'.join(("Values at last order:", "Diff MRC/Theory: %.2f dB" % (diff_mrc,), "Diff EGC/Theory: %.2f dB" % (diff_egc,)))
    plt.text(2.05, 11.5, textstr, fontsize=10, verticalalignment='top', bbox=props)

    plt.grid()
    plt.legend()
    plt.xlabel("Number of Apertures")
    plt.ylabel("SNR combining gain [dB]")

    #plt.savefig("overview_snr_{:.2f}_symbols_{:.2f}_MC_runs_{:.2f}.png".format(snr_mean, sim_par.amount_of_symbols, MC_runs))


    plt.figure(2)
    mean_snr_middle_samples_dB = np.zeros(len(n_dims_list))
    plt.title("middle snr of MC_runs per n_dims")
    for i in range(len(n_dims_list)):
        plt.plot(np.full(MC_runs,fill_value=n_dims_list[i]), result_dict[str(n_dims_list[i])+"_snr_middle_samples_dB"], linestyle="", marker="o", color="coral", alpha=0.01)
        mean_snr_middle_samples_dB[i] = np.mean(result_dict[str(n_dims_list[i])+"_snr_middle_samples_dB"])
    plt.plot(n_dims_list, mean_snr_middle_samples_dB, linestyle="--", marker="o", color="tomato", lw=2)
    plt.grid()
    #plt.savefig("SNRMiddle_snr_{:.2f}_symbols_{:.2f}_MC_runs_{:.2f}.png".format(snr_mean, sim_par.amount_of_symbols, MC_runs))
    plt.show()
    