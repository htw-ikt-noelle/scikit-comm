from re import L
import comm as comm
import numpy as np
import matplotlib.pyplot as plt
import copy

#Annahme: beide signale sind im Mittel gleich stark im SNR
#Fragestellung: Bei wievielen Symbolen ist wieviel Lag möglich? 

amount_of_symbols_list = np.linspace(100,1e6,10)
lag_vec = np.zeros(len(amount_of_symbols_list))
symbol_rate = 1
upsampling = 2
snr_dB_list = np.linspace(0,20,10)
max_detection_point_snr_abs = np.zeros(len(snr_dB_list))
max_detection_point_snr_nonabs = np.zeros(len(snr_dB_list))

roll_vec = 1/np.linspace(100,2,40)

for k in range(len(snr_dB_list)):

    for j in range(len(amount_of_symbols_list)):
        mod_format = 'QAM'
        mod_order = 4
        sig0 = comm.signal.Signal(n_dims=int(1))
        sig0.symbol_rate[0] = symbol_rate
        bits = int(amount_of_symbols_list[j])*2
        sig0.generate_bits(n_bits=bits,seed=1)
        sig0.generate_constellation(format=mod_format,order=mod_order)
        sig0.mapper()
        sig0.pulseshaper(upsampling=upsampling,pulseshape='rrc', roll_off=0.2)

        sig1 = copy.deepcopy(sig0)
        sig2 = copy.deepcopy(sig0)

        sig1.samples[0] = comm.channel.set_snr(sig1.samples[0], snr_dB=snr_dB_list[k], sps=upsampling)
        sig2.samples[0] = comm.channel.set_snr(sig2.samples[0], snr_dB=snr_dB_list[k], sps=upsampling)

        stemp1 = sig1.samples[0]
        stemp2 = sig2.samples[0]

        for i in range(len(roll_vec)):
            samples_to_shift = int(roll_vec[i]*(amount_of_symbols_list[j]*upsampling))
            zero_pad = np.zeros(samples_to_shift)
            stemp1 = np.concatenate([sig1.samples[0], zero_pad])
            stemp2 = np.concatenate([sig2.samples[0], zero_pad])

            stemp2 = np.roll(stemp2,(samples_to_shift))

            stemp1 = stemp1[samples_to_shift:-samples_to_shift]
            stemp2 = stemp2[samples_to_shift:-samples_to_shift]

            stemp1, stemp2, lag = comm.rx.comb_timedelay_compensation(stemp1, stemp2, method="crop", xcorr="abs")
            
            if samples_to_shift != abs((lag)):
                break
        
        print(roll_vec[i]*100)
        lag_vec[j] = roll_vec[i]

    # plt.figure()
    # plt.plot(lag_vec, amount_of_symbols_list, ls="", marker="o", markersize=5)
    # plt.xlim([roll_vec[0], roll_vec[-1]])
    # plt.show()

    max_detection_point_snr_abs[k] = np.mean(lag_vec*100)

for k in range(len(snr_dB_list)):

    for j in range(len(amount_of_symbols_list)):
        mod_format = 'QAM'
        mod_order = 4
        sig0 = comm.signal.Signal(n_dims=int(1))
        sig0.symbol_rate[0] = symbol_rate
        bits = int(amount_of_symbols_list[j])*2
        sig0.generate_bits(n_bits=bits,seed=1)
        sig0.generate_constellation(format=mod_format,order=mod_order)
        sig0.mapper()
        sig0.pulseshaper(upsampling=upsampling,pulseshape='rrc', roll_off=0.2)

        sig1 = copy.deepcopy(sig0)
        sig2 = copy.deepcopy(sig0)

        sig1.samples[0] = comm.channel.set_snr(sig1.samples[0], snr_dB=snr_dB_list[k], sps=upsampling)
        sig2.samples[0] = comm.channel.set_snr(sig2.samples[0], snr_dB=snr_dB_list[k], sps=upsampling)

        stemp1 = sig1.samples[0]
        stemp2 = sig2.samples[0]

        for i in range(len(roll_vec)):
            samples_to_shift = int(roll_vec[i]*(amount_of_symbols_list[j]*upsampling))
            zero_pad = np.zeros(samples_to_shift)
            stemp1 = np.concatenate([sig1.samples[0], zero_pad])
            stemp2 = np.concatenate([sig2.samples[0], zero_pad])

            stemp2 = np.roll(stemp2,(samples_to_shift))

            stemp1 = stemp1[samples_to_shift:-samples_to_shift]
            stemp2 = stemp2[samples_to_shift:-samples_to_shift]

            stemp1, stemp2, lag = comm.rx.comb_timedelay_compensation(stemp1, stemp2, method="crop", xcorr="nonabs")
            
            if samples_to_shift != abs((lag)):
                break
        
        print(roll_vec[i]*100)
        lag_vec[j] = roll_vec[i]

    # plt.figure()
    # plt.plot(lag_vec, amount_of_symbols_list, ls="", marker="o", markersize=5)
    # plt.xlim([roll_vec[0], roll_vec[-1]])
    # plt.show()

    max_detection_point_snr_nonabs[k] = np.mean(lag_vec*100)

plt.figure()
plt.title("maximum lag estimation with given methods \n 4-QAM, sps=2, rrc with roll_off=2 \n one point is given by mean of different amount of symbols")
plt.plot(max_detection_point_snr_abs,snr_dB_list, label="xcorr(abs(y...)", marker="*")
plt.plot(max_detection_point_snr_nonabs,snr_dB_list, label="xcorr(y...)", marker="o")
plt.xlabel("Timedelay in [%] of amount of symbols // upper bound of correct detection")
plt.ylabel("SNR [dB]")
plt.legend()
plt.savefig("detection_treshold_timedelay.png")
plt.show()