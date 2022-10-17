import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import comm as comm
import copy 
import warnings

# ## options for test
# amount_symbols_wanted = int(1e5) #must be higher than 4 for cropping
# timedelay_in_sampels = 100
# symbol_rate = 1
# upsampling = 2
# phase_offset = 0
# snr1_vec = np.concatenate([np.full(10,15), np.linspace(15,30,10), np.full(10,30), np.linspace(30,15,10), np.full(10,15), np.linspace(15,-5,10),np.full(10,-5), np.linspace(-5,15,10), np.full(10,15)])
# roll = 40
# snr2_vec = np.roll(snr1_vec,roll)
# est_phase = np.zeros(len(snr1_vec))
# est_lag_diff_list = np.zeros(len(snr1_vec))
# est_lag_list = np.zeros(len(snr1_vec))

# #### GENERATE SIGNALS
# #generate some signal
# mod_format = 'QAM'
# mod_order = 4
# sig0 = comm.signal.Signal(n_dims=int(1))
# sig0.symbol_rate[0] = symbol_rate
# bits = (amount_symbols_wanted*2)*2
# sig0.generate_bits(n_bits=bits,seed=1)
# sig0.generate_constellation(format=mod_format,order=mod_order)
# sig0.mapper()
# sig0.pulseshaper(upsampling=upsampling,pulseshape='rrc', roll_off=0.2)

# for i in range(0,len(snr1_vec)):
#     #add phase offset to sig 2
#     #sig2.samples[0] = sig2.samples[0]*np.exp(1j*phase_offset_list[i])
#     sig1 = copy.deepcopy(sig0)
#     sig2 = copy.deepcopy(sig0)

#     sig1.samples[0] = comm.channel.set_snr(sig1.samples[0],snr_dB=snr1_vec[i],sps=upsampling)
#     sig2.samples[0] = comm.channel.set_snr(sig2.samples[0],snr_dB=snr2_vec[i],sps=upsampling)

#     # roll and cut sampels (cut -> same legnth but some samples different for reality)
#     sig2.samples[0] = np.roll(sig2.samples[0],timedelay_in_sampels)

#     if timedelay_in_sampels != 0: 
#         sig2.samples[0] = sig2.samples[0][int((amount_symbols_wanted/2)*upsampling):-int((amount_symbols_wanted/2)*upsampling)]
#         sig1.samples[0] = sig1.samples[0][int((amount_symbols_wanted/2)*upsampling):-int((amount_symbols_wanted/2)*upsampling)]
#     else:
#         pass

#     #### COMPENSATE FOR TIME DELAY
#     sig1.samples[0], sig2.samples[0], lag_sr = comm.rx.comb_timedelay_compensation(sig1.samples[0], sig2.samples[0], sr=sig1.sample_rate[0], method="crop", xcorr="abs")
#     est_lag_diff_list[i] = timedelay_in_sampels-abs(lag_sr*upsampling*symbol_rate)
#     est_lag_list[i] = lag_sr*upsampling*symbol_rate

#     # ### COMPENSATE FOR PHASE
#     #print("should be the phase: "+str(np.rad2deg(phase_offset_list[i])))
#     #sig1.samples[0], sig2.samples[0], est_phase[i] = comm.rx.comb_phase_compensation(sig1.samples[0], sig2.samples[0])

# fig, ax1 = plt.subplots()

# ax2 = ax1.twinx()
# ax1.plot(snr1_vec, ls="--", label="SNR sig1", color="seagreen")
# ax1.plot(snr2_vec, ls="--", label="SNR sig2", color="springgreen")

# ax2.plot(est_lag_diff_list, label="Diff in lag estimation", color="brown")

# ax1.legend()
# ax1.set_xlabel('Number of simulation')
# ax1.set_ylabel('set SNR [dB]', color='mediumseagreen')
# ax2.set_ylabel('set_lag - abs(est_lag) in samples', color='brown')
# ax2.set_ylim(-2,30)

# ax2.legend()
# plt.show()

#Annahme: beide signale sind im Mittel gleich stark im SNR
amount_symbols_wanted = int(1e5) #must be higher than 4 for cropping
symbol_rate = 1
upsampling = 2
phase_offset = 0

amount_symbols_wanted = np.linspace(10,1e5,100)
#timedelay_in_samples = np.flipud([1/100,1/10,3/10,1/2])
timedelay_in_samples = np.linspace(1e-3,0.5,50)
est_lagvssymb = np.zeros(len(amount_symbols_wanted))
snr_fix = 10

for i in range(0, len(amount_symbols_wanted)):

    #### GENERATE SIGNALS
    #generate some signal
    mod_format = 'QAM'
    mod_order = 4
    sig0 = comm.signal.Signal(n_dims=int(1))
    sig0.symbol_rate[0] = symbol_rate
    bits = int(amount_symbols_wanted[i])*2
    sig0.generate_bits(n_bits=bits,seed=None)
    sig0.generate_constellation(format=mod_format,order=mod_order)
    sig0.mapper()
    sig0.pulseshaper(upsampling=upsampling,pulseshape='rrc', roll_off=0.2)

    sig1 = copy.deepcopy(sig0)
    sig2 = copy.deepcopy(sig0)

    sig1.samples[0] = comm.channel.set_snr(sig1.samples[0],snr_dB=snr_fix,sps=upsampling)
    sig2.samples[0] = comm.channel.set_snr(sig2.samples[0],snr_dB=snr_fix,sps=upsampling)

    est_lag_diff_list = np.zeros(len(timedelay_in_samples))
    for j in range(0,len(timedelay_in_samples)):
    #add phase offset to sig 2
    #sig2.samples[0] = sig2.samples[0]*np.exp(1j*phase_offset_list[i])
    # roll and cut sampels (cut -> same legnth but some samples different for reality)
        samples_to_shift = int(timedelay_in_samples[j]*len(sig2.samples[0]))
        if samples_to_shift != 0:
            zero_pad = np.zeros(samples_to_shift)
            samples1_temp = np.concatenate([sig1.samples[0], zero_pad])
            samples2_temp = np.concatenate([sig2.samples[0], zero_pad])

            samples2_temp = np.roll(samples2_temp,(samples_to_shift))

            sig1.samples[0] = samples1_temp[samples_to_shift:-samples_to_shift]
            sig2.samples[0] = samples2_temp[samples_to_shift:-samples_to_shift]
        else:
            pass
   
        #### COMPENSATE FOR TIME DELAY
        sig1.samples[0], sig2.samples[0], lag_sr = comm.rx.comb_timedelay_compensation(sig1.samples[0], sig2.samples[0], sr=sig1.sample_rate[0], method="crop", xcorr="abs")

        est_lag_diff_list[j] = samples_to_shift-abs(lag_sr*upsampling*symbol_rate)

    index_first_zero = np.argmax(est_lag_diff_list==0)
    
    est_lagvssymb[i] = timedelay_in_samples[index_first_zero]

    # ### COMPENSATE FOR PHASE
    #print("should be the phase: "+str(np.rad2deg(phase_offset_list[i])))
    #sig1.samples[0], sig2.samples[0], est_phase[i] = comm.rx.comb_phase_compensation(sig1.samples[0], sig2.samples[0])

plt.figure()
plt.plot(amount_symbols_wanted, est_lagvssymb)
plt.show()


# t = np.arange(0,1,1/1000)

# f = [1,2,3]

# A1 = np.sin(2*np.pi*f[0]*t)
# A2 = np.sin(2*np.pi*f[1]*t)
# A3 = np.sin(2*np.pi*f[2]*t)

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# ax.plot(t, A1, f[0])
# ax.plot(t, A2, f[1])
# ax.plot(t, A3, f[2])

# elev, azim, roll = 30, 0, 90
# ax.view_init(azim=-140, elev=30, vertical_axis="y")
# ax.invert_xaxis()

# plt.show()
