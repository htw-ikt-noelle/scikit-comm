import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import comm as comm
import copy 
import warnings

## options for test
amount_symbols_wanted = int(1e2) #must be higher than 4 for cropping
timedelay_in_sampels = 10
symbol_rate = 1
upsampling = 2
phase_offset_list = np.linspace(-2*np.pi,2*np.pi,1000)
est_phase = np.zeros(len(phase_offset_list))

for i in range(0,len(phase_offset_list)):
    #### GENERATE SIGNALS
    #generate some signal
    mod_format = 'QAM'
    mod_order = 4
    sig1 = comm.signal.Signal(n_dims=int(1))
    sig1.symbol_rate[0] = symbol_rate
    bits = (amount_symbols_wanted*2)*2
    sig1.generate_bits(n_bits=bits,seed=1)
    sig1.generate_constellation(format=mod_format,order=mod_order)
    sig1.mapper()
    sig1.pulseshaper(upsampling=upsampling,pulseshape='rrc', roll_off=0.2)
    sig2 = copy.deepcopy(sig1)

    #add phase offset to sig 2
    sig2.samples[0] = sig2.samples[0]*np.exp(1j*phase_offset_list[i])

    # roll and cut sampels (cut -> same legnth but some samples different for reality)
    sig2.samples[0] = np.roll(sig2.samples[0],timedelay_in_sampels)

    if timedelay_in_sampels != 0: 
        sig2.samples[0] = sig2.samples[0][int((amount_symbols_wanted/2)*upsampling):-int((amount_symbols_wanted/2)*upsampling)]
        sig1.samples[0] = sig1.samples[0][int((amount_symbols_wanted/2)*upsampling):-int((amount_symbols_wanted/2)*upsampling)]
    else:
        pass

    t1 = np.arange(0, (len(sig1.samples[0])/sig1.sample_rate[0]),1/sig1.sample_rate[0])
    t2 = np.arange(0, (len(sig2.samples[0])/sig1.sample_rate[0]),1/sig1.sample_rate[0])

    # tit = "Signals before comepensating time + phase"
    # plt.figure()
    # plt.title(tit+str(", constellation"))
    # plt.scatter(sig1.samples[0].real, sig1.samples[0].imag, label="Sig 1")
    # plt.scatter(sig2.samples[0].real, sig2.samples[0].imag, label="Sig 2")
    # plt.legend()
    # plt.grid()
    # plt.show()

    # plt.figure()
    # plt.title(tit+str(", timebase"))
    # plt.plot(t1,sig1.samples[0].real, label="sig1 real")
    # plt.plot(t1,sig1.samples[0].imag, label="sig1 imag",ls="--", marker="o")
    # plt.plot(t2,sig2.samples[0].real, label="sig2 real")
    # plt.plot(t2,sig2.samples[0].imag, label="sig2 imag",ls="--", marker="*")
    # plt.ylabel("Amplitude")
    # plt.xlabel("time [s]")
    # plt.legend()
    # plt.show()    

    #### COMPENSATE FOR TIME DELAY
    sig1.samples[0], sig2.samples[0], lag_sr = comm.rx.comb_timedelay_compensation(sig1.samples[0], sig2.samples[0], sr=sig1.sample_rate[0], method="crop", xcorr="abs")
    print("estimated lag: "+str(lag_sr*upsampling*symbol_rate))
    
    t = np.arange(0, (len(sig1.samples[0])/sig1.sample_rate[0]),1/sig1.sample_rate[0])

    # tit = "(cropped) Signals after time compensation"
    # plt.figure()
    # plt.title(tit+str(", constellation"))
    # plt.scatter(sig1.samples[0].real, sig1.samples[0].imag, label="Sig 1")
    # plt.scatter(sig2.samples[0].real, sig2.samples[0].imag, label="Sig 2")
    # plt.legend()
    # plt.grid()
    # plt.show()

    # plt.figure()
    # plt.title(tit+str(", timebase"))
    # plt.plot(t,sig1.samples[0].real, label="sig1 real")
    # plt.plot(t,sig1.samples[0].imag, label="sig1 imag",ls="--", marker="o")
    # plt.plot(t,sig2.samples[0].real, label="sig2 real")
    # plt.plot(t,sig2.samples[0].imag, label="sig2 imag",ls="--", marker="*")
    # plt.ylabel("Amplitude")
    # plt.xlabel("time [s]")
    # plt.legend()
    # plt.show()

    # ### COMPENSATE FOR PHASE
    print("should be the phase: "+str(np.rad2deg(phase_offset_list[i])))
    sig1.samples[0], sig2.samples[0], est_phase[i] = comm.rx.comb_phase_compensation(sig1.samples[0], sig2.samples[0])

    # tit = "Signals after time + phase compensation"
    # plt.figure()
    # plt.title(tit+str(", constellation"))
    # plt.scatter(sig1.samples[0].real, sig1.samples[0].imag, label="Sig 1")
    # plt.scatter(sig2.samples[0].real, sig2.samples[0].imag, label="Sig 2", marker="*")
    # plt.legend()
    # plt.grid()
    # plt.show()

    # plt.figure()
    # plt.title(tit+str(", timebase"))
    # plt.plot(t,sig1.samples[0].real, label="sig1 real")
    # plt.plot(t,sig1.samples[0].imag, label="sig1 imag",ls="--", marker="o")
    # plt.plot(t,sig2.samples[0].real, label="sig2 real")
    # plt.plot(t,sig2.samples[0].imag, label="sig2 imag",ls="--", marker="*")
    # plt.ylabel("Amplitude")
    # plt.xlabel("time [s]")
    # plt.legend()
    # plt.show()

plt.figure()
plt.plot(np.rad2deg(est_phase), label="estimated offset")
plt.ylabel("Phase in [°]")
plt.xlabel("Nummer Messung")
plt.plot(np.rad2deg(phase_offset_list), label="set phase")
plt.grid()
plt.legend()
plt.show()