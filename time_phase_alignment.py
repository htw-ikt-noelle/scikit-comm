import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import comm as comm
import copy 
import warnings
import scipy.signal as ssignal
import scipy.interpolate as sinterp

def gen_SIMO_samples(sig, max_phase_offset_in_rad=np.pi/3, max_timedelay_in_percent=10, n_apertures=2, repeat=5, cut_to=3, seed=None, subsample_shift=True):
    """
    generate some realistic SIMO signals 
    """
    len_orig = len(sig.samples[0])
    sig.samples[0] = np.tile(sig.samples[0],reps=repeat)
    len_vectors = len(sig.samples[0])
    
    rng = np.random.default_rng(seed=seed)
    time_delay = rng.uniform(0,max_timedelay_in_percent,size=n_apertures)
    time_delay_in_samples = np.around(time_delay*(len_orig/100),0)
    time_delay_in_subsamples = time_delay*(len_orig/100)-time_delay_in_samples

    phase_offset = rng.uniform(0,max_phase_offset_in_rad,size=n_apertures)

    time_delay_in_samples[0] = 0 
    time_delay_in_subsamples[0] = 0
    phase_offset[0] = 0
    
    if subsample_shift != True:
        time_delay_in_subsamples = np.zeros((n_apertures))

    samples_to_proceed = np.zeros((n_apertures, int(len_vectors+max(time_delay_in_samples))), dtype=complex)
    for n in range(n_apertures):
        samples_to_proceed[n,0:len_vectors] = sig.samples[0]
        samples_to_proceed[n,:] = samples_to_proceed[n,:]*np.exp(1j*phase_offset[n])
        samples_to_proceed[n,:] = np.roll(samples_to_proceed[n,:], shift=int(time_delay_in_samples[n]))
        samples_to_proceed[n,:] = comm.filters.time_shift(samples_to_proceed[n,:], sample_rate=sig.sample_rate[0],tau=time_delay_in_subsamples[n]/sig.symbol_rate[0])

    #crop signals to cut_to
    return_samples = np.zeros((n_apertures, len_orig*cut_to), dtype=complex)
    for n in range(n_apertures):
        return_samples[n,:] = samples_to_proceed[n,(int(max(time_delay_in_samples))):(int(len_orig*cut_to+max(time_delay_in_samples)))]
        sig.samples[n] = return_samples[n,:] 

    return_dict = {
        "time_delay_in_samples": time_delay_in_samples,
        "time_delay_in_subsamples": time_delay_in_subsamples,
        "phase_offset": np.rad2deg(phase_offset)
    }

    return sig, return_dict, return_samples

def plot_signal_timebase(sig1, tit=""):
    t = np.arange(0, (len(sig1.samples[0])/sig1.sample_rate[0]),1/sig1.sample_rate[0])
    plt.figure()
    plt.title(tit+str(", constellation"))
    plt.scatter(sig1.samples[0].real, sig1.samples[0].imag, label="Samples 1")
    plt.scatter(sig1.samples[1].real, sig1.samples[1].imag, label="Samples 2")
    plt.legend()
    plt.grid()
    plt.show()

    plt.figure()
    plt.title(tit+str(", timebase"))
    plt.plot(t,sig1.samples[0].real, label="samples 1 real")
    plt.plot(t,sig1.samples[0].imag, label="samples 1 imag",ls="--", marker="o")
    plt.plot(t,sig1.samples[1].real, label="samples 2 real")
    plt.plot(t,sig1.samples[1].imag, label="samples 2 imag",ls="--", marker="*")
    plt.xlim(0,20/sig.sample_rate[0])
    plt.ylabel("Amplitude")
    plt.xlabel("time [s]")
    plt.legend()
    plt.show()

## options for test
amount_symbols_wanted = 20000 #must be higher than 4 for cropping
timedelay_in_sampels = 2 # muss gerade noch eine gerade Zahl sein, da ansonsten subsample shift aussteigt
#subsample_relative_delay = 0
subsample_relative_delay = 0.2 #[UI] - Unit Intervall - Nachrichtentechnik ein Symbolintervall
symbol_rate = 1
upsampling = 2
phase_offset = np.pi/3
n_dims = 2
#phase_offset = np.pi

#### GENERATE SIGNALS
#generate some signal
mod_format = 'QAM'
mod_order = 4
sig1 = comm.signal.Signal(n_dims=n_dims)
sig1.symbol_rate[0] = symbol_rate
bits = int(amount_symbols_wanted*np.log2(mod_order))
sig1.generate_bits(n_bits=bits,seed=None)
sig1.generate_constellation(format=mod_format,order=mod_order)
sig1.mapper()
sig1.pulseshaper(upsampling=upsampling,pulseshape='rrc', roll_off=0.2)

sig, return_dict, _ = gen_SIMO_samples(sig1, subsample_shift=True, max_phase_offset_in_rad=np.pi/3)
#print(return_dict)

tit = "Signals before comepensating time + phase"
#plot_signal_timebase(sig1, tit=tit)

results = comm.rx.sampling_phase_adjustment(sig1.samples[1], sample_rate=sig1.sample_rate[0], symbol_rate=sig1.symbol_rate[0], shift_dir='both')
sig1.samples[1] = results['samples_out']
#print(results['est_shift'])

tit = "Signals after comepensating time (subsample)"
#plot_signal_timebase(sig1, tit=tit)

sig1.samples[0], sig1.samples[1], lag = comm.rx.comb_timedelay_compensation(sig1.samples[0], sig1.samples[1], method="crop", xcorr="abs")
#print(lag)

#tit = "Signals after comepensating time (sample)"
#plot_signal_timebase(sig1, tit=tit)

print("time delay in SUBsamples // soll: {} ist: {}".format(return_dict["time_delay_in_subsamples"][1], results['est_shift']))
print("time delay in samples // soll: {} ist: {}".format(return_dict["time_delay_in_samples"][1], lag))

sig1.samples[0], sig1.samples[1], phase_est = comm.rx.comb_phase_compensation(sig1.samples[0], sig1.samples[1])

tit = "Signals after phase comepensating"
plot_signal_timebase(sig1, tit=tit)


print("phase // soll: {} ist: {}".format(return_dict["phase_offset"][1], np.rad2deg(phase_est)))