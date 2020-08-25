import numpy as np
import comm as comm
import matplotlib.pyplot as plt
import scipy.signal as ssignal
import time

###################### Tx ##################################
# contruct signal
sig_tx = comm.signal.Signal(n_dims=1)
sig_tx.symbol_rate = 50e6 #50e6

# generate bits
sig_tx.generate_bits(n_bits=2**10)

# set constellation (modualtion format)
sig_tx.generate_constellation(order=4)
sig_tx.modulation_info = 'QPSK'

# create symbols
sig_tx.mapper()

# upsampling (to 5*50e6 --> 250 MSa/s)  and pulseshaping
roll_off = 0.1
sig_tx.pulseshaper(upsampling=5, pulseshape='rc', roll_off=[roll_off])

# plot for checks
# sig_tx.plot_constellation(0)
# sig_tx.plot_eye(0)
comm.visualizer.plot_signal(sig_tx.samples[0], sample_rate=sig_tx.sample_rate[0])

# generate DAC samples (analytical signalg at IF)
f_if_nom = 30e6 #30e6
f_granularity = 1 / sig_tx.samples[0].size * sig_tx.sample_rate[0]
f_if = round(f_if_nom / f_granularity) * f_granularity
print('intermediate frequency: {} MHz'.format(f_if/1e6))
t = np.arange(0, np.size(sig_tx.samples[0])) / sig_tx.sample_rate
# t = comm.utils.create_time_axis(sig_tx.sample_rate, np.size(sig_tx.samples[0]))

# sig_tx.plot_spectrum(0)

sig_tx.samples[0] = sig_tx.samples[0] * np.exp(1j * 2 * np.pi * f_if * t)
sig_tx.center_frequency = f_if

# sig_tx.plot_spectrum(0)

# format samples so that driver can handle them (range: +-32768, integer, list)
# maybe scale I and Q equally???
maxVal = np.max(np.abs(np.concatenate((np.real(sig_tx.samples), np.imag(sig_tx.samples)))))
samples = np.asarray(sig_tx.samples) / maxVal
samples = np.concatenate((np.real(samples),np.imag(samples)))

# samples = (np.random.rand(2,int(2560))-0.5)
# print(np.max(samples))
# print(np.min(samples))

# comm.instrument_control.write_samples_AWG33522A(samples, ip_adress='192.168.1.45', sample_rate=[250e6,250e6], 
#                             offset=[0.0, 0.0], amp_pp=[4.0,4.0], channels=[1,2], out_filter=['normal']*2)
# time.sleep(0.3)


######################## Rx #################################################
# get samples from scope
sr, samples = comm.instrument_control.get_samples_DLM2034(channels=[1, 2], address='192.168.1.13')

# differential detection
samples = samples[0] - samples[1]
# comm.visualizer.plot_spectrum(samples, sample_rate=sr)
# comm.visualizer.plot_signal(samples, sample_rate=sr)

# downsampling to up times symbolrate
sr_dsp = sig_tx.symbol_rate[0] * 10
# watch out, that this is really an integer, otherwise the samplerate is asynchronous with the data afterwards!!!
len_dsp = sr_dsp / sr * np.size(samples)
if len_dsp % 1:
    raise ValueError('DSP samplerate results in asynchronous sampling of the data symbols')
samples = ssignal.resample(samples, num=int(len_dsp), window='hamming')
sr = sr_dsp
# comm.visualizer.plot_spectrum(samples, sample_rate=sr)

# IQ-Downmixing and (ideal) lowpass filtering
t = np.arange(0, np.size(samples)) / sr
# t = comm.utils.create_time_axis(sr, np.size(samples))
samples_r = samples *  np.cos(2 * np.pi * f_if * t)
# comm.visualizer.plot_spectrum(samples_r, sample_rate=sr)
fc = sig_tx.symbol_rate[0] / 2 * (1 + roll_off) * 1.1 # cuttoff frequency of filter
fc = fc/(sr/2) # normalization to the sampling frequency
samples_r = comm.filters.ideal_lp(samples_r, fc)
# comm.visualizer.plot_spectrum(samples_r, sample_rate=sr)
samples_i = samples *  np.sin(2 * np.pi * f_if * t)
# comm.visualizer.plot_spectrum(samples_i, sample_rate=sr)
samples_i = comm.filters.ideal_lp(samples_i, fc)
# comm.visualizer.plot_spectrum(samples_i, sample_rate=sr)

# contruct rx signal
sig_rx = comm.signal.Signal(n_dims=1)
sig_rx.symbol_rate = sig_tx.symbol_rate
sig_rx.sample_rate = sr
sig_rx.samples[0] = samples_r - 1j * samples_i
sig_rx.plot_spectrum()
# sig_rx.plot_constellation()

# "standard" coherent complex baseband signal processing

# sampling
start = 7
sps = sig_rx.sample_rate[0] / sig_rx.symbol_rate[0] # CHECK FOR INTEGER SPS!!!
tmp = sig_rx.samples[0][start::int(sps)]
comm.visualizer.plot_constellation(tmp)
# sig_rx.samples[0] = sig_rx.samples[0][start::int(sps)]
# sig_rx.plot_constellation()



### Phase estimation (experimental)
rx_symbols = tmp

## implementing viterbi-viterbi for phase estimation
N_CPE = 15 # must be an odd number for Wiener filter
 # number of CPE_taps
mth_Power = 4 # 4=QAM&QPSK, 2=BPSK,...
## rectangular moving average filter
phi_est = np.roll(comm.filters.moving_average(np.unwrap(mth_Power*np.angle(rx_symbols))/mth_Power, N_CPE, domain='freq'), -N_CPE//2)

## Wiener filter (see Ip2007, Acuna2013)
r = 0.3 # r = sigma²_phi / sigma²_ampl;  r>0; r is the ratio between the magnitude of the phase noise variance sigma²_phi and the additive noise variance sigma²
a = 1+r/2-np.sqrt((1+r/2)**2-1) # alpha
h_Wiener = a*r/(1-a**2) * a**np.arange(N_CPE//2+1) # postive half
h_Wiener = np.concatenate((np.flip(h_Wiener[1:]), h_Wiener)) # make symmetric
h_Wiener = h_Wiener / np.sum(h_Wiener) # normalize to unit sum (make unbiased estimator)
plt.figure(2); plt.stem(np.arange(2*(N_CPE//2)+1)-N_CPE//2,h_Wiener,basefmt='C2-'); plt.show()
phi_est = np.roll(comm.filters.filter_samples(np.unwrap(mth_Power*np.angle(rx_symbols))/mth_Power, h_Wiener, domain='freq'), -N_CPE//2+1)
       
rx_symbols = rx_symbols * np.exp(-1j*(phi_est + np.pi/4))
rx_symbols = rx_symbols[1*N_CPE+1:-N_CPE*1] # crop start and end
comm.visualizer.plot_constellation(rx_symbols)





############### TODO ##################################
# # # rx filter (matched filtering)
# # sig.raised_cosine_filter(roll_off=0.2, root_raised=True)
# # sig.plot_eye(0)

# # # # TODO!!!!
# # # # sampling
# # # for i, samples in enumerate(sig.samples):
# # #     sig.samples[i] = samples[15::16]
# # #     sig.sample_rate[i] = sig.symbol_rate[i]


# # sig.plot_eye(0)
# # sig.plot_spectrum(dimension=0)



