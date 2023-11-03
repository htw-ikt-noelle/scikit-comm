# Scikit-comm: Simulation toolbox for communication systems

This repository contains a collection of DSP routines and algorithms to perform numerical simulations on simple communication systems. Its functionality is divided into transmitter, link, and receiver subsystems and - at its current state - contains very limited functionalities.

This project was initially started by [Lutz Molle](https://www.htw-berlin.de/hochschule/personen/person/?eid=12017) and [Markus Nölle](https://www.htw-berlin.de/hochschule/personen/person/?eid=9586) at the [University of Applied Sciences (HTW), Berlin](https://www.htw-berlin.de/).



## 1. The 'signal' object
The 'signal' object can be seen as the 'heart of the toolbox'. It contains all information to describe a modulated data signal. The object can consist of multiple 'dimensions', while each dimension represents in general a two dimensional (or complex) data signal. The structure looks as follows:

```python
def Signal:    
  self.n_dims
  self.samples
  self.center_frequency
  self.sample_rate
  self.bits
  self.symbols
  self.symbol_rate
  self.modulation_info
  self.constellation
```

Many toolbox functions and methods take this signal object as input or output variables. Others in contrast take only a subset of the signal attributes (e.g. the sampled signal (sig.samples)) as input or output.

## 2. Module / Package structure

Besides the signal object, there are multiple other modules avaialable, which provide different functionalities for the simulation of a communication system:

| module                | Description | examples |
| :---                  | :---        |:---      |
|channel.py                |basic function to emulate a transmission channel| set_snr(), add_phase_noise() |
|filters.py | method to filter a discrete signal | raised_cosine_filter(), moving_average()
|instrument_control.py | methods to communicate with laboratory equipment | get_samples_DLM2034(), write_samples_AWG70002B()
|pre_distortion.py |methods in order to perform identification and pre-distortion of devices / systems| generate_wn_probesignal(), estimate_tf_welch
|rx.py| receiver (dsp) subfunctions| demapper(), sampling_phase_adjustment(), carrier_phase_estimation_VV()
|tx.py| transmitter (dsp) subfuntions | generate_bits(), pulseshaper()
|utils.py| utility functions (mostly used by other methods) | dec_to_bits(), create_time_axis()
|visualizers.py| methods to visualize the data signal | plot_spectrum(), plot_constellation(), plot_eye()




## 3. Code Examples 

### 3.1 Minimal working example

This example generates a 50 MBd QPSK modulated signal and demodulates it afterwards:

```python
import copy
import skcomm as skc

##########################
####### TRANSMITTER ######
##########################

# construct signal
sig_tx = skc.signal.Signal(n_dims=1)
sig_tx.symbol_rate = 50e6 

# generate bits
sig_tx.generate_bits(n_bits=2**12, seed=1)

# set constellation (modulation format)
sig_tx.generate_constellation(format='QAM', order=4)

# create symbols
sig_tx.mapper()

# generate actual sampled signal (pulseshaping)
sig_tx.pulseshaper(upsampling=1, pulseshape='rect')

##########################
####### CHANNEL ##########
##########################

pass

##########################
####### RECEIVER #########
##########################

# rx signal
sig_rx = copy.deepcopy(sig_tx)

# decision
sig_rx.decision()

# demapper
sig_rx.demapper()

# BER counting
ber_res = skc.rx.count_errors(sig_rx.bits[0], sig_rx.samples[0])
```

Multiple, more advanced algorithms and procedures could now be added at transmitter (**from module tx**) and receiver side (**from module rx**). Further, also distortion effects caused by the channel (***module channel***) could be added.

### 3.2 More advanced example using QPSK modulation onto an intermediate frequency

This example demonstrates a 50 MBd QPSK modulation using an intermediate frequency of 30 MHz. Further, pulseshaping, matched filtering and also channel distortions, such as amplitude and phase noise as well as simple digital signal processing (DSP) steps at the receiver are considered. 


```python
import numpy as np
import scipy.signal as ssignal
import skcomm as skc
import copy

############################################################
###################### Tx ##################################
############################################################

# signal parameters
LASER_LINEWIDTH = 1*600 # [Hz]
TX_UPSAMPLE_FACTOR = 5
SNR = 20

# construct signal
sig_tx = skc.signal.Signal(n_dims=1)
sig_tx.symbol_rate = 50e6 

# generate bits
sig_tx.generate_bits(n_bits=2**12, seed=1)

# set constellation (modulation format)
sig_tx.generate_constellation(format='QAM', order=4)

# create symbols
sig_tx.mapper()

# upsampling and pulseshaping
ROLL_OFF = 0.1
sig_tx.pulseshaper(upsampling=TX_UPSAMPLE_FACTOR, pulseshape='rrc', roll_off=[ROLL_OFF])

# generate DAC samples (analytical signal at IF)
f_IF_nom = 1*30e6
f_granularity = 1 / sig_tx.samples[0].size * sig_tx.sample_rate[0]
f_if = round(f_IF_nom / f_granularity) * f_granularity
print('intermediate frequency: {} MHz'.format(f_if/1e6))
t = np.arange(0, np.size(sig_tx.samples[0])) / sig_tx.sample_rate

# upmixing to IF
sig_tx.samples[0] = sig_tx.samples[0] * np.exp(1j * 2 * np.pi * f_if * t)
sig_tx.center_frequency = f_if

############################################################################
################## Link ####################################################
############################################################################
samples = sig_tx.samples[0] 

# repeat rx sequence
sps = int(sig_tx.sample_rate[0] / sig_tx.symbol_rate[0])
ext = 40000*sps + 4000*sps
ratio_base = ext // samples.size
ratio_rem = ext % samples.size        
samples = np.concatenate((np.tile(samples, ratio_base), samples[:ratio_rem]), axis=0)

# add artificial delay 
delay = 10*sps
samples = samples[delay:]

## add amplitude noise
samples = skc.channel.set_snr(samples, snr_dB=SNR, sps=int(sig_tx.sample_rate[0]/sig_tx.symbol_rate[0]), seed=None)

## phase noise emulation
samples = skc.channel.add_phase_noise(samples ,sig_tx.sample_rate[0] , LASER_LINEWIDTH, seed=1)['samples']
sr = sig_tx.sample_rate[0]

# after heterodyne detection and balanced detection
samples = np.real(samples)

#############################################################################
######################## Rx #################################################
#############################################################################
# construct rx signal structure
sig_rx = copy.deepcopy(sig_tx)
sig_rx.samples = samples
sig_rx.sample_rate = sr

# resampling to the same sample rate as at the transmitter
sr_dsp = sig_tx.sample_rate[0]

# watch out, that this is really an integer, otherwise the samplerate is asynchronous with the data afterwards!!!
len_dsp = sr_dsp / sig_rx.sample_rate[0] * np.size(samples)
if len_dsp % 1:
    raise ValueError('DSP samplerate results in asynchronous sampling of the data symbols')
sig_rx.samples = ssignal.resample(sig_rx.samples[0], num=int(len_dsp), window=None)
sig_rx.sample_rate = sr_dsp
sig_rx.plot_spectrum(tit='received spectrum before IF downmixing')

# IQ-Downmixing 
t = skc.utils.create_time_axis(sig_rx.sample_rate[0], np.size(sig_rx.samples[0]))
samples_bb = samples *  np.exp(-1j*2*np.pi*(f_if+1e4*0)*t)
sig_rx.samples[0] = samples_bb

# resample to 2 sps
sps_new = 2
sps = sig_rx.sample_rate[0]/sig_rx.symbol_rate[0]
new_length = int(sig_rx.samples[0].size/sps*sps_new)
sig_rx.samples = ssignal.resample(sig_rx.samples[0], new_length, window='boxcar')
sig_rx.sample_rate = sps_new*sig_rx.symbol_rate[0]

# normalize samples to mean magnitude of original constellation
mag_const = np.mean(abs(sig_rx.constellation[0]))
mag_samples = np.mean(abs(sig_rx.samples[0]))
sig_rx.samples = sig_rx.samples[0] * mag_const / mag_samples

sig_rx.plot_constellation(hist=True, tit='constellation before EQ')

# Rx matched filter
sig_rx.raised_cosine_filter(roll_off=ROLL_OFF,root_raised=True) 

# crop samples here, if necessary
sps = int(sig_rx.sample_rate[0] / sig_rx.symbol_rate[0])
crop = 10*sps
if crop != 0:
    sig_rx.samples = sig_rx.samples[0][crop:-crop]
else:
    sig_rx.samples = sig_rx.samples[0]

# sampling phase adjustment
BLOCK_SIZE = -1 
sig_rx.sampling_clock_adjustment(BLOCK_SIZE)
    
# sampling (if necessary)
START_SAMPLE = 0
sps = sig_rx.sample_rate[0] / sig_rx.symbol_rate[0]
sig_rx.samples = sig_rx.samples[0][START_SAMPLE::int(sps)]

sig_rx.plot_constellation(0, hist=True, tit='constellation after EQ')

# CPE
cpe_results = skc.rx.carrier_phase_estimation_bps(sig_rx.samples[0], sig_rx.constellation[0], 
                                           n_taps=15, n_test_phases=15, const_symmetry=np.pi/2)

sig_rx.samples = cpe_results['samples_corrected']
est_phase = cpe_results['est_phase_noise']
    
skc.visualizer.plot_signal(est_phase, tit='estimated phase noise')
sig_rx.plot_constellation(hist=True, tit='constellation after CPE')

# delay and phase ambiguity estimation and compensation
sig_rx = skc.rx.symbol_sequence_sync(sig_rx, dimension=-1)["sig"]
    
# calc EVM
evm = skc.utils.calc_evm(sig_rx, norm='max')
print("EVM: {:2.2%}".format(evm[0]))

# decision and demapper
sig_rx.decision()
sig_rx.demapper()

# BER counting
ber_res = skc.rx.count_errors(sig_rx.bits[0], sig_rx.samples[0])
print('BER = {}'.format(ber_res['ber']))
```