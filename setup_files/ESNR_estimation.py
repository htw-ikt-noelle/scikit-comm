# ESNR estimation - do not push!

import sys
import os
# if not any(os.path.abspath('..') == p for p in sys.path): 
#     print('adding comm module to path...')
#     sys.path.insert(0, os.path.abspath('..'))
import numpy as np
import matplotlib.pyplot as plt
import comm as comm

# global
TX_UPSAMPLE_FACTOR = 5
SNR = 10 #in dB

# construct signal
sig_tx = comm.signal.Signal(n_dims=1)
sig_tx.symbol_rate = 50e6 

# generate bits
sig_tx.generate_bits(n_bits=2**14, seed=1)

# set constellation (modulation format)
sig_tx.generate_constellation(order=16)
sig_tx.modulation_info = 'QPSK'

# create symbols
sig_tx.mapper()

# sig_tx.samples = sig_tx.symbols[0]
# sig_tx.sample_rate = sig_tx.symbol_rate[0]
# upsampling and pulseshaping
ROLL_OFF = 0.2
sig_tx.pulseshaper(upsampling=TX_UPSAMPLE_FACTOR, pulseshape='rrc', roll_off=[ROLL_OFF])

# add complex noise
sig_tx.set_snr(snr_dB=SNR,seed=None)

# RX matched filter
sig_tx.samples = comm.filters.raised_cosine_filter(sig_tx.samples[0],root_raised=True,roll_off=ROLL_OFF,
                                                   symbol_rate=sig_tx.symbol_rate[0],
                                                   sample_rate=sig_tx.sample_rate[0])


# downsample to 1 sps
sig_tx.samples = sig_tx.samples[0][::TX_UPSAMPLE_FACTOR]

sig_tx.plot_constellation()


# TODO: estimate SNR from (electrical) samples

#### "improved SNR estimation algorithm"
# calc mean
symb_mean = np.mean(np.abs(sig_tx.samples[0]))

# calc variance
symb_var = np.mean((np.abs(sig_tx.samples[0])-symb_mean)**2)


# estimate SNR
snr_est = 10*np.log10((np.abs(symb_mean)**2)/(2*symb_var))

# estimated SNR is higher than actual SNR, the error increases as the actual
# SNR decreases; a modification of the estimation is performed:

if snr_est < 10:
    # introduce alias for sample vector for better legibility
    y = sig_tx.samples[0]
    
    # calc z (Qun & Jian, Eq. 8)
    z = snr_est
    #z = (np.mean(np.real(y)**2)+np.mean(np.imag(y)**2)) / (np.abs(np.mean(np.real(y)))**2+np.abs(np.mean(np.imag(y)))**2)
    
    # calc modified SNR, which is valid instead of snr_est
    #snr_est_mod = 1e4*(-0.041292958452235*(z**5)+2.66418532072905*(z**4)-6.86724072350538*(z**3)+8.84039993634297*(z**2)-5.68658561155135*z+1.464045795143920)
    snr_est_mod = np.sqrt((z-2.5)*39.2)-7

#### M2M4 estimation algorithm

# not entirely clear on how to obtain expected values for m2 & m4
m2 = np.mean(np.abs(sig_tx.samples[0]**2))
m4 = np.mean(np.abs(sig_tx.samples[0]**4))

rho_numerator = 0.5 * np.sqrt(6*(m2**2)-(2*m4))
rho_denominator = m2 - (0.5 * np.sqrt(6*(m2**2)-(2*m4)))
rho_hat = rho_numerator/rho_denominator

#### theta_2 estimator for QPSK 
# ref.: Pauluzzi & Beaulieu, "Comparison of Four SNR Estimators for QPSK Modulations", 2000

X_i = np.real(sig_tx.samples[0])
Y_i = np.imag(sig_tx.samples[0])
theta_2_numerator = (np.abs(X_i) - np.abs(Y_i))**2
theta_2_denominator = (X_i**2) + (Y_i**2)
theta_2_lin = len(sig_tx.samples[0])*((np.sum(theta_2_numerator/theta_2_denominator))**-1)
theta_2_dB = 10*np.log10(theta_2_lin)



#### Zuo, Liu QAM SNR estimator

# lookup table for estimation coefficients
coeff_lut = np.asarray([[16,32,64,128,256], # order
                        [17,113,13,905,4369], # a
                        [8,87,8,776,2856]]) # b

# I and Q components of received signal
y_i = np.real(sig_tx.samples[0])
y_q = np.imag(sig_tx.samples[0])

# calc lambda_MQAM
lmbd_numerator = np.mean((y_i**2+y_q**2)**2) - np.mean(y_i**2 + y_q**2)**2
lmbd_denominator = (np.mean(y_i**2+y_q**2))**2
lmbd_mqam = lmbd_numerator / lmbd_denominator

# pull coefficients from LUT
lut_idx = np.where(coeff_lut[0]==len(sig_tx.constellation[0]))
a_mqam = coeff_lut[1][lut_idx]
b_mqam = coeff_lut[2][lut_idx]

# calc SNR estimate
mqam_snr_estimate_numerator = np.sqrt((1-lmbd_mqam)*a_mqam/(a_mqam+b_mqam)) - (lmbd_mqam-1)
mqam_snr_estimate_denominator = lmbd_mqam - (b_mqam/(a_mqam + b_mqam))
mqam_snr_estimate = 10*np.log10(mqam_snr_estimate_numerator / mqam_snr_estimate_denominator)


#### Xu, Li, Zheng QAM SNR estimator
# I-/Q-component based estimator using a fifth-order polynomial to approximate SNR
# in the range of [-5,20] dB with coefficient lookup table

r_kI = np.real(sig_tx.samples[0])
r_kQ = np.imag(sig_tx.samples[0])

z_hat = (np.mean(r_kI**2)+np.mean(r_kQ**2)) / (np.mean(np.abs(r_kI))**2+np.mean(np.abs(r_kQ))**2)

# lookup table of polynomial coefficients
# TODO: implement coeff lists for other constellations
C_16 = 1e6 * np.asarray([-0.06503489292716,0.45823466671427,-1.29109195371317,1.81839487545758,-1.28034019700542,0.36060357620798])

# dictionary eith supported QAM orders and coefficient vectors in the form [C0,C1,C2,C3,C4,C5]
coeff_dict = {'qam_order' : [16,32,64,128,256],
              'coeff' : [[0.36060357620798,-1.28034019700542,1.81839487545758,-1.29109105371317,0.45823466671427,-0.06503489292716],
                         [0.53715056289170,-1.85210885301961,2.55864235321160,-1.76993734603623,0.61298700208470,-0.08502242157078],
                         [1.81625572448046,-6.24952901412163,8.60050533607873,-5.91706608901663,2.03511551491328,-0.27993710478023],
                         [0.64033054858630,-2.17678215614423,2.95932583860006,-2.01114864174439,0.68323069211818,-0.09282225372024],
                         [0.33595278506244,-1.15419807009244,1.58563212231193,-1.08880229086714,0.37369521988006,-0.05128588224013]]}
Xu_estimator = 10*np.log10(C_16[0]*(z_hat**5)+C_16[1]*(z_hat**4)+C_16[2]*(z_hat**3)+C_16[3]*(z_hat**2)+C_16[4]*(z_hat)+C_16[5])