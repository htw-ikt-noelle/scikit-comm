
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc, gamma
from scipy.integrate import quad

class config_sim():
    def __init__(self) -> None:
        self.laser_linewidth = 100e3
        self.DAC_sr = 16e9
        self.symbol_rate = 12.8e9
        self.n_symbols = 10000 
        self.roll_off = 0.1
        self.n_dims = 1
        self.TX_seed_bits = None
        self.modulation_format = "QAM"
        self.modulation_order = 4
        self.pulseshape = "rrc"
        self.snr_dB = [20]
        self.repeat_samples = 44000
        self.crop_factor_matched_filter = 10
        self.noise_distribution = "AWGN"

def berTheory_rayleighFading_mrcComb(Nr, EbN0dB):
    """
    convert Eb/NO to BER for rayleigh fading enviroment.

    Caution! There is some mismatch between matlab and python if Nr = 1 see powerpoint
    about simulation for documentation.
    """

    # prealighnment of array
    mrcBER = np.zeros(len(EbN0dB))

    # curve generation
    Q = lambda arg: 0.5*erfc(arg/np.sqrt(2))
    for snr_idx, snr_val in enumerate(10**(EbN0dB/10)):
        integrand = lambda arg: Q(np.sqrt(2*snr_val*arg)) * arg**(Nr-1) * np.exp(-arg)
        mrcBER[snr_idx] = 1/gamma(Nr) * quad(integrand,0,np.inf)[0]

    return mrcBER
