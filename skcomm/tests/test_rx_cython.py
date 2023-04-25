import unittest
import copy

import numpy as np
from numpy.testing import assert_array_almost_equal

from .. import rx
from .. import signal


class BlindAdaptiveEqualizerCython(unittest.TestCase):
    """
    Test class for blind adaptive equalizer cython implementation
    """
    
    def test_equal_output(self):
        # generate random input samples
        sig_in = signal.Signal()
        sig_in.symbol_rate = 1
        sig_in.sample_rate = 2
        sig_in.samples = np.random.randn(1000) + 1j * np.random.randn(1000) 
        sig_in.constellation = np.asarray([np.mean(sig_in.samples[0])])
        # run python implementation
        res_py = rx.blind_adaptive_equalizer(copy.deepcopy(sig_in), n_taps=11, mu_cma=1e-3, 
                                                 mu_rde=1e-3, mu_dde=1e-3, decimate=True, 
                                                 return_info=False, stop_adapting=-1, 
                                                 start_rde=200*1, start_dde=200*1,
                                                 compiled=False)
        sig_out_py = res_py['sig'].samples[0]
        # run cython implementation
        res_cy = rx.blind_adaptive_equalizer(copy.deepcopy(sig_in), n_taps=11, mu_cma=1e-3, 
                                                 mu_rde=1e-3, mu_dde=1e-3, decimate=True, 
                                                 return_info=False, stop_adapting=-1, 
                                                 start_rde=200*1, start_dde=200*1,
                                                 compiled=True)
        sig_out_cy = res_cy['sig'].samples[0]
        # test for equality
        assert_array_almost_equal(sig_out_py, sig_out_cy, decimal=12)
        
