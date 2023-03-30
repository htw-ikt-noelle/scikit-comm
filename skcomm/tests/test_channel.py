import unittest

import numpy as np
from numpy.testing import (
        assert_, assert_raises, assert_equal, assert_warns,
        assert_no_warnings, assert_array_equal, assert_array_almost_equal,
        suppress_warnings
        )

from .. import channel

class TestSetSNR(unittest.TestCase):
    """
    Test class for set_snr
    """
    
    def test_high_snr(self):
        sig_in = np.array([1, 1, 1, 1])
        sig_out = channel.set_snr(sig_in, snr_dB=100)
        assert_array_almost_equal(sig_in, sig_out, decimal=2)
        
        
class TestAddPhaseNoise(unittest.TestCase):
    """
    Test class for add_phase_noise
    """
    
    def test_low_pn(self):
        sig_in = np.array([1+0.j, 1+0.j, 1+0.j, 1+0.j])
        sig_out = channel.add_phase_noise(sig_in, linewidth=0.0)['samples']
        assert_array_equal(sig_in, sig_out)
        
    def test_no_amp_change(self):
        sig_in = np.array([1+0.j, 1+0.j, 1+0.j, 1+0.j])
        sig_out = channel.add_phase_noise(sig_in, linewidth=1.0)['samples']
        assert_array_almost_equal(np.abs(sig_in), np.abs(sig_out), decimal=10)
        
        


