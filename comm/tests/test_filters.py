import unittest
import numpy as np
from numpy.testing import (
        assert_, assert_raises, assert_equal, assert_warns,
        assert_no_warnings, assert_array_equal, assert_array_almost_equal,
        suppress_warnings
        )

import os
import sys
sys.path.insert(0, os.path.abspath('..\..'))
import comm as comm


class TestMovingAverage(unittest.TestCase):
    """
    Test class for moving average filer
    """
    
    def test_dirac_even_freq(self):
        sig_in = np.array([1.0, 0, 0, 0, 0, 0, 0, 0])
        sig_out = comm.filters.moving_average(sig_in, average=4, domain='freq')
        assert_array_almost_equal(sig_out, np.array([0.25, 0.25, 0.0, 0.0, 0.0, 0.0, 0.25, 0.25]), decimal=6)
        
    def test_dirac_odd_freq(self):
        sig_in = np.array([1.0, 0, 0, 0, 0, 0, 0, 0])
        sig_out = comm.filters.moving_average(sig_in, average=5, domain='freq')
        assert_array_almost_equal(sig_out, np.array([0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.2, 0.2]), decimal=6)
        
    def test_dirac_even_time(self):
        sig_in = np.array([1.0, 0, 0, 0, 0, 0, 0, 0])
        sig_out = comm.filters.moving_average(sig_in, average=4, domain='time')
        assert_array_almost_equal(sig_out, np.array([0.25, 0.25, 0.25, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), decimal=6)
        
    def test_dirac_odd_time(self):
        sig_in = np.array([1.0, 0, 0, 0, 0, 0, 0, 0])
        sig_out = comm.filters.moving_average(sig_in, average=5, domain='time')
        assert_array_almost_equal(sig_out, np.array([0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), decimal=6)
        