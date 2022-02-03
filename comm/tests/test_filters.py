import comm as comm
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


class TestFilterSamples(unittest.TestCase):
    """
    Test class for filter_samples function
    """

    def test_time_domain(self):
        samples = np.random.randn(1, 10)
        filt = np.asarray([1])[:, np.newaxis]
        samples_out = comm.filters.filter_samples(samples, filt, domain='time')
        assert_array_almost_equal(samples, samples_out, decimal=10)

    def test_freq_domain(self):
        samples = np.random.randn(1, 10)
        filt = np.full_like(samples, 1)
        samples_out = comm.filters.filter_samples(samples, filt, domain='freq')
        assert_array_almost_equal(samples, samples_out, decimal=10)
        
class TestRaisedCosineFilter(unittest.TestCase):
    """
    Test class for raised_cosine_filter
    """
    
    def test_freq_domain(self):
        samples = np.ravel(np.asarray([1,0,0,1,1,0,1,0]))
        samples_out = comm.filters.raised_cosine_filter(samples,roll_off=1.0, domain = 'freq')
        


class TestMovingAverage(unittest.TestCase):
    """
    Test class for moving average filter
    """

    def test_dirac_even_freq(self):
        sig_in = np.array([1.0, 0, 0, 0, 0, 0, 0, 0])
        sig_out = comm.filters.moving_average(sig_in, average=4, domain='freq')
        assert_array_almost_equal(sig_out, np.array(
            [0.25, 0.25, 0.0, 0.0, 0.0, 0.0, 0.25, 0.25]), decimal=6)

    def test_dirac_odd_freq(self):
        sig_in = np.array([1.0, 0, 0, 0, 0, 0, 0, 0])
        sig_out = comm.filters.moving_average(sig_in, average=5, domain='freq')
        assert_array_almost_equal(sig_out, np.array(
            [0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.2, 0.2]), decimal=6)

    def test_dirac_even_time(self):
        sig_in = np.array([1.0, 0, 0, 0, 0, 0, 0, 0])
        sig_out = comm.filters.moving_average(sig_in, average=4, domain='time')
        assert_array_almost_equal(sig_out, np.array(
            [0.25, 0.25, 0.25, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), decimal=6)

    def test_dirac_odd_time(self):
        sig_in = np.array([1.0, 0, 0, 0, 0, 0, 0, 0])
        sig_out = comm.filters.moving_average(sig_in, average=5, domain='time')
        assert_array_almost_equal(sig_out, np.array(
            [0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), decimal=6)

class TestIdealLP(unittest.TestCase):
    """
    Test class for ideal_lp function
    """
    
    def test_dirac(self):
        samples = np.ravel(np.zeros((1,1000),dtype='complex'))
        samples[0] = 1
        samples_out = comm.filters.ideal_lp(samples,fc=1)['samples_out']
        assert_array_almost_equal(samples,samples_out,decimal=10)
        
        
class TestTimeShift(unittest.TestCase):
    """
    Test class for time_shift function
    """
    
    def full_roll(self):
        samples = np.random.randn(10)
        tau = len(samples)
        samples_shifted = comm.filters.time_shift(samples,tau=tau)
        assert_array_equal(samples,samples_shifted)