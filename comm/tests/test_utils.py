import unittest
import numpy as np
from numpy.testing import (
        assert_, assert_raises, assert_equal, assert_almost_equal, assert_warns,
        assert_no_warnings, assert_array_equal, assert_array_almost_equal,
        suppress_warnings
        )

import os
import sys
sys.path.insert(0, os.path.abspath('..\..'))
import comm as comm


class TestGenerateConstellation(unittest.TestCase):
    """
    Test class for generate_constellation
    """
    
    def test_completeness(self):
        # QPSK
        ref = set([-1.-1.j, -1.+1.j,  1.-1.j,  1.+1.j])
        const = comm.utils.generate_constellation(format='QAM', order=4)
        assert_equal(len(ref.difference(const)), 0)
        # 16QAM
        ref = set([-3.-3.j, -3.-1.j, -3.+3.j, -3.+1.j, -1.-3.j, -1.-1.j, -1.+3.j,
                   -1.+1.j,  3.-3.j,  3.-1.j,  3.+3.j,  3.+1.j,  1.-3.j,  1.-1.j,
                   1.+3.j,  1.+1.j])
        const = comm.utils.generate_constellation(format='QAM', order=16)
        assert_equal(len(ref.difference(const)), 0)
        # 32 QAM        
        ref = set([-3.-5.j, -1.-5.j, -3.+5.j, -1.+5.j, -5.-3.j, -5.-1.j, -5.+3.j,
                   -5.+1.j, -1.-3.j, -1.-1.j, -1.+3.j, -1.+1.j, -3.-3.j, -3.-1.j,
                   -3.+3.j, -3.+1.j,  3.-5.j,  1.-5.j,  3.+5.j,  1.+5.j,  5.-3.j,
                   5.-1.j,  5.+3.j,  5.+1.j,  1.-3.j,  1.-1.j,  1.+3.j,  1.+1.j,
                   3.-3.j,  3.-1.j,  3.+3.j,  3.+1.j])
        const = comm.utils.generate_constellation(format='QAM', order=32)
        assert_equal(len(ref.difference(const)), 0)
        
    def test_errors(self):
        # verifies errors raised by function
        with assert_raises(TypeError):
            # order type
            comm.utils.generate_constellation(order=4.)  
        with assert_raises(ValueError):
            # valid order numbers
            comm.utils.generate_constellation(order=5)
            comm.utils.generate_constellation(order=2)
            # valid modulation formats
            comm.utils.generate_constellation(format='DPSK')
        
        
class TestTimeAxis(unittest.TestCase):
    """
    Test class for creat_time_axis
    """
    
    def test_length(self):
        l = 100
        t = comm.utils.create_time_axis(sample_rate=1.0, n_samples=l)
        assert_equal(len(t), l)
        l = 101
        t = comm.utils.create_time_axis(sample_rate=1.0, n_samples=l)
        assert_equal(len(t), l)
    
    def test_end_points_and_delta(self):
        l = 100
        sr = 10
        t = comm.utils.create_time_axis(sample_rate=sr, n_samples=l)
        assert_equal(t[0], 0)
        assert_equal(t[-1], 1/sr*(l-1))
        assert_equal(t[1]-t[0], 1/sr)
        
        l = 101
        sr = 10
        t = comm.utils.create_time_axis(sample_rate=sr, n_samples=l)
        assert_equal(t[0], 0)
        assert_equal(t[-1], 1/sr*(l-1))
        assert_equal(t[1]-t[0], 1/sr)
        
        
class TestBitsToDec(unittest.TestCase):
    """
    Test class for bits_to_dec
    """
    
    def test_errors(self):
        # test errors thrown by function
        with assert_raises(ValueError):
            # one-dimensionality
            comm.utils.bits_to_dec(np.asarray([[True,False],[False,True]]),4)
            # n_bits is integer multiple of m
            comm.utils.bits_to_dec(np.asarray([True,False,True,True]),3)
    
    def test_reconversion(self):
        # convert and reconvert bit sequence and test for equality
        dummy = comm.tx.generate_bits(n_bits=20)
        dec_dummy = comm.utils.bits_to_dec(dummy,4)
        bits = comm.utils.dec_to_bits(dec_dummy,4)
        assert_equal(bits,dummy)
        
class TestEstimateSnrSpectrum(unittest.TestCase):
    """
    Test class for estimate_snr_spectrum
    """
    
    def test_correct_snr(self):
        
        x = np.arange(-10, 11, 1)
        # noise = 0.5
        y = np.ones_like(x) * 0.5
        #  signal + noise = 6
        y[int(x.size/2)-1:int(x.size/2)+2] = 6
        noise_range = np.asarray([-4.2, -8.0, 4.1, 8.4])
        sig_range = np.asarray([-0.9, 0.6])
        noise_bw = 1
        order = 1

        snr = comm.utils.estimate_snr_spectrum(x, y, sig_range, noise_range, 
                                               order=order, noise_bw=noise_bw, 
                                               scaling='lin', plotting=False)
        
        assert_almost_equal(snr, 10.0, decimal=10)
    
        
        