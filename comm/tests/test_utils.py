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


class TestGenerateConstelation(unittest.TestCase):
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
        
        
        
        