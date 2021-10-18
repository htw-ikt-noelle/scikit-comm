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


class TestGenerateBits(unittest.TestCase):
    """
    Test class for generate_bits
    """
    
    def test_dimensions(self):
        bits = comm.tx.generate_bits(n_bits=4, type='random', seed=1)
        assert_equal(bits.ndim, 1)
        assert_equal(bits.shape[0], 4)
        
        
class TestMapper(unittest.TestCase):
    """
    Test class for mapper
    """
    
    def test_four_point_const(self):
        bits = np.asarray([0,0,0,1,1,1,1,0]) # 0 1 3 2
        constellation = np.asarray([1+0j, 0+1j, -1+0j, 0-1j])
        mapped = comm.tx.mapper(bits, constellation)
        assert_equal(mapped, np.asarray([1+0j, 0+1j, 0-1j, -1+0j]))
        