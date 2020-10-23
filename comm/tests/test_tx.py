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
        
        
