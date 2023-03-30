import unittest
import numpy as np
from numpy.testing import (
        assert_, assert_raises, assert_equal, assert_warns,
        assert_no_warnings, assert_array_equal, assert_array_almost_equal,
        suppress_warnings
        )

from .. import tx


class TestGenerateBits(unittest.TestCase):
    """
    Test class for generate_bits
    """
    
    def test_dimensions(self):
        bits = tx.generate_bits(n_bits=4, type='random', seed=1)
        assert_equal(bits.ndim, 1)
        assert_equal(bits.shape[0], 4)
        
        
class TestMapper(unittest.TestCase):
    """
    Test class for mapper
    """
    
    def test_four_point_const(self):
        bits = np.asarray([0,0,0,1,1,1,1,0]) # 0 1 3 2
        constellation = np.asarray([1+0j, 0+1j, -1+0j, 0-1j])
        mapped = tx.mapper(bits, constellation)
        assert_equal(mapped, np.asarray([1+0j, 0+1j, 0-1j, -1+0j]))
        
    def test_length(self):
        bits = np.asarray([0,0,0,1,1,1,1,0]) # 0 1 3 2
        constellation = np.asarray([1+0j, 0+1j, -1+0j, 0-1j])
        mapped = tx.mapper(bits, constellation)
        # test length of vector of complex symbols
        assert_equal(len(mapped),len(bits)/np.log2(len(constellation)))
        
class TestPulseshaper(unittest.TestCase):
    """
    Test class for pulseshaper
    """
    
    def test_upsampling(self):
        # test correct length of upsampled vector
        samples_dummy = np.random.rand(100)
        up = 2
        # upsample by default value (2) and apply pulseshaping
        rect_dummy = tx.pulseshaper(samples_dummy,upsampling=up,pulseshape='rect')
        rc_dummy = tx.pulseshaper(samples_dummy,upsampling=up,pulseshape='rc')
        rrc_dummy = tx.pulseshaper(samples_dummy,upsampling=up,pulseshape='rrc')
        # test lengths
        assert_equal(len(rect_dummy),len(samples_dummy)*up)
        assert_equal(len(rc_dummy),len(samples_dummy)*up)
        assert_equal(len(rrc_dummy),len(samples_dummy)*up)