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


class TestDemapper(unittest.TestCase):
    """
    Test class for demapper
    """
    
    def test_four_point_const(self):
        bits = np.asarray([0,0,0,1,1,1,1,0]) # 0 1 3 2
        samples = np.asarray([1+0j, 0+1j, 0-1j, -1+0j])
        constellation = np.asarray([1+0j, 0+1j, -1+0j, 0-1j])
        demapped = comm.rx.demapper(samples, constellation)
        assert_equal(demapped, bits)
        
class TestDecision(unittest.TestCase):
    """
    Test class for decision
    """
    
    def test_four_point_const(self):
        samples = np.asarray([2+0j, 1+0.99j, 0.01-0j, 0.5-0.49j])
        constellation = np.asarray([1+0j, 0+1j, -1+0j, 0-1j])
        result = np.asarray([1+0j, 1+0j, 1+0j, 1+0j]) 
        dec = comm.rx.decision(samples, constellation)
        assert_equal(result, dec)
        
    def test_scaling(self):
        constellation = np.asarray([1,2,3,4])
        samples1 = constellation*4
        samples2 = constellation*0.1
        result = constellation
        dec1 = comm.rx.decision(samples1, constellation)
        dec2 = comm.rx.decision(samples2, constellation)
        assert_equal(result, dec1)
        assert_equal(result, dec2)
        
class TestSamplingPhaseClockAdjustment(unittest.TestCase):
    """
    Test class for sampling_phase_adjustment and sampling_clock_adjumstment
    """
    
    def test_four_point_const(self):
        sr = 10.0
        symbr = 1.0
        phase = 1.0
        n_samples = 100
        t = comm.utils.create_time_axis(sample_rate=10.0, n_samples=n_samples)
        samples = np.cos(2*np.pi*symbr/2*t + phase)
        shift = comm.rx.sampling_phase_adjustment(samples, sample_rate=sr, symbol_rate=symbr, shift_dir='delay')['est_shift']
        assert_array_almost_equal(shift, -phase/np.pi/symbr, decimal=10)
        shift = comm.rx.sampling_phase_adjustment(samples, sample_rate=sr, symbol_rate=symbr, shift_dir='advance')['est_shift']
        assert_array_almost_equal(shift, 1.0-(phase/np.pi/symbr), decimal=10)
        shift = comm.rx.sampling_phase_adjustment(samples, sample_rate=sr, symbol_rate=symbr, shift_dir='both')['est_shift']
        assert_array_almost_equal(shift, -phase/np.pi/symbr, decimal=10)
        shifts = comm.rx.sampling_clock_adjustment(samples, sample_rate=sr, symbol_rate=symbr, block_size=int(n_samples/(sr/symbr*2)))['est_shift']
        assert_array_almost_equal(shifts, np.asarray([-phase/np.pi/symbr, -phase/np.pi/symbr]), decimal=10)
        
        
        
