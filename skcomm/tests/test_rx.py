import unittest
import numpy as np
from numpy.testing import (
        assert_, assert_raises, assert_equal, assert_warns,
        assert_no_warnings, assert_array_equal, assert_array_almost_equal,
        suppress_warnings
        )

from .. import rx
from .. import utils
from .. import signal


class TestDemapper(unittest.TestCase):
    """
    Test class for demapper
    """
    
    def test_four_point_const(self):
        bits = np.asarray([0,0,0,1,1,1,1,0]) # 0 1 3 2
        samples = np.asarray([1+0j, 0+1j, 0-1j, -1+0j])
        constellation = np.asarray([1+0j, 0+1j, -1+0j, 0-1j])
        demapped = rx.demapper(samples, constellation)
        assert_equal(demapped, bits)
        
class TestDecision(unittest.TestCase):
    """
    Test class for decision
    """
    
    def test_four_point_const(self):
        samples = np.asarray([2+0j, 1+0.99j, 0.01-0j, 0.5-0.49j])
        constellation = np.asarray([1+0j, 0+1j, -1+0j, 0-1j])
        result = np.asarray([1+0j, 1+0j, 1+0j, 1+0j]) 
        dec = rx.decision(samples, constellation)
        assert_equal(result, dec)
        
    def test_scaling(self):
        constellation = np.asarray([1,2,3,4])
        samples1 = constellation*4
        samples2 = constellation*0.1
        result = constellation
        dec1 = rx.decision(samples1, constellation)
        dec2 = rx.decision(samples2, constellation)
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
        t = utils.create_time_axis(sample_rate=10.0, n_samples=n_samples)
        samples = np.cos(2*np.pi*symbr/2*t + phase)
        shift = rx.sampling_phase_adjustment(samples, sample_rate=sr, symbol_rate=symbr, shift_dir='delay')['est_shift']
        assert_array_almost_equal(shift, -phase/np.pi/symbr, decimal=10)
        shift = rx.sampling_phase_adjustment(samples, sample_rate=sr, symbol_rate=symbr, shift_dir='advance')['est_shift']
        assert_array_almost_equal(shift, 1.0-(phase/np.pi/symbr), decimal=10)
        shift = rx.sampling_phase_adjustment(samples, sample_rate=sr, symbol_rate=symbr, shift_dir='both')['est_shift']
        assert_array_almost_equal(shift, -phase/np.pi/symbr, decimal=10)
        shifts = rx.sampling_clock_adjustment(samples, sample_rate=sr, symbol_rate=symbr, block_size=int(n_samples/(sr/symbr*2)))['est_shift']
        assert_array_almost_equal(shifts, np.asarray([-phase/np.pi/symbr, -phase/np.pi/symbr]), decimal=10)
        

class TestSymbolSequenceSync(unittest.TestCase):
    """
    Test class for symbol_sequence_sync
    """
    
    def test_delays_and_phase_shifts(self):
        # generate test signal
        sig = signal.Signal(n_dims=2)
        sig.constellation = [np.asarray([-1, 1]), np.asarray([-1-1j,-1+1j, 1-1j, 1+1j])]
        sig.symbols = [sig.constellation[0][np.random.randint(low=0,high=2,size=100)], sig.constellation[1][np.random.randint(low=0,high=4,size=100)]]
        # 1) no shift, no phase error
        sig.samples = sig.symbols
        tmp = rx.symbol_sequence_sync(sig, dimension=-1)        
        assert_equal(tmp[0],
                {'symbol_delay_est':0, 'phase_est':0.0}
                )
        assert_equal(tmp[1],
                {'symbol_delay_est':0, 'phase_est':-0.0}
                )
        # 2) shift, shift+phase error
        sig.samples = [np.roll(sig.symbols[0],shift=-10), np.roll(sig.symbols[1], shift=-5)*np.exp(1j*np.pi)]
        tmp = rx.symbol_sequence_sync(sig, dimension=-1)        
        assert_equal(tmp[0],
                {'symbol_delay_est':10, 'phase_est':0.0}
                )
        assert_equal(tmp[1],
                {'symbol_delay_est':5, 'phase_est':-np.pi}
                )
        # 3) shift+phase error, shift+conj+phase phase error
        sig.samples = [np.roll(sig.symbols[0]*np.exp(-1j*np.pi/2),shift=-10), np.conj(np.roll(sig.symbols[1], shift=-10)*np.exp(1j*np.pi))]
        tmp = rx.symbol_sequence_sync(sig, dimension=-1)        
        assert_equal(tmp[0],
                {'symbol_delay_est':10, 'phase_est':np.pi/2}
                )
        assert_equal(tmp[1],
                {'symbol_delay_est':10, 'phase_est':np.pi}
                )

