import numpy as np
from numpy.testing import (
        assert_, assert_raises, assert_equal, assert_warns,
        assert_no_warnings, assert_array_equal, assert_array_almost_equal,
        suppress_warnings
        )

# TODO: CHECK IMPORT!!!!!
from ... import comm


class TestSetSNR:
    
    def test_high_snr(self):
        sig_in = np.array([1, 1, 1, 1])
        sig_out = comm.channel.set_snr(sig_in, snr_dB=100)
        assert_array_almost_equal(sig_in, sig_out, decimal=10)
        

