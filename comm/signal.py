import numpy as np
import matplotlib.pyplot as plt
from . import visualizer
from . import tx
from . import utils
from . import channel
from . import filters
from . import rx


class Signal():
    """
    Overall Signal definition.

    samples: list of ndarrays, list of length n_dims, each element containing
    a complex ndarray of size (nsamples,) representing the complex
    samples of the signal
    
    center_frequency: list of scalars of length n_dims, float, [Hz]
    
    sample_rate: list of scalars of length n_dims, float, [Hz], positive
    
    bits: list of ndarrays, list of length n_dims, each element
    containing an ndarray of size (nbits,) representing the logical
    binary information per complex signal dimension
    
    symbols: list of ndarrays, complex, list of length n_dims, each element
    containing an ndarray of size (nsymbols,) representing the complex
    modulation symbols per complex signal dimension
    
    symbol_rate: list of scalars of length n_dims, float, [Hz], positive
    
    modulation_info: list of stings, list of length n_dims, each element
    containing a descriptive name of the used constellation per complex
    signal dimension
    
    constellation: list of ndarrays, list of length n_dims, each element containing
    a complex ndarray of size (nConstellationPoints,) representing representing the
    complex modulation symbols, while the order specifies the mapping
    between bits and modulation symbols (see comm.tx.mapper() for details)

    """

    def __init__(self, n_dims=1):
        """
        Initialize signal structure.

        Initialize signal with n_dims, while all signal dimensions have equal
        default values.

        Parameters
        ----------
        n_dims : int, optional
            Number of complex dimensions of the signal. The default is 1.

        Returns
        -------
        None.

        """
        self.n_dims = n_dims
        self.samples = [np.empty(0, dtype=complex)] * n_dims
        self.center_frequency = [0.0] * n_dims
        self.sample_rate = [1.0] * n_dims
        self.bits = [np.empty(0, dtype=bool)] * n_dims
        self.symbols = [np.empty(0, dtype=complex)] * n_dims
        self.symbol_rate = [1.0] * n_dims
        self.modulation_info = [''] * n_dims
        self.constellation = [np.empty(0, dtype=complex)] * n_dims


    def _check_attribute(self, value):
        """
        Check if attribute is of valid type (or can be converted to a valid type).

        Attribute can be either be:
            * a list of lenth n_dims:
                -> set attribute of every signal dimension accordingly
            * an integer, float, string, ndarray or None of dimension 1:
                -> set attributes of all signal dimensions to this single value
            * an ndarray containing n_dims rows:
                -> set attribute of each dimension to one row of ndarray

        Parameters
        ----------
        value : list, integer, float, string, ndarray, None
            The value to be set in the signal structure.

        Returns
        -------
        value : list
            The value as a list of length n_dims.

        """

        # check for list...
        if isinstance(value, list):
            # ...and correct dimension
            if len(value) == self.n_dims:
                # simple case
                value = value
            else:
                raise ValueError('Signal attributes have to be lists of length n_dims...');
        # try to convert to list
        else:
            # ...for arrays
            if isinstance(value, np.ndarray):
                if (value.ndim == 1):
                    # set all dimensions at once by generating list of correct
                    # length having the ndarray in each element
                    value = [value] * self.n_dims
                elif ((value.ndim == 2) and (value.shape[0] == self.n_dims)):
                    # generate a list in which ever entry contains one row of
                    # the given ndarray
                    value = list(value)
                else:
                    raise ValueError('attribute has to be a ndarray of dimension 1 or has to be of shape (n_dims,X)...')
            # ...for single vlaues                    
            elif (isinstance(value, (int, float, str, bool)) or (value == None)):
                # set all dimensions at once by generating list of correct
                # length form salar integers, floats, strings, bool or None
                value = [value] * self.n_dims
            else:
                raise TypeError('Cannot reasonably convert attribute type to list...')

        return value


    @property
    def n_dims(self):
        return self._n_dims

    @n_dims.setter
    def n_dims(self, value):
        if (not isinstance(value,int)) or (value < 1):
            raise ValueError('n_dims must be an integer and >= 1...')
        self._n_dims = value

    @property
    def samples(self):
        return self._samples

    @samples.setter
    def samples(self, value):
        value = self._check_attribute(value)
        self._samples = value

    @property
    def bits(self):
        return self._bits

    @bits.setter
    def bits(self, value):
        value = self._check_attribute(value)
        self._bits = value

    @property
    def center_frequency(self):
        return self._center_frequency

    @center_frequency.setter
    def center_frequency(self, value):
        value = self._check_attribute(value)
        self._center_frequency = value

    @property
    def sample_rate(self):
        return self._sample_rate

    @sample_rate.setter
    def sample_rate(self, value):
        value = self._check_attribute(value)
        self._sample_rate = value

    @property
    def symbols(self):
        return self._symbols

    @symbols.setter
    def symbols(self, value):
        value = self._check_attribute(value)
        self._symbols = value

    @property
    def symbol_rate(self):
        return self._symbol_rate

    @symbol_rate.setter
    def symbol_rate(self, value):
        value = self._check_attribute(value)
        self._symbol_rate = value

    @property
    def modulation_info(self):
        return self._modulation_info

    @modulation_info.setter
    def modulation_info(self, value):
        value = self._check_attribute(value)
        self._modulation_info = value

    @property
    def constellation(self):
        return self._constellation

    @constellation.setter
    def constellation(self, value):
        value = self._check_attribute(value)
        self._constellation = value



    def generate_bits(self, n_bits=2**15, type='random', seed=None):
        """
        Generate an array of size (n_bits,) binary values.

        For detailed documentation see comm.tx.generate_bits.     
        """
        n_bits = self._check_attribute(n_bits)
        type = self._check_attribute(type)
        seed = self._check_attribute(seed)

        for i, (b, t, s) in enumerate(zip(n_bits, type, seed)):
            self.bits[i] = tx.generate_bits(n_bits=b, type=t, seed=s)


    def set_snr(self, snr_dB=10, seed=None):
        """
        Set the SNR of the signal.
        
        For detailed documentation see comm.channel.set_snr.       
        """
        snr_dB = self._check_attribute(snr_dB)
        seed = self._check_attribute(seed)


        for i, (sn, se) in enumerate(zip(snr_dB, seed)):
            sps = self.sample_rate[i] / self.symbol_rate[i]
            self.samples[i] = channel.set_snr(self.samples[i], sn, sps, se)


    def mapper(self):
        """
        Generate sig.symbols from sig.bits and sig.constellation.

        For detailed documentation see comm.tx.mapper.     
        """
        for i, (b, c) in enumerate(zip(self.bits, self.constellation)):
            self.symbols[i] = tx.mapper(bits=b, constellation=c)
    
    def demapper(self):
        """
        Demap samples to bits using a given constellation alphabet.
        
        For detailed documentation see comm.rx.demapper.        
        """
        for i, (s, c) in enumerate(zip(self.samples, self.constellation)):
            self.samples[i] = rx.demapper(samples=s, constellation=c)
    

    def decision(self):
        """
        Decide samples to a given constellation alphabet.
        
        For detailed documentation see comm.rx.decision.        
        """
        for i, (s, c) in enumerate(zip(self.samples, self.constellation)):
            self.samples[i] = rx.decision(samples=s, constellation=c)


    def raised_cosine_filter(self, roll_off=0.0, root_raised=False, **kargs):
        """
        Filter samples with a raised cosine filter.

        For detailed documentation see comm.filters.raised_cosine_filter.      
        """
        roll_off = self._check_attribute(roll_off)
        root_raised = self._check_attribute(root_raised)

        for i, (s, sr, symr, ro, rr) in enumerate(zip(self.samples,
                                                      self.sample_rate,
                                                      self.symbol_rate,
                                                      roll_off, root_raised)):
            self.samples[i] = filters.raised_cosine_filter(samples=s,
                                                           sample_rate=sr,
                                                           symbol_rate=symr,
                                                           roll_off=ro,
                                                           root_raised=rr, **kargs)
            
    def sampling_phase_adjustment(self):
        """
        Estimate the sampling phase offset and compensate for it.
        
        For detailed documentation see comm.rx.sampling_phase_adjustment.
        """
        for i, (s, sr, symr) in enumerate(zip(self.samples,
                                                      self.sample_rate,
                                                      self.symbol_rate)):
            
            self.samples[i] = rx.sampling_phase_adjustment(samples=s,
                                                           sample_rate=sr,
                                                           symbol_rate=symr)['samples_out']
            
    def sampling_clock_adjustment(self, block_size=500):
        """
        Estimate the sampling clock offset and compensate for it.
        
        For detailed documentation see comm.rx.sampling_clock_adjustment.
        """
        block_size = self._check_attribute(block_size)
        
        for i, (s, sr, symr, bs) in enumerate(zip(self.samples,
                                                      self.sample_rate,
                                                      self.symbol_rate, block_size)):
            
            self.samples[i] = rx.sampling_clock_adjustment(samples=s,
                                                           sample_rate=sr,
                                                           symbol_rate=symr,
                                                           block_size=bs)['samples_out']
        
        


    def generate_constellation(self, format='QAM', order=4):
        """
        Set sig.constellation and sig.modulation_info.

        For detailed documentation see comm.utils.generate_constellation.    
        """
        format = self._check_attribute(format)
        order = self._check_attribute(order)

        for i, (f, o) in enumerate(zip(format, order)):
            self.constellation[i] = utils.generate_constellation(format=f, order=o)
            self.modulation_info[i] = str(o) + "-" + str(f)


    def pulseshaper(self, upsampling=2.0, pulseshape='rc', roll_off=0.2):
        """
        Upsample and pulseshape the modulated symbols and write them to samples.

        For detailed documentation see comm.tx.pulseshaper
        """
        upsampling = self._check_attribute(upsampling)
        pulseshape = self._check_attribute(pulseshape)
        roll_off = self._check_attribute(roll_off)


        for i, (u, p, r) in enumerate(zip(upsampling, pulseshape, roll_off)):
            self.samples[i] = tx.pulseshaper(self.symbols[i], u, p, r)
            self.sample_rate[i] = u * self.symbol_rate[i]



    def plot_spectrum(self, dimension=0, **kwargs):
        """
        Plot spectum of the signal samples of a given dimension.

        For further documentation see comm.visualizer.plot_spectrum.
        """
        visualizer.plot_spectrum(self.samples[dimension], self.sample_rate[dimension], **kwargs)

    def plot_constellation(self, dimension=0, decimation=1, **kwargs):
        """
        Plot constellation of signal samples of a given dimension.

        For further documentation see comm.visualizer.plot_constellation.
        """
        visualizer.plot_constellation(self.samples[dimension], decimation=decimation, **kwargs)


    def plot_eye(self, dimension=0, offset=0, **kwargs):
        """
        Plot eye diagramm of signal samples of a given dimension.

        For further documentation see comm.visualizer.plot_eye
        """

        visualizer.plot_eye(self.samples[dimension], self.sample_rate[dimension],
                            self.symbol_rate[dimension], offset, **kwargs)
