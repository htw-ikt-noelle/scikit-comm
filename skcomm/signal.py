import copy

import numpy as np

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
    between bits and modulation symbols (see skcomm.tx.mapper() for details)

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


    def __iter__(self):
        for attr, value in vars(self).items():
            yield attr, value

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
                    value = [copy.deepcopy(value) for dim in range(self.n_dims)]
                elif ((value.ndim == 2) and (value.shape[0] == self.n_dims)):
                    # generate a list in which every entry contains one row of
                    # the given ndarray
                    value = list(value)
                else:
                    raise ValueError('attribute has to be a ndarray of dimension 1 or has to be of shape (n_dims,X)...')
            # ...for single vlaues                    
            elif (isinstance(value, (int, float, str, bool)) or (value == None)):
                # set all dimensions at once by generating list of correct
                # length form salar integers, floats, strings, bool or None
                value = [copy.deepcopy(value) for dim in range(self.n_dims)]
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

    def get_dimensions(self, dims=[0]):
        """
        Gets all attributes form specific signal dimensions.

        This method returns a new skcomm.signal.Signal object
        containing all attributes of the signal dimensions specified by dims.

        This method could also be used to reorder (or duplicate) the dimensions 
        of a signal object. 
        Examples: 
        1)
        signal.get_dimensions(dims=[0,0]) results in a two-dimensional 
        signal object containing two identical signal dimensions.

        2) signal.get_dimensions(dims=[1,0]) reoders the dimensions of
        a two-dimensional signal object.
        
        Parameters
        ----------
        dims : list of int
            Specifies which signal dimensions should be returned. The defauls is [0] 
        
        Returns
        -------
        tmp : skcomm.signal.Signal
            Signal containing all attributes of specified signal dimensions.
        """
        if not isinstance(dims,list):
            raise ValueError('dims needs to be of type list')

        if (len(dims) > self.n_dims) or (len(dims)==0):
            raise ValueError(f'ndims needs to be a list of length between 1 and {self.n_dims}')
        
        if any(np.asarray(dims)<0) or any(np.asarray(dims)>(self.n_dims-1)):
            raise ValueError(f'The requested dimensions need to be between 0 and {self.n_dims-1}')
        
        tmp = Signal(n_dims=len(dims))
        for sour_attr, sour_value in self:
            if not (sour_attr == '_n_dims'):                
                vars(tmp)[sour_attr] = [sour_value[i] for i in dims]
        return tmp
    
    def add_dimension(self, sig, dim=0):
        """
        Adds a signal dimension to the signal.

        The content of a one dimensional signal object is added to the existing
        signal object. The position at which the additional dimension will be added can be 
        specified with the parameter dim.
        
        Parameters
        ----------
        sig : skc.signal.Signal
            One dimensional signal object which will be added to the signal.
        dim : int
            Specifies at which position the additional dimension will be inserted. dim is
            the index of the dimension before which to insert, so dim=0 inserts the dimension
            at the front while dim=self.n_dims appends the dimension at the end of the signal
            object. The defauls is 0.
        
        """
        if not isinstance(dim,int):
            raise ValueError('dim needs to be an integer')
        
        if not isinstance(sig,Signal):
            raise ValueError('sig needs to be of type skc.signal.Signal')
        
        if sig.n_dims != 1:
            raise ValueError('signal needs to be a one-dimensional signal')        
        
        self.n_dims += 1

        for sour_attr, sour_value in sig:
            if not (sour_attr == '_n_dims'):                                                   
                vars(self)[sour_attr].insert(dim,sour_value[0])
        


    def generate_bits(self, n_bits=2**15, type='random', seed=None):
        """
        Generate an array of shape (n_bits,) binary values.

        For detailed documentation see skcomm.tx.generate_bits.     
        """
        n_bits = self._check_attribute(n_bits)
        type = self._check_attribute(type)
        seed = self._check_attribute(seed)

        for i, (b, t, s) in enumerate(zip(n_bits, type, seed)):
            self.bits[i] = tx.generate_bits(n_bits=b, type=t, seed=s)


    def set_snr(self, snr_dB=10, seed=None):
        """
        Set the SNR of the signal.
        
        For detailed documentation see skcomm.channel.set_snr.       
        """
        snr_dB = self._check_attribute(snr_dB)
        seed = self._check_attribute(seed)


        for i, (sn, se) in enumerate(zip(snr_dB, seed)):
            sps = self.sample_rate[i] / self.symbol_rate[i]
            self.samples[i] = channel.set_snr(self.samples[i], sn, sps, se)


    def mapper(self):
        """
        Generate sig.symbols from sig.bits and sig.constellation.

        For detailed documentation see skcomm.tx.mapper.     
        """
        for i, (b, c) in enumerate(zip(self.bits, self.constellation)):
            self.symbols[i] = tx.mapper(bits=b, constellation=c)
    
    def demapper(self):
        """
        Demap samples to bits using a given constellation alphabet.
        
        For detailed documentation see skcomm.rx.demapper.        
        """
        for i, (s, c) in enumerate(zip(self.samples, self.constellation)):
            self.samples[i] = rx.demapper(samples=s, constellation=c)
    

    def decision(self):
        """
        Decide samples to a given constellation alphabet.
        
        For detailed documentation see skcomm.rx.decision.        
        """
        for i, (s, c) in enumerate(zip(self.samples, self.constellation)):
            self.samples[i] = rx.decision(samples=s, constellation=c)


    def raised_cosine_filter(self, roll_off=0.0, root_raised=False, **kargs):
        """
        Filter samples with a raised cosine filter.

        For detailed documentation see skcomm.filters.raised_cosine_filter.      
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
        
        For detailed documentation see skcomm.rx.sampling_phase_adjustment.
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
        
        For detailed documentation see skcomm.rx.sampling_clock_adjustment.
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

        For detailed documentation see skcomm.utils.generate_constellation.    
        """
        format = self._check_attribute(format)
        order = self._check_attribute(order)

        for i, (f, o) in enumerate(zip(format, order)):
            self.constellation[i] = utils.generate_constellation(format=f, order=o)
            self.modulation_info[i] = str(o) + "-" + str(f)


    def pulseshaper(self, upsampling=2.0, pulseshape='rc', roll_off=0.2):
        """
        Upsample and pulseshape the modulated symbols and write them to samples.

        For detailed documentation see skcomm.tx.pulseshaper
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

        For further documentation see skcomm.visualizer.plot_spectrum.
        """
        results = visualizer.plot_spectrum(self.samples[dimension], sample_rate=self.sample_rate[dimension], **kwargs)
        return results

    def plot_constellation(self, dimension=0, decimation=1, **kwargs):
        """
        Plot constellation of signal samples of a given dimension.

        For further documentation see skcomm.visualizer.plot_constellation.
        """
        visualizer.plot_constellation(self.samples[dimension], decimation=decimation, **kwargs)


    def plot_eye(self, dimension=0, boundaries=[None, None], **kwargs):
        """
        Plot eye diagramm of signal samples of a given dimension.

        For further documentation see skcomm.visualizer.plot_eye
        """

        visualizer.plot_eye(self.samples[dimension], self.sample_rate[dimension],
                            self.symbol_rate[dimension],
                            boundaries=boundaries, **kwargs)
        
        
    def plot_signal(self, dimension=0, boundaries=[None, None], **kwargs):
        """
        Plot the signal samples of a given dimension as a function of time.

        For further documentation see skcomm.visualizer.plot_signal
        """
        
        visualizer.plot_signal(self.samples[dimension], self.sample_rate[dimension],
                               boundaries=boundaries, **kwargs)
