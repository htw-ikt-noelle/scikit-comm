import numpy as np
import matplotlib.pyplot as plt
from . import visualizer
from . import tx
from . import utils

       
class Signal():
    """ Overall Signal definition.        
              
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
        
        if (not isinstance(n_dims,int)) or (n_dims < 1):
            raise ValueError('n_dims must be an integer and >= 1...')
        
        self.samples = [np.empty(0, dtype=complex)] * n_dims
        self.center_frequency = [0.0] * n_dims
        self.sample_rate = [1.0] * n_dims
        self.bits = [np.empty(0, dtype=bool)] * n_dims
        self.symbols = [np.empty(0, dtype=complex)] * n_dims
        self.symbol_rate = [1.0] * n_dims
        self.modulation_info = [''] * n_dims
        self.constellation = [np.empty(0, dtype=complex)] * n_dims
        
    def generate_bits(self, n_bits=[2**15], type=['random'], seed=[None]):
        """
        Generate an array of size (n_bits,) binary values.

        Parameters
        ----------
        n_bits : TYPE, optional
            DESCRIPTION. The default is [2**15].
        type : TYPE, optional
            DESCRIPTION. The default is ['random'].
        seed : TYPE, optional
            DESCRIPTION. The default is [None].

        Raises
        ------
        TypeError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        n_dims = len(self.bits)
        
        if not (isinstance(n_bits,list) and isinstance(type,list) and isinstance(seed,list)):
            raise TypeError('input parameters have to be lists...')
            
        if not ((len(n_bits) == 1) or (len(n_bits) == n_dims)):
            raise TypeError('length of n_bits should be 1 or n_dims...')
        elif len(n_bits) == 1:
            n_bits = n_bits * n_dims            
            
        if not ((len(type) == 1) or (len(type) == n_dims)):
            raise TypeError('length of type should be 1 or n_dims...')
        elif len(type) == 1:
            type = type * n_dims
            
        if not ((len(seed) == 1) or (len(seed) == n_dims)):
            raise TypeError('length of seed should be 1 or n_dims...')
        elif len(seed) == 1:
            seed = seed * n_dims    

        for i, (b, t, s) in enumerate(zip(n_bits, type, seed)):
            self.bits[i] = tx.generate_bits(n_bits=b, type=t, seed=s)
            
            
            
            
    def set_snr(self, snr_dB=[10], seed=[None]):
        
        
        n_dims = len(self.bits)
        
        if not (isinstance(snr_dB,list) and isinstance(seed,list)):
            raise TypeError('input parameters have to be lists...')
            
        if not ((len(snr_dB) == 1) or (len(snr_dB) == n_dims)):
            raise TypeError('format of snr_dB should be 1 or n_dims...')
        elif len(snr_dB) == 1:
            snr_dB = snr_dB * n_dims            
            
        if not ((len(seed) == 1) or (len(seed) == n_dims)):
            raise TypeError('seed of type should be 1 or n_dims...')
        elif len(seed) == 1:
            seed = seed * n_dims
            
        for i, (sn, se) in enumerate(zip(snr_dB, seed)):
            sps = self.sample_rate[i] / self.symbol_rate[i]
            self.samples[i] = channel.set_snr(self.samples[i], sn, sps, se)        
            
            
    def mapper(self):
        """
        Generate sig.symbols from sig.bits and sig.constellation.

        Returns
        -------
        None.

        """
        
        for i, (b, c) in enumerate(zip(self.bits, self.constellation)):
            self.symbols[i] = tx.mapper(bits=b, constellation=c)
            
            
    def generate_constellation(self, format=['QAM'], order=[4]):
        """
        Set sig.constellation and sig.modulation_info.

        Parameters
        ----------
        format : TYPE, optional
            DESCRIPTION. The default is ['QAM'].
        order : TYPE, optional
            DESCRIPTION. The default is [4].

        Raises
        ------
        TypeError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        n_dims = len(self.bits)
        
        if not (isinstance(format,list) and isinstance(order,list)):
            raise TypeError('input parameters have to be lists...')
            
        if not ((len(format) == 1) or (len(format) == n_dims)):
            raise TypeError('format of n_bits should be 1 or n_dims...')
        elif len(format) == 1:
            format = format * n_dims            
            
        if not ((len(order) == 1) or (len(order) == n_dims)):
            raise TypeError('order of type should be 1 or n_dims...')
        elif len(order) == 1:
            order = order * n_dims
            
        for i, (f, o) in enumerate(zip(format, order)):
            self.constellation[i] = utils.generate_constellation(format=f, order=o)
            self.modulation_info[i] = f
            
        
    def pulseshaper(self, upsampling=[2], pulseshape=['rc'], roll_off=[0.2]):
        """
        Upsample and pulseshape the modulated symbols and write them to samples.
        
        For detailed documentation see comm.tx.pulseshaper()

        Parameters
        ----------
        upsampling : TYPE, optional
            DESCRIPTION. The default is [2].
        pulseshape : TYPE, optional
            DESCRIPTION. The default is ['rc'].
        roll_off : TYPE, optional
            DESCRIPTION. The default is [0.2].

        Returns
        -------
        None.

        """

        n_dims = len(self.symbols)
        
        if not (isinstance(upsampling,list) and isinstance(pulseshape,list) and
                isinstance(roll_off,list)):
            raise TypeError('input parameters have to be lists...')
            
        if not ((len(upsampling) == 1) or (len(upsampling) == n_dims)):
            raise TypeError('upsampling of n_bits should be 1 or n_dims...')
        elif len(upsampling) == 1:
            upsampling = upsampling * n_dims            
            
        if not ((len(pulseshape) == 1) or (len(pulseshape) == n_dims)):
            raise TypeError('pulseshape of type should be 1 or n_dims...')
        elif len(pulseshape) == 1:
            pulseshape = pulseshape * n_dims
            
        if not ((len(roll_off) == 1) or (len(roll_off) == n_dims)):
            raise TypeError('roll_off of type should be 1 or n_dims...')
        elif len(roll_off) == 1:
            roll_off = roll_off * n_dims
            
        for i, (u, p, r) in enumerate(zip(upsampling, pulseshape, roll_off)):
            self.samples[i] = tx.pulseshaper(self.symbols[i], u, p, r)
            self.sample_rate[i] = u
            
            
            
    def mapper(self):
        """
        Generate sig.symbols from sig.bits and sig.constellation.

        Returns
        -------
        None.

        """
        
        for i, (b, c) in enumerate(zip(self.bits, self.constellation)):
            self.symbols[i] = tx.mapper(bits=b, constellation=c)
        


    def plot_spectrum(self, dimension=0):
        visualizer.plot_spectrum(self.samples[dimension], self.sample_rate[dimension])
