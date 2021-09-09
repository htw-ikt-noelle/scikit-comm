import numpy as np
import matplotlib.pyplot as plt

def generate_gray_code(nbits):
    """
    'Recursively' generates a Gray-coded bit sequence with [nbits] bits.

    Parameters
    ----------
    nbits : int
        number of bits per element

    Returns
    -------
    gray_alphabet : list of str
        list of Gray-coded binary representations of all numbers from 0 to (2**nbits)-1.

    """
    L1 = ['0','1']
    # generate n-bit Gray code alphabet
    for i in range(1,nbits):
        # reverse L1
        L2 = L1[::-1]
        # prefix 0 to L1 list elements
        L1 = ['0' + i for i in L1]
        # prefix 1 to L2 list elements
        L2 = ['1' + i for i in L2]
        # concatenate lists
        L1 = L1+L2
    
    gray_alphabet = L1
    
    return gray_alphabet

def generate_qam_constellation(order):
    """
    Generates a Gray-coded [order]-QAM constellation.

    Parameters
    ----------
    order : int
        number of symbols in the constellation.

    Returns
    -------
    const_order : list of complex
        list of complex symbols to map bits to.

    """
    # algorithm for generating gray-coded qam constellation alphabet
    
    # first implementation for qam with even powers of 2 (4,16,64,256...)

    # recursive one-dimensional gray-code implementation
    
    # derive number of bits encoded in one symbol from QAM order
    n = int(np.log2(order))
    # generate Gray code with n bits
    L1 = generate_gray_code(n)
    
    # for odd-bit QAM constellations:
    if int(np.log2(order))%2:
        # so far, rectangular QAM constellation is generated
        # TODO: implement symmetrical QAM (see paper by Smith, 1975)
        
        # generate separate Gray codes for I and Q, with [(n/2)+1] bits in I
        # and [n/2] bits in Q branch
        gray_I = [i[:-(-len(i) // 2)] for i in L1]
        gray_Q = [i[-(-len(i) // 2):] for i in L1]
        
        # convert into indices
        gray_I_idx = [int(i,base=2) for i in gray_I]
        gray_Q_idx = [int(i,base=2) for i in gray_Q]
        
        # generate evenly spaced values for I and Q branches
        values_I = 2*np.linspace(-1,1,2**(int(-(-n // 2))))
        values_Q = np.linspace(-1,1,2**(int(n // 2)))
        
        # reorder values with Gray code
        # values_I = values_I[[int(i,base=2) for i in gray_I]]
        # values_Q = values_Q[[int(i,base=2) for i in gray_Q]]
        
        # promote to matrix with np.tile() to facilitate building complex signals
        # values_I = np.tile(values_I,(len(gray_Q),1))
        # values_Q = np.transpose(np.tile(values_Q,(len(gray_I),1)))
        
        # index (reordered) evenly spaced values and build complex symbols
        const = (values_I[gray_I_idx] + 1j*values_Q[gray_Q_idx])
        const = const.flatten() # ATTN: not correct, continue here!
        # const = values_I[gray_I_idx,gray_Q_idx] + 1j*values_Q[gray_I_idx,gray_Q_idx]
        # const_order = const
        
        # bring bits and complex symbols into correct order
        tmp = [int(i,base=2) for i in L1]
            
        constellation = []
        bits = []
        for i in range(0,len(tmp)):
            constellation.append(const[tmp.index(i)])
            bits.append(L1[tmp.index(i)])
        # raise ValueError('Currently, only QAM orders equal to even powers of 2 are supported!')
    
    # for even-bit QAM constellations:
    else:
        # # for even bit qam constellations
        # # split Gray coded list elements (strings) in half
        # L1re = [i[:len(i)//2] for i in L1]
        # L1im = [i[len(i)//2:] for i in L1]
        
        # # convert into indices
        # L1re_idx = [int(i,base=2) for i in L1re]
        # L1im_idx = [int(i,base=2) for i in L1im]
        
        # # generate evenly spaced values
        # values = np.linspace(-1,1,int(np.sqrt(order)))
        # # generate Gray code with order n/2
        # G1 = generate_gray_code(int(n/2))
        # # reorder values with gray code
        # values = values[[int(i,base=2) for i in G1]]
        
        # # index (reordered) evenly spaced values and build complex symbols
        # const = values[L1re_idx] + 1j*values[L1im_idx]
        # # bring bits and complex symbols into correct order
        # tmp = [int(i,base=2) for i in L1]
        
        gray_I = generate_gray_code(int(n/2))
        gray_Q = generate_gray_code(int(n/2))
        
        gray_I_dec = np.array([int(i,base=2) for i in gray_I])
        gray_Q_dec = np.array([int(i,base=2) for i in gray_Q])
            
        # first half of bits denote I position, latter half of bits denote Q
        # use np.meshgrid??
        
        x_I = np.arange(int(np.sqrt(order)))
        y_Q = np.arange(int(np.sqrt(order)))
        xx,yy = np.meshgrid(x_I,y_Q)
        
        # build matrix of decimal values whose binary representations have
        # a Hamming distance of 1 in both vertical and horizontal direction
        bits = (gray_I_dec[xx]<<2) + gray_Q_dec[yy]
        # convert to binary for control purposes
        bits_bin = np.full_like(bits,0)
        for i in range(0,np.size(bits,axis=1)):
            for j in range(0,np.size(bits,axis=0)):
                bits_bin[i,j] =  np.binary_repr(bits[i,j], width=int(np.sqrt(order)))
        
        # generate evenly space values for I and Q and build matrix of complex
        # symbols
        values_I = np.linspace(-1,1,int(np.sqrt(order)))
        values_Q = np.linspace(-1,1,int(np.sqrt(order)))
        II,QQ = np.meshgrid(values_I,values_Q)
        
        symbols = II + 1j*QQ
        
        constellation = []
        tmp = []
        # bits = []
        # for i in range(0,int(np.log2(order))):
            # tmp = zip(np.where(bits == i))
            # constellation.append(symbols[zip(np.where(bits == i))])
            # bits.append(L1[tmp.index(i)])
            # constellation.append(const[L1_tmp.index(i)])
            # bits.append(L1[L1_idx.index(i)])
        
        # TODO: append symbols to constellation vector whose indices match
        # those of incrementing number in bits
        # --> bits[i,j] == 0 --> constellation.append(symbols[i,j])
        # --> bits[k,l] == 1 --> constellation.append(symbols[k,l])
        # ...
        # --> bits[x,y] == order --> constellation.append(symbols[x,y])
        for i in range(order):
            tmp.append(np.argwhere(bits == i))
        
        for i in tmp:
            constellation.append(symbols[tuple(tmp)])
    
    return constellation,bits


# =============================================================================
# generate constellation 
# TODO: fix imperfect Gray code for orders above 16!
order = 16
gray_symbols, gray_bits = generate_qam_constellation(order)

# plot constellation
fig, ax = plt.subplots()
ax.scatter(np.real(gray_symbols), np.imag(gray_symbols))
# label constellation points with their associated bit sequences
for i, txt in enumerate(gray_bits):
    ax.annotate(txt, (np.real(gray_symbols)[i], np.imag(gray_symbols)[i]))

# calculate packing coefficient C_p
C_p = np.sum(np.abs(gray_symbols)**2)/len(gray_symbols)
print('Packing coefficient C_p = {}'.format(C_p))
# =============================================================================