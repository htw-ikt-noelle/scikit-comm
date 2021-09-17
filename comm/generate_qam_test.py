import numpy as np
import matplotlib.pyplot as plt

def generate_qam_constellation(order):
    """
    Generates a Gray-coded [order]-QAM constellation.

    Parameters
    ----------
    order : int
        number of symbols in the constellation.
        
    Raises
    ------
    ValueError
        If 'order' is not a power of 2.
    TypeError
        If 'order' is not an integer.

    Returns
    -------
    constellation : array of complex
        Array of complex symbols to map bits to.
    bits: array of str
        Array of binary representations of constellation indices as strings
        for labeling purposes.

    """
    # check for a reasonable order parameter
    if np.log2(order)%1:
        raise ValueError('Order must be a power of 2.')
    if type(order) != int:
        raise TypeError('Order parameter must be passed as integer.')
    # derive number of bits encoded in one symbol from QAM order
    n = int(np.log2(order))
    
    # for odd-bit QAM constellations:
    if int(np.log2(order))%2:
        ### build rectangular matrix of gray-coded bit words
        # generate separate Gray codes for I and Q using XOR with [(n/2)+1]
        # bits in I and [n/2] bits in Q branch
        width_bit = int(-(-n // 2))
        height_bit = int(n // 2)
        gray_I_dec = np.arange(2**width_bit) ^ ((np.arange(2**width_bit))>>1)
        gray_Q_dec = np.arange(2**height_bit) ^ ((np.arange(2**height_bit))>>1)
        # generate indices for I and Q values to build matrix to project 
        # constellation points onto later
        x_I = np.arange(int(2**(-(-n // 2))))
        y_Q = np.arange(int(2**(n // 2)))
        # combine into meshgrid
        xx,yy = np.meshgrid(x_I,y_Q)
        # build matrix of decimal values whose binary representations have
        # a Hamming distance of 1 in both vertical and horizontal direction by
        # shifting bits of I-branch Gray code left by (n/2) and adding bits of 
        # Q-branch Gray code
        bits = (gray_I_dec[xx]<<int(n // 2)) + gray_Q_dec[yy]

        ### build symmetrical matrix of complex symbols
        # calc dimension for symmetrical constellation matrix
        cross_dim = int((2**(-(- n // 2)) + 2**(n // 2)) / 2)
        # generate evenly spaced values with euclidean distance of 2
        cross_values = np.linspace(-cross_dim + 1, cross_dim - 1, cross_dim)
        # generate meshgrid
        cross_I, cross_Q = np.meshgrid(cross_values,cross_values)
        # build complex symbols
        cross_symbols = cross_I + 1j*cross_Q
        # cut away corners
        cut = int(cross_dim / 6)
        cross_symbols[:cut,:cut] = 0
        cross_symbols[-cut:,:cut] = 0
        cross_symbols[:cut,-cut:] = 0
        cross_symbols[-cut:,-cut:] = 0
        # copy matrix for assigning gray-coded decimal values to
        bits_symm = np.full_like(cross_symbols,0,dtype='int')
        # write 'middle block' of rectangular constellation into symmetrical
        # matrix
        bits_symm[cut:-cut,:] = bits[:,cut:-cut]
        # manipulate the 8 'end blocks' of rectangular constellation and 
        # write them into new positions in the symmetrical matrix
        # top left block
        bits_symm[:cut,cut:2*cut] = np.flipud(bits[:cut,:cut])
        # upper middle left block
        bits_symm[:cut,2*cut:3*cut] = np.fliplr(bits[cut:2*cut,:cut])
        # lower middle left block
        bits_symm[-cut:,2*cut:3*cut] = np.fliplr(bits[-(2*cut):-cut,:cut])
        # bottom left block
        bits_symm[-cut:,cut:2*cut] = np.flipud(bits[-cut:,:cut])
        # top right block
        bits_symm[:cut,-(2*cut):-cut] = np.flipud(bits[:cut,-cut:])
        # upper middle right block
        bits_symm[:cut,-(3*cut):-(2*cut)] = np.fliplr(bits[cut:2*cut,-cut:])
        # lower middle right block
        bits_symm[-cut:,-(3*cut):-(2*cut)] = np.fliplr(bits[-(2*cut):-cut,-cut:])
        # bottom right block
        bits_symm[-cut:,-(2*cut):-cut] = np.flipud(bits[-cut:,-cut:])
        
        ### manipulate and reshape symmetrical matrix into array of connstellation points
        # flatten matrices out and delete entries at indices where
        # cross_symbols == 0 (the corners that were cut away)
        bits_symm = np.delete(bits_symm.flatten(),np.argwhere(cross_symbols.flatten()==0))
        bits_bin = np.full_like(bits_symm,0,dtype='<U8')
        for i in range(len(bits_symm)):
            bits_bin[i] = np.binary_repr(bits_symm[i],width=int(n))
        cross_symbols = np.delete(cross_symbols.flatten(),np.argwhere(cross_symbols.flatten()==0))
        # initialize lists for return values
        constellation = []
        bits_tmp = []
        # iterate over flattened symbols and bits matrices and append complex
        # constellation points and binary number labels to respective lists
        for i in range(order):
            constellation.append(cross_symbols[np.argwhere(bits_symm == i)][0][0])
            bits_tmp.append(bits_bin[np.argwhere(bits_symm == i)][0][0])
        # convert into arrays
        bits = np.asarray(bits_tmp)
        constellation = np.asarray(constellation)
        
    # for even-bit QAM constellations:
    else:
        ### for even bit qam constellations
        # generate individual Gray codes for I and Q branch
        width_bit = int(-(-n // 2))
        height_bit = int(n // 2)
        gray_I_dec = np.arange(2**width_bit) ^ ((np.arange(2**width_bit))>>1)
        gray_Q_dec = np.arange(2**height_bit) ^ ((np.arange(2**height_bit))>>1)
        # generate indices for I and Q values to build matrix to project 
        # constellation points onto later
        x_I = np.arange(int(np.sqrt(order)))
        y_Q = np.arange(int(np.sqrt(order)))
        # combine into meshgrid
        xx,yy = np.meshgrid(x_I,y_Q)
        # build matrix of decimal values whose binary representations have
        # a Hamming distance of 1 in both vertical and horizontal direction by
        # shifting bits of I-branch Gray code left by (n/2) and adding bits of 
        # Q-branch Gray code
        bits = (gray_I_dec[xx]<<int(n/2)) + gray_Q_dec[yy]
        # convert to binary for control purposes
        # change dtype if needed in case more than 8 bits are encoded per symbol
        bits_bin = np.full_like(bits,0,dtype='<U8')
        for i in range(0,np.size(bits,axis=1)):
            for j in range(0,np.size(bits,axis=0)):
                bits_bin[i,j] =  np.binary_repr(bits[i,j], width=int(n))
        # generate evenly space values for I and Q and build matrix of complex
        # symbols
        values_I = np.linspace(-(np.sqrt(order))+1,np.sqrt(order)-1,int(np.sqrt(order)))
        values_Q = np.linspace(-(np.sqrt(order))+1,np.sqrt(order)-1,int(np.sqrt(order)))
        II,QQ = np.meshgrid(values_I,values_Q)
        symbols = II + 1j*QQ
        # initialize lists for return values
        constellation = []
        bits_tmp = []
        # iterate over flattened symbols and bits matrices and append complex
        # constellation points and binary number labels to respective lists
        for i in range(order):
            constellation.append(symbols.flatten()[np.argwhere(bits.flatten() == i)][0][0])
            bits_tmp.append(bits_bin.flatten()[np.argwhere(bits.flatten() == i)][0][0])
        # convert into arrays
        bits = np.asarray(bits_tmp)
        constellation = np.asarray(constellation)
    
    return constellation,bits

def generate_psk_constellation(order):
    """
    Generates a Gray coded PSK constellation of the specified order

    Parameters
    ----------
    order : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # check for a reasonable order parameter
    if np.log2(order)%1:
        raise ValueError('Order must be a power of 2.')
    if type(order) != int:
        raise TypeError('Order parameter must be passed as integer.')
    # derive number of bits encoded in one symbol from QAM order
    n = int(np.log2(order))
    # generate Gray code
    gray = np.arange(2**n) ^ (np.arange(2**n)>>1)
    gray_bin = np.full_like(gray,0,dtype='<U8')
    for i in range(len(gray)):
        gray_bin[i] = np.binary_repr(gray[i],width=int(n))
    # build constellation points
    symbols = np.asarray([np.exp(1j*2*np.pi*i/order) for i in range(len(gray))]) 
    # reorder symbols and label vector
    constellation = []
    bits = []
    for i in range(order):
        constellation.append(symbols.flatten()[np.argwhere(gray.flatten()==i)][0][0])
        bits.append(gray_bin.flatten()[np.argwhere(gray.flatten()==i)][0][0])
    return constellation, bits

def Smith_bit_manip(value,nbits):
    """
    Performs bit manipulation on an input integer with odd number of bits
    according to the paper by Smith.

    Parameters
    ----------
    value : int
        Integer value to potentially have its bits manipulated.
    nbits : int
        Number of bits encoded per symbol.

    Returns
    -------
    None.

    """
    # check for reasonable input parameters
    if type(nbits) != int or type(value) != int:
        raise TypeError('Parameters must be passed as integer.')
    if not nbits%2 or nbits<5:
        raise ValueError('nbits must be an odd integer >= 5')
    if value < 0 or value > 2**nbits:
        raise ValueError('Value must be within the interval [0, 2**nbits]')
    # add an extra zero bit to input word
    value = value << 1
    # obtain floored half of input bit length
    n = nbits // 2
    # model: two half-words with n+1 bits each (for total length of 2n+2 bits)
    x_n_complement = ((value >> (n+2)) ^ 1) & 1
    x_nPlus1 = (value >> (n+1)) & 1
    y_n = (value >> 1) & 1
    y_n_complement = ((value >> 1) ^ 1) & 1
    # if X_n* and X_(n+1) =/= 1
    if not (x_n_complement & x_nPlus1):
        print('no change performed')
    # if X_n* and X_(n+1) and Y_n = 1
    elif x_n_complement & x_nPlus1 & y_n:
        print('complement x_n, x_nPlus1, y_nPlus1')
        value = value ^ ((1<<(n+2)) + (1<<(n+1)) + 1)
    # if X_n and X_(n+1) and Y_n* = 1
    elif x_n_complement & x_nPlus1 & y_n_complement:
        print('complement x_nPlus1, y_n, y_nPlus1')
        value = value ^ ((1<<(n+1)) + (1<<1) + 1)
    # if none of these conditions apply (shouldn't happen!)
    else:
        raise ValueError('Unexpected bit sequence encountered.')
    return value
# =============================================================================
### QAM constellation 
order = 32
gray_symbols, gray_bits = generate_qam_constellation(order)

# plot constellation
fig, ax = plt.subplots()
ax.scatter(np.real(gray_symbols), np.imag(gray_symbols))
# label constellation points with their associated bit sequences
for i, txt in enumerate(gray_bits):
    ax.annotate(txt, (np.real(gray_symbols)[i], np.imag(gray_symbols)[i]))

# calculate packing coefficient C_p
# C_p = np.sum(np.abs(gray_symbols)**2)/len(gray_symbols)
# print('Packing coefficient C_p = {}'.format(C_p))

### PSK constellation
psk_symbols, psk_bits = generate_psk_constellation(order)
# plot constellation
fig, ax = plt.subplots()
ax.scatter(np.real(psk_symbols), np.imag(psk_symbols))
# label constellation points with their associated bit sequences
for i, txt in enumerate(psk_bits):
    ax.annotate(txt, (np.real(psk_symbols)[i], np.imag(psk_symbols)[i]))

### Bit manipulation test    
a = 23
abin = bin(a)
b = Smith_bit_manip(a,5)
print(b)
# =============================================================================