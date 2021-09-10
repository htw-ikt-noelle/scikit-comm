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
    constellation : list of complex
        List of complex symbols to map bits to.
    bits: list of unicode
        List of binary representations of constellation indices in unicode
        for labeling purposes.

    """
    # algorithm for generating gray-coded qam constellation alphabet
    
    # check for a reasonable order parameter
    if np.log2(order)%1:
        raise ValueError('Order must be a power of 2.')
    if type(order) != int:
        raise TypeError('Order parameter must be passed as integer.')
    # derive number of bits encoded in one symbol from QAM order
    n = int(np.log2(order))
    
    # for odd-bit QAM constellations:
    if int(np.log2(order))%2:
        # so far, rectangular QAM constellation is generated
        # TODO: implement symmetrical QAM (see paper by Smith, 1975)
        
        # generate separate Gray codes for I and Q, with [(n/2)+1] bits in I
        # and [n/2] bits in Q branch
        gray_I = generate_gray_code(int(-(-n // 2)))
        gray_Q = generate_gray_code(int(n // 2))
        # convert into indices
        gray_I_dec = np.array([int(i,base=2) for i in gray_I])
        gray_Q_dec = np.array([int(i,base=2) for i in gray_Q])
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
        # convert to binary for control purposes
        # change dtype if needed in case more than 8 bits are encoded per symbol
        bits_bin = np.full_like(bits,0,dtype='<U8')
        for i in range(0,np.size(bits,axis=0)):
            for j in range(0,np.size(bits,axis=1)):
                bits_bin[i,j] =  np.binary_repr(bits[i,j], width=int(n))
        # generate evenly spaced values for I and Q branches
        values_I = 2*np.linspace(-1,1,2**(int(-(-n // 2))))
        values_Q = np.linspace(-1,1,2**(int(n // 2)))
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
        # reassing to bits variable
        bits = bits_tmp
        
    # for even-bit QAM constellations:
    else:
        ### for even bit qam constellations
        # generate individual Gray codes for I and Q branch
        gray_I = generate_gray_code(int(n/2))
        gray_Q = generate_gray_code(int(n/2))
        # convert into decimals to allow for bit manipulation later on
        gray_I_dec = np.array([int(i,base=2) for i in gray_I])
        gray_Q_dec = np.array([int(i,base=2) for i in gray_Q])
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
        values_I = np.linspace(-1,1,int(np.sqrt(order)))
        values_Q = np.linspace(-1,1,int(np.sqrt(order)))
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
        # reassing to bits variable
        bits = bits_tmp
    
    return constellation,bits


# =============================================================================
# generate constellation 
order = 32
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