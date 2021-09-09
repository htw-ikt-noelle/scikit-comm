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
        L1re = [i[:-(-len(i)//2)] for i in L1]
        L1im = [i[-(-len(i)//2):] for i in L1]
        # raise ValueError('Currently, only QAM orders equal to even powers of 2 are supported!')
    
    # for even-bit QAM constellations:
    else:
        # for even bit qam constellations
        # split Gray coded list elements (strings) in half
        L1re = [i[:len(i)//2] for i in L1]
        L1im = [i[len(i)//2:] for i in L1]
        
        # convert into indices
        L1re_idx = [int(i,base=2) for i in L1re]
        L1im_idx = [int(i,base=2) for i in L1im]
        
        # generate evenly spaced values
        values = np.linspace(-1,1,int(np.sqrt(order)))
        # generate Gray code with order n/2
        G1 = ['0','1']
        for i in range(1,int(n/2)):
            # reverse L1
            G2 = G1[::-1]
            # prefix 0 to L1 list elements
            G1 = ['0' + i for i in G1]
            # prefix 1 to L2 list elements
            G2 = ['1' + i for i in G2]
            # concatenate lists
            G1 = G1+G2
        # reorder values with gray code
        values = values[[int(i,base=2) for i in G1]]
        
        # index (reordered) evenly spaced values and build complex symbols
        const = values[L1re_idx] + 1j*values[L1im_idx]
        # bring bits and complex symbols into correct order
        tmp = [int(i,base=2) for i in L1]
            
        const_order = []
        for i in range(0,len(tmp)):
            const_order.append(const[tmp.index(i)])
    
    return const_order

res = generate_qam_constellation(16)
# print(res)
plt.scatter(np.real(res),np.imag(res))

#packing coefficient C_p
C_p = np.sum(np.abs(res)**2)/len(res)
print('Packing coefficient C_p = {}'.format(C_p))