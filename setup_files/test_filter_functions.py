# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 14:59:23 2020

@author: UET-Labor
"""


import numpy as np
import comm as comm
import matplotlib.pyplot as plt
import scipy.signal as ssignal


# x = np.zeros(10)
# x[0] = 1 
# M_ave = 5

# plt.stem(x,markerfmt='C0o')
# y = comm.filters.moving_average(x, M_ave, domain='freq')
# # y = np.roll(y, -M_ave//2+1)
# plt.stem(y,markerfmt='C1o')
# plt.show()

# plt.stem(x,markerfmt='C0o')
# y = comm.filters.moving_average(x, M_ave, domain='time')
# # y = np.roll(y, -M_ave//2+1)
# plt.stem(y,markerfmt='C1o')
# plt.show()

# x = np.zeros(11)
# x[0] = 1 


# plt.stem(x,markerfmt='C0o')
# y = comm.filters.raised_cosine_filter(x, sample_rate=2.0, domain='freq', length=-1)
# plt.stem(y,markerfmt='C1o')
# plt.show()

# plt.stem(x,markerfmt='C0o')
# y = comm.filters.raised_cosine_filter(x,  sample_rate=2.0, domain='time', length=-1)
# plt.stem(y,markerfmt='C1o')
# plt.show()

# x = np.zeros(50)
# x[0] = 1 


# plt.stem(x,markerfmt='C0o')
# y = comm.filters.windowed_sinc(x, fc=0.1, order=11, window=None, domain='freq')
# plt.stem(y,markerfmt='C1o')
# plt.show()

# plt.stem(x,markerfmt='C0o')
# y = comm.filters.windowed_sinc(x, fc=0.1, order=11, window=None, domain='time')
# plt.stem(y,markerfmt='C1o')
# plt.show()


# x = np.zeros(100)
# x[0] = 1 


# plt.stem(x,markerfmt='C0o')
# y, real_fc = comm.filters.ideal_lp(x, 0.11)
# plt.stem(y,markerfmt='C1o')
# plt.show()