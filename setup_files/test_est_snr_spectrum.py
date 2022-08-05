import math

import numpy as np
from scipy import optimize
import scipy.interpolate as sinter
import scipy.signal as ssignal
import scipy.special as sspecial
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial

dx = 1e9
x = np.arange(193.0e12, 193.5e12+dx, dx)

# linear!!!!
y = ssignal.windows.gaussian(x.size, 30) + np.linspace(0, 0.2, x.size)

sig_range = np.asarray([193.2e12, 193.3e12])
noise_range = np.asarray([193.1e12, 193.05e12, 193.4e12, 193.45e12])

order = 5
##############################################
sig_range.sort()
sig_range = np.expand_dims(sig_range, axis=-1)
sig_range_idx = np.argmin(np.abs(x-sig_range), axis=-1)
sig_range=np.squeeze(sig_range)

noise_range.sort()

noise_range = np.expand_dims(noise_range, axis=-1)
noise_range_idx = np.argmin(np.abs(x-noise_range), axis=-1)
noise_range = np.squeeze(noise_range)

# noise_range_lft = noise_range[:2]
# noise_range_lft = np.expand_dims(noise_range_lft, axis=-1)
# noise_range_lft_idx = np.argmin(np.abs(x-noise_range_lft), axis=-1)
# noise_range_lft = np.squeeze(noise_range_lft)

# noise_range_rgt = noise_range[2:]
# noise_range_rgt = np.expand_dims(noise_range_rgt, axis=-1)
# noise_range_rgt_idx = np.argmin(np.abs(x-noise_range_rgt), axis=-1)
# noise_range_rgt = np.squeeze(noise_range_rgt)

y_sig_n2 = y[sig_range_idx[0]:sig_range_idx[1]]
dx_sig = np.diff(x[sig_range_idx[0]:sig_range_idx[1]+1])
p_sig_n2 = np.sum(y_sig_n2*dx_sig)

x_n = np.append(x[noise_range_idx[0]:noise_range_idx[1]], x[noise_range_idx[2]:noise_range_idx[3]])
y_n = np.append(y[noise_range_idx[0]:noise_range_idx[1]], y[noise_range_idx[2]:noise_range_idx[3]])
c = Polynomial.fit(x_n, y_n, order)
xx, yy = c.linspace(n=x.size, domain=[x[0], x[-1]])


plt.plot(x_n,y_n,'o')
plt.plot(xx,yy,'r-')
plt.show()




plt.plot(x,y)
plt.plot(x[sig_range_idx], y[sig_range_idx], 'o')
plt.plot(x[noise_range_idx[0]:noise_range_idx[1]], y[noise_range_idx[0]:noise_range_idx[1]], 'r--')
plt.plot(x[noise_range_idx[2]:noise_range_idx[3]], y[noise_range_idx[2]:noise_range_idx[3]], 'r--')








# p = np.poly1d([1.0, 0.0, 0.0])
# P = np.polyint(p)

# a = np.polyval(P, 2) - np.polyval(P, -2)

# print(a)
