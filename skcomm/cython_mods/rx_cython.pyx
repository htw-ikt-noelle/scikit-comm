import numpy as np
import math
import cython
from libc.math cimport atan2
# from libcpp.complex cimport real
# import cmath

# speed gain can be increased when deactivating boundary checks for NumPy arrays
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True) # to use C divisions (e.g. for modulo operations)
def _bae_loop(double complex [:] samples_in, samples_out, double complex [:] h, 
              int n_taps, int sps, int n_CMA, double mu_cma, int n_RDE, 
              double mu_rde, double [:] radii, double mu_dde, int stop_adapting, 
              double complex [:] sig_constellation, double r, int shift, 
              bint return_info, h_tmp, eps_tmp):
    
    # define C variables
    cdef:
        # integer types
        Py_ssize_t n_out = samples_out.size
        Py_ssize_t len_constellation = len(sig_constellation)
        Py_ssize_t len_h = len(h)
        Py_ssize_t len_radii = len(radii)        
        Py_ssize_t n_update = 0
        Py_ssize_t x, sample, idx_min
        # double types
        double err, err_min, mu, r_min
        # complex types
        double complex eps, const_point, out
        # memory view of Pythons NumPy arrays        
        double complex [:] samples_out_view = samples_out
        double complex [:,:] h_tmp_view = h_tmp
        double complex [:] eps_tmp_view = eps_tmp        
    
    # equalizer loop
    # Since Cython does not know if step size is positive or negative at compile 
    # time, the standard "for-in-range loop" (for sample in range(0, n_out, shift))
    # is not converted to C!!!
    # Fallback: use "old" for loop syntax OR convert for loop into while loop
    # see https://github.com/cython/cython/issues/1106
    for sample from 0 <= sample < n_out by shift:          
        # filter the signal for each desired output sample (convolution)
        # see [1], eq. (5)
        out = 0.0
        for x in range(len_h):
            out += h[x] * samples_in[sample + n_taps - x] 
        samples_out_view[sample] = out
        
        # for each symbol, calculate error signal... 
        if (sample % sps == 0):
            # in CMA operation case
            if sample <= n_CMA:
                # calc error, see [1], eq. (26)
                eps = samples_out_view[sample] * (abs(samples_out_view[sample])**2 - r) 
                mu = mu_cma
            # in DDE operation case
            elif sample > (n_CMA + n_RDE):
                # decision (find closest point of original constellation)
                err_min = abs(samples_out_view[sample] - sig_constellation[0])
                idx_min = 0
                for x in range(1, len_constellation):
                    err = abs(samples_out_view[sample] - sig_constellation[x])
                    if err < err_min:
                        err_min = err
                        idx_min = x                                
                const_point = sig_constellation[idx_min]
                eps = (samples_out_view[sample] - const_point)
                mu = mu_dde
            # in RDE operation case
            else:
                # decision (find closest radius of original constellation)
                err_min = abs(abs(samples_out_view[sample])**2 - radii[0])
                idx_min = 0
                for x in range(len_radii):
                    err = abs(abs(samples_out_view[sample])**2 - radii[x])
                    if err < err_min:
                        err_min = err
                        idx_min = x
                r_min = radii[idx_min]    
                eps = samples_out_view[sample] * (abs(samples_out_view[sample])**2 - r_min)                         
                mu = mu_rde
            
            # ...and update impulse response, if necessary
            if (sample/sps <= stop_adapting):
                # update impulse response, see [1], eq (28)                
                for x in range(len_h):                    
                    h[x] = h[x] - mu * samples_in[sample + n_taps - x].conjugate()  * eps
            
            # save return info, if necessary
            if return_info:
                h_tmp_view[n_update,:] = h.copy()
                eps_tmp_view[n_update] = eps
            
            n_update += 1
                
    return samples_out, h_tmp, eps_tmp

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True) # to use C divisions (e.g. for modulo operations)
def _cpe_bps_loop(double complex [:] samples_norm, samples_out, dec_samples, 
                  est_phase_noise, double [:] errors, int n_blocks, int n_taps, 
                  double complex [:] rotations, 
                  double complex [:] constellation):
    
    cdef:
        Py_ssize_t block, idx_rot, idx_const, sample, min_err_rot_idx
        Py_ssize_t len_rot = len(rotations)
        Py_ssize_t len_const = len(constellation)
        
        double err, err_min, min_err_rot
        
        double complex rotated_sample
        
        # memory view of Pythons NumPy arrays
        double complex [:,:] samples_out_view = samples_out
        double [:] est_phase_noise_view = est_phase_noise
        double complex [:,:] dec_samples_view = dec_samples
    
    for block in range(n_blocks):  
        min_err_rot = 1e10
        min_err_rot_idx = 0
        for idx_rot in range(len_rot):
            errors[idx_rot] = 0.0            
            # decide nearest constellation points for each sample in block for particular test phase
            for sample in range(n_taps):
                # rotate each sample in block by test phase
                rotated_sample = samples_norm[block*n_taps+sample] * rotations[idx_rot]
                err_min = abs(rotated_sample-constellation[0])
                err_idx = 0
                for idx_const in range(1,len_const):
                    err = abs(rotated_sample-constellation[idx_const])
                    if err < err_min:
                        err_min = err
                        err_idx = idx_const
                dec_samples_view[sample, idx_rot] = constellation[err_idx]
                errors[idx_rot] += abs(rotated_sample - dec_samples_view[sample, idx_rot])**2
            if errors[idx_rot] < min_err_rot:
                min_err_rot = errors[idx_rot]
                min_err_rot_idx = idx_rot
                                    
        samples_out_view[block, :] = dec_samples_view[:,min_err_rot_idx]        
        # est_phase_noise_view[block] = np.angle(rotations[min_err_rot_idx])
        est_phase_noise_view[block] = atan2(rotations[min_err_rot_idx].imag, rotations[min_err_rot_idx].real)
        # est_phase_noise_view[block] = arg(rotations[min_err_rot_idx])
        
    
    return samples_out, est_phase_noise