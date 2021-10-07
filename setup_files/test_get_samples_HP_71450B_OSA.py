import sys
import os
if not any(os.path.abspath('..') == p for p in sys.path): 
    print('adding comm module to path...')
    sys.path.insert(0, os.path.abspath('.'))
    print(os.path.abspath('.'))
import numpy as np  # Version 1.20.1
import matplotlib.pyplot as plt # Version 3.3.4
import comm as comm

# Test with one trace
test_one_trace = comm.instrument_control.get_samples_HP_71450B_OSA(traces = ['A'], GPIB_address='4', log_mode=False)

# Test with two traces
test_two_traces = comm.instrument_control.get_samples_HP_71450B_OSA(traces = ['A','B'], GPIB_address='4', log_mode=False)

# Test with three traces
test_three_traces = comm.instrument_control.get_samples_HP_71450B_OSA(traces = ['A','B','C'], GPIB_address='4', log_mode=False)
